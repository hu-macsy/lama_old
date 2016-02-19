/**
 * @file GPICommunicator.cpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief GPICommunicator.cpp
 * @author Thomas Brandes, Lauretta Schubert
 * @date 06.05.2014
 * @since 1.1.0
 */

// hpp
#include <scai/dmemo/gpi/GPICommunicator.hpp>

// others
#include <scai/dmemo/gpi/GPISyncToken.hpp>
#include <scai/dmemo/gpi/GPIMemManager.hpp>
#include <scai/dmemo/gpi/GPIMemory.hpp>
#include <scai/tasking/NoSyncToken.hpp>
#include <scai/dmemo/gpi/GPIUtils.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/bind.hpp>

// tracing
#include <scai/tracing.hpp>

// logging
#include <scai/logging.hpp>

#include <iostream>
#include <algorithm>
#include <unistd.h>

#include <GASPI.h>

using namespace std;

namespace scai
{

using common::TypeTraits;
using common::shared_ptr;
using common::unique_ptr;

namespace dmemo
{

/* ---------------------------------------------------------------------------------- */

const int GPICommunicator::defaultTag = 1;

const gaspi_timeout_t GPICommunicator::timeout = GASPI_BLOCK;
const gaspi_group_t   GPICommunicator::group   = GASPI_GROUP_ALL;

SCAI_LOG_DEF_LOGGER( GPICommunicator::logger, "Communicator.GPICommunicator" )

/* ---------------------------------------------------------------------------------- */

static shared_ptr<const GPICommunicator> theGPICommunicator;

CommunicatorPtr GPICommunicator::create()
{
    if ( !theGPICommunicator )
    {
        SCAI_LOG_INFO( logger, "create new GPICommunicator" )

        // create a new instance of GPICommunicator, will call GPI_Init

        theGPICommunicator = shared_ptr<const GPICommunicator>( new GPICommunicator() );
    }

    return theGPICommunicator;
}

communicator::CommunicatorKind GPICommunicator::createValue()
{
    return communicator::GPI;
}

/* ---------------------------------------------------------------------------------- */

GPICommunicator::GPIGuard::GPIGuard()
{
}

GPICommunicator::GPIGuard::~GPIGuard()
{
    if ( theGPICommunicator.use_count() > 1 )
    {
        SCAI_LOG_WARN( logger,
                       "GPICommunicator has " << theGPICommunicator.use_count() - 1 
                        << " remaining references, seems that not all distributed data structures have been freed" )
    }
}

GPICommunicator::GPIGuard GPICommunicator::guard;

/* ---------------------------------------------------------------------------------- */

GPICommunicator::GPICommunicator( )
    : CRTPCommunicator<GPICommunicator>( communicator::GPI ), 
      mQueueID( 0 ), 
      mThreadSafetyLevel( Communicator::Funneled )
{
    SCAI_TRACE_SCOPE( false )   // switch off tracing in this scope as it might call this constructor again

    // set gaspi_printf as print routine for logging

    logging::GenLogger::myPrintf = &gaspi_printf;

    SCAI_LOG_DEBUG( logger, "GPICommunicator(): call init" )

    gaspi_timeout_t init_timeout = 10000;   // in milli seconds

    SCAI_GASPI_CALL( gaspi_proc_init( init_timeout ) )

    SCAI_GASPI_CALL( gaspi_proc_num( &mSize ) )
    SCAI_GASPI_CALL( gaspi_proc_rank( &mRank ) )

    SCAI_LOG_INFO( logger, "GASPI proc " << mRank << " of " << mSize << " started")

    gaspi_number_t num_notifications;

    SCAI_GASPI_CALL( gaspi_notification_num ( &num_notifications ) )

    mMaxNotifications = num_notifications;

    SCAI_ASSERT_ERROR( mMaxNotifications >= 2 * mSize, 
                       "# notifications = " << mMaxNotifications << " not sufficient" 
                       << ", need at least 2 * " << mSize )

    // the following call is for convenience to get correct node rank by hostnames

    setNodeData();   // determine mNodeRank, mNodeSize, also uses some communication pattern
}

/* ---------------------------------------------------------------------------------- */

#define GPI_MAX_PROCESSOR_NAME 512

void GPICommunicator::setNodeData()
{
    char nodeName[GPI_MAX_PROCESSOR_NAME];  // name of node for this processor

    memset( nodeName,'\0', GPI_MAX_PROCESSOR_NAME );

    // GPI is available only on Linux, so use this Linux routine

    gethostname( nodeName, GPI_MAX_PROCESSOR_NAME );

    SCAI_LOG_DEBUG( logger, "Processor " << mRank << " runs on node " << nodeName )

    char* allNodeNames = new char[ GPI_MAX_PROCESSOR_NAME * mSize ];

    if( allNodeNames == NULL )
    {
        COMMON_THROWEXCEPTION( "Can't alloc enough memory for all node names." )
    }

    memset( allNodeNames, '\0', GPI_MAX_PROCESSOR_NAME * mSize * sizeof( char) );

    const PartitionId root = 0;

    gatherImpl( allNodeNames, GPI_MAX_PROCESSOR_NAME, root, nodeName );

    bcastImpl( allNodeNames, GPI_MAX_PROCESSOR_NAME * mSize, root );

    mNodeSize = 0;
    mNodeRank = mSize;  // illegal value to verify that it will be set

    for ( int i = 0; i < mSize; ++i )
    {
        if ( strcmp( &allNodeNames[ i * GPI_MAX_PROCESSOR_NAME ], nodeName ) )
        {
            continue;   // processor i is not on same node
        }

        // Processor i is on same node as this processor

        if ( i == mRank )
        {
            mNodeRank = mNodeSize;
        }

        ++mNodeSize;
    }

    delete[] allNodeNames;

    SCAI_ASSERT_ERROR( mNodeSize > 0, "Serious problem encountered to get node size" )

    SCAI_ASSERT_ERROR( mNodeRank < mNodeSize, "Serious problem encountered to get node size" )

    SCAI_LOG_INFO( logger, *this << ": runs on " << nodeName 
                           << ", node rank " << mNodeRank << " of " << mNodeSize )
}

/* ---------------------------------------------------------------------------------- */

GPICommunicator::~GPICommunicator()
{
    SCAI_LOG_INFO( logger, *this << ": ~GPICommunicator" )

    GPIMemManager::freeAll();

    // do not call synchronize here as instrumentation of that routine causes problems

    SCAI_GASPI_CALL( gaspi_barrier( GASPI_GROUP_ALL, timeout ) )
    SCAI_GASPI_CALL( gaspi_proc_term( timeout ) )
}

/* ---------------------------------------------------------------------------------- */

bool GPICommunicator::isEqual( const Communicator& other ) const
{
    bool equal = false;
    const GPICommunicator* otherGPI = dynamic_cast<const GPICommunicator*>( &other );

    if ( otherGPI )
    {
        equal = mQueueID == otherGPI->mQueueID;
    }

    return equal;
}

/* ---------------------------------------------------------------------------------- */

Communicator::ThreadSafetyLevel GPICommunicator::getThreadSafetyLevel() const
{
    return mThreadSafetyLevel;
}

PartitionId GPICommunicator::getSize() const
{
    return mSize;
}

PartitionId GPICommunicator::getRank() const
{
    return mRank;
}

PartitionId GPICommunicator::getNodeSize() const
{
    return mNodeSize;
}

PartitionId GPICommunicator::getNodeRank() const
{
    return mNodeRank;
}

/* ---------------------------------------------------------------------------------- */
/*              getGPIType                                                            */
/* ---------------------------------------------------------------------------------- */

template<>
inline gaspi_datatype_t GPICommunicator::getGPIType<float>()
{
    return GASPI_TYPE_FLOAT;
}

template<>
inline gaspi_datatype_t GPICommunicator::getGPIType<double>()
{
    return GASPI_TYPE_DOUBLE;
}

template<>
inline gaspi_datatype_t GPICommunicator::getGPIType<int>()
{
    return GASPI_TYPE_INT;
}

template<>
inline gaspi_datatype_t GPICommunicator::getGPIType<unsigned long>()
{
    return GASPI_TYPE_ULONG;
}

template<>
inline gaspi_datatype_t GPICommunicator::getGPIType<LongDouble>()
{
    COMMON_THROWEXCEPTION( "LongDouble communication not suppported by GASPI" )
    return GASPI_TYPE_ULONG;
}

template<>
inline gaspi_datatype_t GPICommunicator::getGPIType<ComplexFloat>()
{
    COMMON_THROWEXCEPTION( "ComplexFloat communication not suppported by GASPI" )
    return GASPI_TYPE_FLOAT;
}

template<>
inline gaspi_datatype_t GPICommunicator::getGPIType<ComplexDouble>()
{
    COMMON_THROWEXCEPTION( "ComplexDouble communication not suppported by GASPI" )
    return GASPI_TYPE_DOUBLE;
}

template<>
inline gaspi_datatype_t GPICommunicator::getGPIType<ComplexLongDouble>()
{
    COMMON_THROWEXCEPTION( "ComplexLongDouble communication not suppported by GASPI" )
    return GASPI_TYPE_DOUBLE;
}

void GPICommunicator::wait() const
{
    // SCAI_REGION( "Communicator.GPI.wait" )
    SCAI_GASPI_CALL ( gaspi_wait ( mQueueID, timeout ) )
}

/* ---------------------------------------------------------------------------------- */
/*                                     all to all                                     */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::all2all( IndexType recvSizes[], const IndexType sendSizes[] ) const
{
    SCAI_REGION( "Communicator.GPI.all2all" )

    SCAI_ASSERT_ERROR( sendSizes != 0, " invalid sendSizes " )
    SCAI_ASSERT_ERROR( recvSizes != 0, " invalid recvSizes " )


    // mSize will be the same on all processors

    SegmentData<int> sendSegment( this, mSize, const_cast<IndexType*>( sendSizes ) );  // locally used for sending
    SegmentData<int> recvSegment( this, mSize, recvSizes );  // other procs do remote write here

    SCAI_LOG_DEBUG( logger, *this << ": all2all exchange of sizes" )

    // Step 1: notify each receive partner that segment is available and give him the write offset

    for ( PartitionId fromP = 0; fromP < mSize; ++fromP )
    {
        if ( fromP == mRank )
        {
            continue;
        }

        const IndexType recvOffset = fromP + recvSegment.getOffset();

        notify( recvSegment.getID(), fromP, mRank, recvOffset );
    }

    sendSegment.assign( sendSizes, mSize );

    // Step 2: wait for notifications and then write the remote data

    for ( PartitionId toP = 0; toP < mSize; ++toP )
    {
        // blocking wait until we get the remote offset that also indicates availability

        const IndexType quantity = 1;
        const IndexType localOffset = toP;

        if ( toP != mRank )
        {
            const IndexType remoteOffset = notifyWait( recvSegment.getID(), toP );

            remoteWrite( sendSegment, localOffset, toP, recvSegment, remoteOffset, quantity );

            SCAI_LOG_DEBUG( logger, *this << ": write " << toP << " " << sendSegment[ localOffset ] )

            // notify partner and tell him the number of written items

            notify( recvSegment.getID(), toP, mRank + mSize, 1 );
        }
        else
        {
            // recvSegment[mRank] = sendSegment[mRank];
            localWrite( sendSegment, localOffset, recvSegment, mRank, quantity );
        }
    }

    wait();  // wait that all writes of sendSegment have been done

    // Step 3: wait until my receive partners have me written my needed values

    for ( PartitionId fromP = 0; fromP < mSize; ++fromP )
    {
        if ( fromP == mRank )
        {
            continue;
        }

        const IndexType quantity = notifyWait( recvSegment.getID(), fromP + mSize );

        // notifiaction value is written size that we match against what is expected

        SCAI_ASSERT_EQUAL_ERROR( 1, quantity );
    }

    recvSegment.copyTo( recvSizes, mSize );  // copy from recvSegment back to recvSizes

    SCAI_LOG_DEBUG( logger, *this << ": all2all done" )
}

/* ---------------------------------------------------------------------------------- */


template<typename ValueType>
void GPICommunicator::all2allvImpl( 
    ValueType* recvBuffer[], 
    IndexType recvCount[], 
    ValueType* sendBuffer[], 
    IndexType sendCount[] ) const
{
    SCAI_ASSERT_UNEQUAL( recvBuffer, 0, "illegal" )
    SCAI_ASSERT_UNEQUAL( recvCount, 0, "illegal" )
    SCAI_ASSERT_UNEQUAL( sendBuffer, 0, "illegal" )
    SCAI_ASSERT_UNEQUAL( sendCount, 0, "illegal" )

    COMMON_THROWEXCEPTION( "Not available yet" )
}

/* ---------------------------------------------------------------------------------- */
/*                         notify                                                     */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::notify( const gaspi_segment_id_t segID, 
                              PartitionId target, 
                              PartitionId pos, 
                              IndexType val ) const
{
    SCAI_LOG_DEBUG( logger, *this << ": notify processor " << target << " at id " << pos << " with val = " << val )

    SCAI_ASSERT_ERROR( pos < mMaxNotifications, "notifiation id " << pos << " illegal, max = " << mMaxNotifications )
    SCAI_ASSERT_ERROR( target < mSize, "target rank = " << target << " illegal, size = " << mSize )
    SCAI_ASSERT_ERROR( target != mRank, *this << ": notify for myself is not supported yet" )
    SCAI_ASSERT_ERROR( val >= 0, "notify val must not be negative" )

    // Be careful: val must not be 0

    gaspi_notification_t notVal = static_cast<gaspi_notification_t>( val + 1 );

    gaspi_rank_t            rank   = target;
    gaspi_notification_id_t notID  = pos;

    SCAI_GASPI_CALL( gaspi_notify( segID, rank, notID, notVal, mQueueID, GASPI_BLOCK ) );
}

/* ---------------------------------------------------------------------------------- */
/*                         notifyWait                                                 */
/* ---------------------------------------------------------------------------------- */

IndexType GPICommunicator::notifyWait( const gaspi_segment_id_t segID, PartitionId pos ) const
{
    SCAI_LOG_DEBUG( logger, *this << ": wait for notification at " << pos );

    SCAI_ASSERT_ERROR( pos < mMaxNotifications, "notifiation id " << pos << " illegal, max = " << mMaxNotifications )

    gaspi_notification_id_t begin_id = pos;

    gaspi_notification_id_t notID;

    SCAI_GASPI_CALL ( gaspi_notify_waitsome( segID, begin_id, 1 , &notID, GASPI_BLOCK) );

    gaspi_notification_t notifyVal;

    SCAI_GASPI_CALL( gaspi_notify_reset( segID, notID, &notifyVal) );

    SCAI_LOG_DEBUG( logger, *this << ": got notification at " << pos << ", val = " << ( notifyVal - 1 ) );

    return static_cast<IndexType>( notifyVal - 1 );
}

/* ---------------------------------------------------------------------------------- */
/*                                exchange by plan                                    */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void GPICommunicator::exchangeByPlanImpl(
    T recvData[],
    const CommunicationPlan& recvPlan,
    const T sendData[],
    const CommunicationPlan& sendPlan ) const
{
    SCAI_REGION( "Communicator.GPI.exchangeByPlan" )

    SCAI_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    SCAI_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )

    SCAI_LOG_DEBUG( logger, *this << ": exchange by plan for values of type " << TypeTraits<T>::id()
                             << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )

    int noReceives = recvPlan.size();
    int noSends = sendPlan.size();

    int totalSendSize = sendPlan.totalQuantity();
    int totalRecvSize = recvPlan.totalQuantity();

    // setup segment data

    SegmentData<T> srcDataSegment( this, totalSendSize, const_cast<T*>( sendData ) );
    SegmentData<T> dstDataSegment( this, totalRecvSize, recvData );

    SCAI_LOG_DEBUG( logger, *this << ": exchangeByPlan, send = " << noSends << ", " << totalSendSize 
                              << ", recv = " << noReceives << ", " << totalRecvSize )

    // Step 1: notify each receive partner that segment is available and give him the write offset

    for ( PartitionId i = 0; i < noReceives; ++i )
    {
        const PartitionId source = recvPlan[i].partitionId;
        const IndexType   offset = dstDataSegment.getOffset() + recvPlan[i].offset;
        notify( dstDataSegment.getID(), source, mRank, offset );
    }

    // now fill my buffer with the send data

    srcDataSegment.assign( sendData, totalSendSize );  // copy only if new segment data was allocated

    // Step 2: wait for notifications and then write the remote data

    for ( PartitionId i = 0; i < noSends; ++i )
    {
        IndexType quantity = sendPlan[i].quantity;
        PartitionId toP    = sendPlan[i].partitionId;

        // blocking wait until we get the remote offset that also indicates availability

        const IndexType remoteOffset = notifyWait( dstDataSegment.getID(), toP );
        const IndexType localOffset  = sendPlan[i].offset;

        remoteWrite( srcDataSegment, localOffset, toP, dstDataSegment, remoteOffset, quantity );

        // notify partner and tell him the number of written items

        notify( dstDataSegment.getID(), toP, mRank + mSize, quantity );
    }

    // Step 3: wait until my receive partners have me written my needed values

    for ( PartitionId i = 0; i < noReceives; ++i )
    {
        const PartitionId fromP    = recvPlan[i].partitionId;
        const IndexType   quantity = recvPlan[i].quantity;

        const IndexType val = notifyWait( dstDataSegment.getID(), fromP + mSize );

        // notifiaction value is written size that we match against what is expected

        SCAI_ASSERT_EQUAL_ERROR( quantity, val );
    }

    wait();

    // cpy for dst vector value

    dstDataSegment.copyTo( recvData, totalRecvSize ); // copy from segment to ptr

    // synchronize();   // make sure that my data is available until all procs have read it

    SCAI_LOG_DEBUG( logger, *this << ": free segments" )
}

/* ---------------------------------------------------------------------------------- */
/*                                exchange by plan async                              */
/* ---------------------------------------------------------------------------------- */

template<typename T>
tasking::SyncToken* GPICommunicator::exchangeByPlanAsyncImpl(
    T recvData[],
    const CommunicationPlan& recvPlan,
    const T sendData[],
    const CommunicationPlan& sendPlan ) const
{
    SCAI_REGION( "Communicator.GPI.exchangeByPlanAsync" )

    SCAI_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    SCAI_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )

    // Some global information from the communication plans

    int noReceives = recvPlan.size();
    int noSends = sendPlan.size();

    int totalSendSize = sendPlan.totalQuantity();
    int totalRecvSize = recvPlan.totalQuantity();

    SCAI_LOG_INFO( logger, *this << ": async exchange for values of type " << TypeTraits<T>::id() 
                           << ", send to " << sendPlan.size() << " processors " << totalSendSize << " values"
                           << ", recv from " << recvPlan.size() << " processors " << totalRecvSize << " values" )

    // setup segment data

    shared_ptr<SegmentData<T> > srcDataSegment( new SegmentData<T>( this, totalSendSize ) );
    shared_ptr<SegmentData<T> > dstDataSegment( new SegmentData<T>( this, totalRecvSize ) );

    {
        // Step 1: notify each receive partner that segment is available and give him the write offset

        const IndexType segOffset = dstDataSegment->getOffset();
        const gaspi_segment_id_t segId = dstDataSegment->getID();

        for ( PartitionId i = 0; i < noReceives; ++i )
        {
            const PartitionId source = recvPlan[i].partitionId;
            const IndexType   offset = segOffset + recvPlan[i].offset;
            notify( segId, source, mRank, offset );
        }
    }

    srcDataSegment->assign( sendData, totalSendSize );

    // Step 2: wait for notifications and then write the remote data

    for ( PartitionId i = 0; i < noSends; ++i )
    {
        IndexType quantity = sendPlan[i].quantity;
        PartitionId toP    = sendPlan[i].partitionId;

        // blocking wait until we get the remote offset that also indicates availability

        const IndexType remoteOffset = notifyWait( dstDataSegment->getID(), toP );
        const IndexType localOffset  = sendPlan[i].offset;

        remoteWrite( *srcDataSegment, localOffset, toP, *dstDataSegment, remoteOffset, quantity );

        // notify partner and tell him the number of written items

        notify( dstDataSegment->getID(), toP, mRank + mSize, quantity );
    }

    // need an GPI Sync Token that waits on queue

    unique_ptr<GPISyncToken> pSyncToken( new GPISyncToken( dstDataSegment->getID(), noReceives ) );

    // Step 3: deal with remote write notifications

    for ( PartitionId i = 0; i < noReceives; ++i )
    {
        const PartitionId fromP = recvPlan[i].partitionId;
        pSyncToken->pushNotification( fromP + mSize );
    }

    // Keep SegmentData until synchronization of token

    pSyncToken->pushToken( srcDataSegment );
    pSyncToken->pushToken( dstDataSegment );

    // routines to be executed after the asnchronous write
    // To be careful: routines should be executed before accesses are freed

    pSyncToken->pushRoutine( common::bind( &SegmentData<T>::copyTo, dstDataSegment.get(), recvData, totalRecvSize ) );
    // Also wait that all remote writes have been finished
    pSyncToken->pushRoutine( common::bind( &GPICommunicator::wait, this ) );

    return pSyncToken.release();
}

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

template<typename T>
IndexType GPICommunicator::shiftImpl(
    T recvVals[],
    const IndexType recvSize,
    const PartitionId source,
    const T sendVals[],
    const IndexType sendSize,
    const PartitionId dest ) const
{
    SCAI_REGION( "Communicator.GPI.shift" )

    SCAI_ASSERT_ERROR( source != mRank, "source must not be this partition" )
    SCAI_ASSERT_ERROR( dest != mRank, "dest must not be this partition" )

    SCAI_ASSERT_ERROR( sendSize <= recvSize, "too large send" )

    SCAI_LOG_INFO( logger,
                    *this << ": shift<" << TypeTraits<T>::id() 
                          << ">, recv from " << source << " max " << recvSize << " values " 
                          << ", send to " << dest << " " << sendSize << " values." )

    SegmentData<T> srcDataSegment( this, recvSize, const_cast<T*>( sendVals ) );
    SegmentData<T> dstDataSegment( this, recvSize, recvVals );

    SCAI_LOG_DEBUG( logger, "shift: srcDataSegment = " << srcDataSegment )
    SCAI_LOG_DEBUG( logger, "shift: dstDataSegment = " << dstDataSegment )

    // Give my source process the offset where it can write the data
    {
        const IndexType offset = dstDataSegment.getOffset();
        notify( dstDataSegment.getID(), source, mRank, offset );
    }

    {
        // SCAI_REGION( "Communicator.GPI.shift.assign" )
        srcDataSegment.assign( sendVals, sendSize );  // skipped if sendVals is segment data
    }

    IndexType realRecvSize;

    {
        // SCAI_REGION( "Communicator.GPI.shift.remWrite" )

        const IndexType remoteOffset = notifyWait( dstDataSegment.getID(), dest );
        const IndexType localOffset  = 0;
        remoteWrite( srcDataSegment, localOffset, dest, dstDataSegment, remoteOffset, sendSize );
    }
    {
        // SCAI_REGION( "Communicator.GPI.shift.notification" )

        notify( dstDataSegment.getID(), dest, mRank + mSize, sendSize );
        realRecvSize = notifyWait( dstDataSegment.getID(), source + mSize );
    }

    // cpy the real receive size and the received data

    {
        // SCAI_REGION( "Communicator.GPI.shift.copyOut" )

        dstDataSegment.copyTo( recvVals, recvSize );  // skipped if recvVals is segment data
    }

    wait();  // Queues should be empty before returning

    SCAI_LOG_DEBUG( logger, "shift ready, recvSize = " << realRecvSize )

    return realRecvSize;
}

/* ---------------------------------------------------------------------------------- */
/*              shiftAsync                                                            */
/* ---------------------------------------------------------------------------------- */

template<typename T>
tasking::SyncToken* GPICommunicator::shiftAsyncImpl(
    T recvVals[],
    const PartitionId source,
    const T sendVals[],
    const PartitionId dest,
    const IndexType size ) const
{
    // SCAI_REGION( "Communicator.GPI.shiftAsyncImpl" )

    SCAI_LOG_DEBUG( logger,
                    *this << ": recv from " << source << ", send to " << dest << ", both " << size << " values." )

    SCAI_ASSERT_ERROR( source != mRank, "source must not be this partition" )
    SCAI_ASSERT_ERROR( dest != mRank, "dest must not be this partition" )

    shared_ptr<SegmentData<T> > srcSegment( new SegmentData<T>( this, size ) );
    shared_ptr<SegmentData<T> > dstSegment( new SegmentData<T>( this, size ) );

    {
        const IndexType offset = dstSegment->getOffset();
        notify( dstSegment->getID(), source, mRank, offset );
    }

    srcSegment->assign( sendVals, size );

    // need an GPI communicator with 1 notification

    unique_ptr<GPISyncToken> pSyncToken( new GPISyncToken( dstSegment->getID(), 1 ) );

    const IndexType remoteOffset = notifyWait( dstSegment->getID(), dest );
    const IndexType localOffset  = 0;
    remoteWrite( *srcSegment, localOffset, dest, *dstSegment, remoteOffset, size );
    notify( dstSegment->getID(), dest, mRank + mSize, size );

    // wait on receive is pushed in SyncToken

    pSyncToken->pushNotification( source + mSize );

    // routines to be executed after the wait

    pSyncToken->pushRoutine( common::bind( &SegmentData<T>::copyTo, dstSegment.get(), recvVals, size ) );

    // Keep SegmentData until synchronization of token

    pSyncToken->pushToken( srcSegment );
    pSyncToken->pushToken( dstSegment );

    return pSyncToken.release();
}

/* ---------------------------------------------------------------------------------- */
/*              sum                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename T>
T GPICommunicator::sumImpl( const T value ) const
{
    SCAI_REGION( "Communicator.GPI.sum" )

    T localSum = value;   // might be used also as buffer

    T globalSum;

    gaspi_datatype_t commType = getGPIType<T>();

    SCAI_GASPI_CALL( gaspi_allreduce( &localSum, &globalSum, 1, GASPI_OP_SUM, commType, group, timeout ) )

    SCAI_LOG_DEBUG( logger, "sum: my value = " << value << ", sum = " << globalSum )

    return globalSum;
}

/* ---------------------------------------------------------------------------------- */
/*                    min reduction                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename T>
T GPICommunicator::minImpl( const T value ) const
{
    SCAI_REGION( "Communicator.GPI.min" )

    SCAI_LOG_DEBUG( logger, "minImpl: local value = " << value )

    T localMin = value;   // might be used also as buffer

    T globalMin = 0;

    gaspi_datatype_t commType = getGPIType<T>();

    SCAI_GASPI_CALL( gaspi_allreduce( &localMin, &globalMin, 1, GASPI_OP_MIN, commType, group, timeout ) )

    SCAI_LOG_DEBUG( logger, "minImpl: global value = " << globalMin )

    return globalMin;
}

/* ---------------------------------------------------------------------------------- */
/*                    max reduction                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename T>
T GPICommunicator::maxImpl( const T value ) const
{
    SCAI_REGION( "Communicator.GPI.max" )

    T localMax = value;   // might be used also as buffer

    T globalMax = 0;

    SCAI_LOG_DEBUG( logger, "maxImpl: local value = " << value )

    gaspi_datatype_t commType = getGPIType<T>();

    SCAI_GASPI_CALL( gaspi_allreduce( &localMax, &globalMax, 1, GASPI_OP_MAX, commType, group, timeout ) );

    SCAI_LOG_DEBUG( logger, "maxImpl: global value = " << globalMax )

    return globalMax;
}

/* ---------------------------------------------------------------------------------- */
/*      synchronize                                                                   */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::synchronize() const
{
    SCAI_REGION( "Communicator.GPI.synchronize" )

    SCAI_LOG_DEBUG( logger, *this << ": before barrier" )

    SCAI_GASPI_CALL( gaspi_barrier( GASPI_GROUP_ALL, timeout ) )

    SCAI_LOG_DEBUG( logger, *this << ": after barrier" )
}

/* ---------------------------------------------------------------------------------- */
/*      bcast                                                                         */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void GPICommunicator::bcastImpl( T val[], const IndexType n, const PartitionId root ) const
{
    SCAI_REGION( "Communicator.GPI.bcast" )

    if ( n <= 0 )
    {
        synchronize();
        return;
    }

    SegmentData<T> segment( this, n );

    SCAI_LOG_DEBUG( logger, *this << ": bcast<" << TypeTraits<T>::id() << ">, root = " << root << ", n = " << n )

    int distance = 1;

    if ( mRank == root )
    {
        segment.assign( val, n );
    }

    // SCAI_ASSERT_EQUAL_ERROR( 0, root )   // not yet supported for other root

    // Note: If root processor is not 0 an 'implicit' shift of '-root' is done

    while ( distance < mSize )    /* log NP (base 2) loop */
    {
        distance = 2 * distance;
    }

    const PartitionId gRank = ( mRank + mSize - root ) % mSize;

    while ( distance > 1 )
    { 
        int steph = distance;
        distance = distance / 2;

        if ( ( gRank % steph ) == 0 )
        { 
            PartitionId partner = gRank + distance;

            if ( partner < mSize ) 
            {
                // send partner data

                PartitionId target = ( partner + root ) % mSize;

                SCAI_LOG_DEBUG( logger, *this << ": bcast send to " << target )
                
                const IndexType remoteOffset = notifyWait( segment.getID(), target );
                const IndexType localOffset  = 0;
                remoteWrite( segment, localOffset, target, segment, remoteOffset, n );
                notify( segment.getID(), target, mRank + mSize, n );
            }
        }

        if ( ( gRank % steph ) == distance )
        {
            PartitionId partner = gRank - distance;

            // receive partner data

            PartitionId source = ( partner + root ) % mSize;

            SCAI_LOG_DEBUG( logger, *this << ": bcast recv from " << source )

            const IndexType offset = segment.getOffset();
            notify( segment.getID(), source, mRank, offset );
            const IndexType val = notifyWait( segment.getID(), source + mSize );
            SCAI_ASSERT_EQUAL_ERROR( n, val );
        }
    }

    SCAI_LOG_DEBUG( logger, *this << ": done bcast with n = " << n << ", root = " << root );

    if ( root != mRank )
    {
        segment.copyTo( val, n ); // copy from segment to ptr
    }

    wait();  // make sure that all data has been transferred
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myvals, n, root, allvals )                                           */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void GPICommunicator::scatterImpl( T myvals[], const IndexType n, const PartitionId root, const T allvals[] ) const
{
    SCAI_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
    SCAI_LOG_DEBUG( logger, *this << ": scatter of " << n << " elements, root = " << root )

    if ( n < 1 )
    {
        return;
    }

    SegmentData<T> srcSegment( this, n * mSize );  // for allVals
    SegmentData<T> dstSegment( this, n );          // for myVals

    if ( mRank == root )
    {
        srcSegment.assign( allvals, n * mSize); // copy from ptr to segment

        IndexType offset = srcSegment.getOffset();

        for ( PartitionId pid = 0; pid < mSize; ++pid )
        {
            if ( pid != mRank )
            {
                notify( srcSegment.getID(), pid, mRank, offset );
                SCAI_LOG_DEBUG( logger, *this << ": notified " << pid << " with offset " << offset )
            }
            offset += n;
        }
    }

    // each processor reads its contribution from root processor and notifies

    if ( root != mRank )
    {
        IndexType remoteOffset = notifyWait( srcSegment.getID(), root );
        remoteRead( dstSegment, 0, root, srcSegment, remoteOffset, n );
        notify( srcSegment.getID(), root, mRank + mSize, n );
    }
    else
    {
        remoteRead( dstSegment, 0, root, srcSegment, mRank * n + srcSegment.getOffset(), n );
    }

    wait();  // make sure that remoteRead has finished

    dstSegment.copyTo( myvals, n );  // each proc copies its values from dst

    for ( int i = 0; i < n; ++i )
    {
        SCAI_LOG_DEBUG( logger, "dstSegment[" << i << "] = " << dstSegment[i] )
        SCAI_LOG_DEBUG( logger, "myvals[" << i << "] = " << myvals[i] )
    }

    // root waits for all processor until data has been read

    if ( root == mRank )
    {
        for ( PartitionId pid = 0; pid < mSize; ++pid )
        {
            if ( pid == mRank )
            {
                continue;
            }

            IndexType writtenSize = notifyWait( srcSegment.getID(), pid + mSize );
            SCAI_ASSERT_EQUAL_ERROR( n, writtenSize )
        }
    }

    SCAI_LOG_DEBUG( logger, *this << ": done scatter of " << n << " elements, root = " << root )
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myvals, n, root, allvals, sizes )                                    */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void GPICommunicator::scatterVImpl(
    T myvals[],
    const IndexType n,
    const PartitionId root,
    const T allvals[],
    const IndexType sizes[] ) const
{
    SCAI_REGION( "Communicator.GPI.scatterV" )

    SCAI_ASSERT_ERROR( root < getSize(), "illegal root, root = " << root )

    int totalSize = 1;

    // Array sizes might be available only on root processor

    if ( mRank == root )
    {
        totalSize = 0;

        for ( int i = 0; i < mSize; ++i )
        {
            totalSize += sizes[i];
        }
    }

    SCAI_LOG_DEBUG( logger, *this << ": scatterV<" << TypeTraits<T>::id() << ">, n = "
                            << n << ", total = " << totalSize )

    SegmentData<T> srcSegment( this, totalSize );   // for allVals
    SegmentData<T> dstSegment( this, n );           // for myVals

    IndexType rootOffset = 0;  // only used if mRank == root and n > 0

    if ( mRank == root )
    {
        srcSegment.assign( allvals, totalSize); // copy from ptr to segment

        IndexType offset = srcSegment.getOffset();

        for ( PartitionId pid = 0; pid < mSize; ++pid )
        {
            int size = sizes[pid];

            if ( size == 0 )
            {
                continue;
            }

            if ( pid != mRank )
            {
                notify( srcSegment.getID(), pid, mRank, offset );
                SCAI_LOG_DEBUG( logger, *this << ": notified " << pid << " with offset " << offset )
            }
            else
            {
                rootOffset = offset;
            }

            offset += size;
        }
    }

    // each processor reads its contribution from root processor and notifies

    if ( n > 0 )
    {
        if ( root != mRank )
        {
            IndexType remoteOffset = notifyWait( srcSegment.getID(), root );
            remoteRead( dstSegment, 0, root, srcSegment, remoteOffset, n );
            notify( srcSegment.getID(), root, mRank + mSize, n );
        }
        else
        {
            remoteRead( dstSegment, 0, root, srcSegment, rootOffset, n );
        }

        wait();  // make sure that remoteRead has finished

        dstSegment.copyTo( myvals, n );  // each proc copies its values from dst
    }

    // root waits for all processor until data has been read

    if ( root == mRank )
    {
        for ( PartitionId pid = 0; pid < mSize; ++pid )
        {
            if ( pid == mRank )
            {
                continue;
            }

            IndexType size = sizes[pid];

            if ( size > 0 )
            {
                IndexType writtenSize = notifyWait( srcSegment.getID(), pid + mSize );
                SCAI_ASSERT_EQUAL_ERROR( size, writtenSize )
            }
        }
    }

    SCAI_LOG_DEBUG( logger, *this << ": scatterV<" << TypeTraits<T>::id() << ">, ready" )
}

template<typename T>
void GPICommunicator::remoteWrite( const SegmentData<T>& localSegment, const IndexType localOffset, const PartitionId remP,  
                                   SegmentData<T>& remSegment, const IndexType remoteOffset, const IndexType size ) const
{
    SCAI_LOG_DEBUG( logger, *this << " remoteWrite " << remSegment  << ":" << remoteOffset
                             << " @ p = " << remP << ", size = " << size
                             << " = local " << localSegment << ":" << localOffset 
                             << ", size = " << size )

    // SegmentData might have an additonal offset to be considered

    const gaspi_offset_t localSegOffset = ( localSegment.getOffset() + localOffset ) * sizeof( T );
    const gaspi_offset_t remSegOffset   = remoteOffset * sizeof( T );

    // gaspi_write needs segment ids

    const gaspi_segment_id_t localID = localSegment.getID();
    const gaspi_segment_id_t remID   = remSegment.getID();

    gaspi_pointer_t ptr;

    SCAI_GASPI_CALL ( gaspi_segment_ptr ( localID, &ptr ) )

    gaspi_pointer_t s = static_cast<char*>( ptr ) + localSegOffset;

    SCAI_GASPI_CALL ( gaspi_segment_ptr ( remID, &ptr ) )

    gaspi_pointer_t t = static_cast<char*>( ptr ) + remSegOffset;

    SCAI_LOG_DEBUG( logger, *this << ": localAddr = " << s << ", remAddr = " << t )

    SCAI_GASPI_CALL( gaspi_write( localID, localSegOffset, remP, remID, remSegOffset, size * sizeof( T ), mQueueID, timeout ) )
}

template<typename T>
void GPICommunicator::localWrite( const SegmentData<T>& srcSegment, const IndexType srcOffset,
                                   SegmentData<T>& dstSegment, const IndexType dstOffset, const IndexType size ) const
{
    SCAI_LOG_DEBUG( logger, *this << " dest " << dstSegment  << ":" << dstOffset
                             << ", size = " << size
                             << " = local " << srcSegment << ":" << srcOffset
                             << ", size = " << size )

    const T* srcPtr = srcSegment.get() + srcOffset;
    T* dstPtr = dstSegment.get() + dstOffset;
    memcpy( dstPtr, srcPtr, size * sizeof( T ) );
}

template<typename T>
void GPICommunicator::remoteRead( SegmentData<T>& localSegment, const IndexType localOffset, const PartitionId remP,
                                  const SegmentData<T>& remSegment, const IndexType remoteOffset, const IndexType size ) const
{
    SCAI_LOG_DEBUG( logger, *this << ": remoteRead<" << TypeTraits<T>::id() 
                              << ">, local " << localSegment << ":" << localOffset 
                              << " = remote " << remSegment  << ":" << remoteOffset
                              << " @ p = " << remP << ", size = " << size)

    const gaspi_offset_t localSegOffset = ( localSegment.getOffset() + localOffset ) * sizeof( T );
    const gaspi_offset_t remSegOffset = remoteOffset * sizeof( T );  

    // Important: remoteOffset must have added remoteSegment.getOffset() of the remote processor
    
    const gaspi_segment_id_t localID = localSegment.getID();
    const gaspi_segment_id_t remID = remSegment.getID();

    SCAI_GASPI_CALL( gaspi_read( localID, localSegOffset, remP, remID, remSegOffset, size * sizeof( T ), mQueueID, timeout ) )
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allvals, n, root, myvals )                                            */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void GPICommunicator::gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const
{
    SCAI_REGION( "Communicator.GPI.gather" )

    SCAI_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )

    SCAI_LOG_DEBUG( logger, *this << ": gather<" << TypeTraits<T>::id() << 
                              " of " << n << " elements, root = " << root )

    if ( n < 1 )
    {
        return;
    }

    SegmentData<T> srcSegment( this, n );          // keeps myvals to be accessed by root
    SegmentData<T> dstSegment( this, n * mSize );  // keeps allvals on root

    srcSegment.assign( myvals, n ); // copy from ptr to segment

    if ( root == mRank )
    {
        IndexType offset = dstSegment.getOffset();  

        for ( PartitionId pid = 0; pid < mSize; ++pid )
        {
            if ( pid != mRank )
            {
                notify( dstSegment.getID(), pid, mRank, offset );
                SCAI_LOG_DEBUG( logger, *this << ": notified " << pid << " with offset " << offset )
            }
            offset += n;
        }
    }

    // each processor write its contribution to root processor and notifies

    if ( mRank != root )
    {
        SCAI_LOG_DEBUG( logger, *this << ": wait for ready send data to root" )
        IndexType remoteOffset = notifyWait( dstSegment.getID(), root );
        remoteWrite( srcSegment, 0, root, dstSegment, remoteOffset, n );
        notify( dstSegment.getID(), root, mRank + mSize, n );
        SCAI_LOG_DEBUG( logger, *this << ": sent to root" )
    }
    else
    {
        SCAI_LOG_DEBUG( logger, *this << " as root just copy my data" )
        localWrite( srcSegment, 0, dstSegment, mRank * n, n );
    }

    // root waits for all processor until data is available

    if ( root == mRank )
    {
        for ( PartitionId pid = 0; pid < mSize; ++pid )
        {
            if ( pid == mRank )
            {
                continue;
            }

            IndexType writtenSize = notifyWait( dstSegment.getID(), pid + mSize );
            SCAI_ASSERT_EQUAL_ERROR( n, writtenSize )
        }

        dstSegment.copyTo( allvals, n * mSize );
    }
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allvals, n, root, myvals, sizes )                                     */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void GPICommunicator::gatherVImpl(
    T allvals[],
    const IndexType n,
    const PartitionId root,
    const T myvals[],
    const IndexType sizes[] ) const
{
    SCAI_REGION( "Communicator.GPI.gatherV" )

    SCAI_ASSERT_ERROR( root < getSize(), "illegal root, root = " << root )
    SCAI_LOG_DEBUG( logger, *this << ": gather, root = " << root )

    int totalSize = 1;   // default value for non-root processors 

    if ( mRank == root )
    {
        totalSize = 0;

        for ( int i = 0; i < mSize; ++i )
        {
            totalSize += sizes[i];
        }
    }

    SCAI_LOG_DEBUG( logger, *this << ": gatherV<" << TypeTraits<T>::id() << ">, n = "
                               << n << ", total = " << totalSize )

    // Note: segment data is now safe against different values of n, totalSize

    SegmentData<T> srcSegment( this, n );
    SegmentData<T> dstSegment( this, totalSize );   // allValues on root

    srcSegment.assign( myvals, n ); // copy from ptr to segment

    IndexType rootOffset = 0;

    if ( root == mRank )
    {
        IndexType offset = dstSegment.getOffset();

        for ( PartitionId pid = 0; pid < mSize; ++pid )
        {
            int size = sizes[pid];

            if ( size > 0 )
            {
                if ( pid != mRank )
                {
                    notify( dstSegment.getID(), pid, mRank, offset );
                }
                else
                {
                    rootOffset = offset;
                }
            }
            offset += size;
        }
    }

    // each processor write its contribution to root processor and notifies

    if ( n > 0 )
    {
        if ( mRank != root )
        {
            IndexType remoteOffset = notifyWait( dstSegment.getID(), root );
            remoteWrite( srcSegment, 0, root, dstSegment, remoteOffset, n );
            notify( dstSegment.getID(), root, mRank + mSize, n );
        }
        else
        {
            // Note: rootOffset is physical offset, so do not use localWrite
            remoteWrite( srcSegment, 0, root, dstSegment, rootOffset, n );
        }
    }

    // root waits for all processor until data is available

    if ( root == mRank )
    {
        for ( PartitionId pid = 0; pid < mSize; ++pid )
        {
            if ( pid == mRank ) 
            {
                continue;
            }

            int size = sizes[pid];

            if ( size > 0 )
            {
                IndexType writtenSize = notifyWait( dstSegment.getID(), pid + mSize );
                SCAI_ASSERT_EQUAL_ERROR( size, writtenSize )
            }
        }

        dstSegment.copyTo( allvals, totalSize );
    }

    wait();
}

/* ---------------------------------------------------------------------------------- */
/*           maxloc                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void GPICommunicator::maxlocImpl( T& val, IndexType& location, PartitionId root ) const
{
    // get max values
    T myVal = val;
    val = max( val );

    IndexType itsMyValue = -1;
    if( myVal == val )
    {
        itsMyValue = location;
    }

    vector<IndexType> locations( mSize );

    gatherImpl( locations.data(), 1, root, &itsMyValue );

    if( mRank == root )
    {
        // find first property not -1
        for( int i = 0; i < mSize; ++i )
        {
            if( locations[i] != -1 )
            {
                location = i;
                break;
            }
        }
    }

    bcast( &location, 1, root );
}

/* ---------------------------------------------------------------------------------- */
/*           swap                                                                     */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void GPICommunicator::swapImpl( T val[], const IndexType n, PartitionId partner ) const
{
    // this is now safe as there is no global synchronization 

    if ( partner == mRank )
    {
        return;
    }

    SegmentData<T> srcSegment( this, n );   // used remote, will be read
    SegmentData<T> dstSegment( this, n );   // used locally
  
    SCAI_LOG_DEBUG( logger, *this << ": swap with processor " << partner << ", n = " << n )

    // Tell my partner where to write

    const IndexType offset = dstSegment.getOffset();
    notify( dstSegment.getID(), partner, mRank, offset );

    srcSegment.assign( val, n );

    // wait for ready to write of my partner
    const IndexType remoteOffset = notifyWait( dstSegment.getID(), partner );
    const IndexType localOffset  = 0;
    remoteWrite( srcSegment, localOffset, partner, dstSegment, remoteOffset, n );
    notify( dstSegment.getID(), partner, mRank + mSize, n );

    // wait for acknowledge that my data has been written
    const IndexType nn = notifyWait( dstSegment.getID(), partner + mSize );
    SCAI_ASSERT_EQUAL_ERROR( n, nn );

    dstSegment.copyTo( val, n );
}

/* ---------------------------------------------------------------------------------- */
/*          getCommunicationContext                                                   */
/* ---------------------------------------------------------------------------------- */

hmemo::ContextPtr GPICommunicator::getCommunicationContext( const hmemo::_HArray& ) const
{
    // for comparison: 
    // return ContextFactory::getContext( Context::Host );
    
    return hmemo::Context::getHostPtr();
}

/* ---------------------------------------------------------------------------------- */
/*           writeAt                                                                  */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::writeAt( std::ostream& stream ) const
{
    // Info about rank and size of the communicator is very useful
    stream << "GPI(" << mRank << ":" << mSize << ")";
}

} // namespace dmemo

} // namespace scai
