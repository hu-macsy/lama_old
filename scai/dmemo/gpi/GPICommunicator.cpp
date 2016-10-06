/**
 * @file GPICommunicator.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief GPICommunicator.cpp
 * @author Thomas Brandes, Lauretta Schubert
 * @date 06.05.2014
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
#include <scai/common/Settings.hpp>

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

Communicator::CommunicatorKind GPICommunicator::createValue()
{
    return GPI;
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
    : Communicator( GPI ),
      mQueueID( 0 ),
      mThreadSafetyLevel( Communicator::Funneled )
{
    SCAI_TRACE_SCOPE( false )   // switch off tracing in this scope as it might call this constructor again

    // set gaspi_printf as print routine for logging, no more needed for GPI 2.1.3.0

    // logging::GenLogger::myPrintf = &gaspi_printf;

    SCAI_LOG_DEBUG( logger, "GPICommunicator(): call init" )
    gaspi_timeout_t init_timeout = 10000;   // in milli seconds
    SCAI_GASPI_CALL( gaspi_proc_init( init_timeout ) )
    SCAI_GASPI_CALL( gaspi_proc_num( &gSize ) )
    SCAI_GASPI_CALL( gaspi_proc_rank( &gRank ) )
    mRank = gRank;
    mSize = gSize;
    SCAI_LOG_INFO( logger, "GASPI proc " << mRank << " of " << mSize << " started" )
    gaspi_number_t num_notifications;
    SCAI_GASPI_CALL( gaspi_notification_num ( &num_notifications ) )
    mMaxNotifications = num_notifications;
    SCAI_ASSERT_ERROR( mMaxNotifications >= 2 * mSize,
                       "# notifications = " << mMaxNotifications << " not sufficient"
                       << ", need at least 2 * " << mSize )

    // the following call is for convenience to get correct node rank by hostnames
    // Note: mRank, mSize must be set correctly before

    setNodeData();   // determine mNodeRank, mNodeSize, also uses some communication pattern

    std::ostringstream commVal;
    commVal << *this;
    common::Settings::putEnvironment( "SCAI_COMM", commVal.str().c_str() );
    common::Settings::putEnvironment( "SCAI_RANK", mRank );
}

/* ---------------------------------------------------------------------------------- */

#define GPI_MAX_PROCESSOR_NAME 512

/* --------------------------------------------------------------- */

void GPICommunicator::getProcessorName( char* name ) const
{
    size_t len = maxProcessorName();

    memset( name, 0, len * sizeof( char ) );

    gethostname( name,len );
}

/* --------------------------------------------------------------- */

size_t GPICommunicator::maxProcessorName() const
{
    return GPI_MAX_PROCESSOR_NAME;
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

/* ---------------------------------------------------------------------------------- */
/*              getGPIType                                                            */
/* ---------------------------------------------------------------------------------- */

gaspi_datatype_t GPICommunicator::getGPIType( common::scalar::ScalarType stype )
{
    switch ( stype )
    {
        case common::scalar::INT                 : return GASPI_TYPE_INT;
        case common::scalar::FLOAT               : return GASPI_TYPE_FLOAT;
        case common::scalar::DOUBLE              : return GASPI_TYPE_DOUBLE;
        case common::scalar::UNSIGNED_LONG       : return GASPI_TYPE_ULONG;

        default:
             COMMON_THROWEXCEPTION( "No GPI Type specified for " << stype )
             return GASPI_TYPE_DOUBLE;
    }
}

/* ---------------------------------------------------------------------------------- */
/*                                     wait                                           */
/* ---------------------------------------------------------------------------------- */

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

    SegmentData sendSegment( common::scalar::INDEX_TYPE, this, mSize, const_cast<IndexType*>( sendSizes ) );  // locally used for sending
    SegmentData recvSegment( common::scalar::INDEX_TYPE, this, mSize, recvSizes );  // other procs do remote write here

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
            // no more possible: SCAI_LOG_DEBUG( logger, *this << ": write " << toP << " " << sendSegment[ localOffset ] )
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
        SCAI_ASSERT_EQ_ERROR( 1, quantity, "serious mismatch" );
    }

    recvSegment.copyTo( recvSizes, mSize );  // copy from recvSegment back to recvSizes
    SCAI_LOG_DEBUG( logger, *this << ": all2all done" )
}

/* ---------------------------------------------------------------------------------- */

void GPICommunicator::all2allvImpl(
    void* recvBuffer[],
    IndexType recvCount[],
    void* sendBuffer[],
    IndexType sendCount[],
    common::scalar::ScalarType ) const
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
    SCAI_GASPI_CALL ( gaspi_notify_waitsome( segID, begin_id, 1 , &notID, GASPI_BLOCK ) );
    gaspi_notification_t notifyVal;
    SCAI_GASPI_CALL( gaspi_notify_reset( segID, notID, &notifyVal ) );
    SCAI_LOG_DEBUG( logger, *this << ": got notification at " << pos << ", val = " << ( notifyVal - 1 ) );
    return static_cast<IndexType>( notifyVal - 1 );
}

/* ---------------------------------------------------------------------------------- */
/*                                exchange by plan                                    */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::exchangeByPlanImpl(
    void* recvData,
    const CommunicationPlan& recvPlan,
    const void* sendData,
    const CommunicationPlan& sendPlan,
    common::scalar::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.GPI.exchangeByPlan" )
    SCAI_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    SCAI_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )
    SCAI_LOG_DEBUG( logger, *this << ": exchange by plan for values of type " << stype
                    << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )

    int noReceives = recvPlan.size();
    int noSends = sendPlan.size();
    int totalSendSize = sendPlan.totalQuantity();
    int totalRecvSize = recvPlan.totalQuantity();

    // setup segment data
    SegmentData srcDataSegment( stype, this, totalSendSize, const_cast<void*>( sendData ) );
    SegmentData dstDataSegment( stype, this, totalRecvSize, recvData );

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
        SCAI_ASSERT_EQ_ERROR( quantity, val, "notification mismatch" );
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

tasking::SyncToken* GPICommunicator::exchangeByPlanAsyncImpl(
    void* const recvData,
    const CommunicationPlan& recvPlan,
    const void* const sendData,
    const CommunicationPlan& sendPlan,
    const common::scalar::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.GPI.exchangeByPlanAsync" )
    SCAI_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    SCAI_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )
    // Some global information from the communication plans
    int noReceives = recvPlan.size();
    int noSends = sendPlan.size();
    int totalSendSize = sendPlan.totalQuantity();
    int totalRecvSize = recvPlan.totalQuantity();

    SCAI_LOG_INFO( logger, *this << ": async exchange for values of type " << stype
                   << ", send to " << sendPlan.size() << " processors " << totalSendSize << " values"
                   << ", recv from " << recvPlan.size() << " processors " << totalRecvSize << " values" )

    // setup segment data

    shared_ptr<SegmentData> srcDataSegment( new SegmentData( stype, this, totalSendSize ) );
    shared_ptr<SegmentData> dstDataSegment( new SegmentData( stype, this, totalRecvSize ) );

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
    pSyncToken->pushRoutine( common::bind( &SegmentData::copyTo, dstDataSegment.get(), recvData, totalRecvSize ) );
    // Also wait that all remote writes have been finished
    pSyncToken->pushRoutine( common::bind( &GPICommunicator::wait, this ) );
    return pSyncToken.release();
}

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

IndexType GPICommunicator::shiftImpl(
    void* recvVals,
    const IndexType recvSize,
    const PartitionId source,
    const void* sendVals,
    const IndexType sendSize,
    const PartitionId dest,
    common::scalar::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.GPI.shift" )

    SCAI_ASSERT_ERROR( source != mRank, "source must not be this partition" )
    SCAI_ASSERT_ERROR( dest != mRank, "dest must not be this partition" )
    SCAI_ASSERT_ERROR( sendSize <= recvSize, "too large send" )

    SCAI_LOG_INFO( logger,
                   *this << ": shift<" << stype << ">" << 
                   ", recv from " << source << " max " << recvSize << " values " << 
                   ", send to " << dest << " " << sendSize << " values." )

    SegmentData srcDataSegment( stype, this, recvSize, const_cast<void*>( sendVals ) );
    SegmentData dstDataSegment( stype, this, recvSize, recvVals );

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

tasking::SyncToken* GPICommunicator::shiftAsyncImpl(
    void* recvVals,
    const PartitionId source,
    const void* sendVals,
    const PartitionId dest,
    const IndexType size,
    common::scalar::ScalarType stype ) const
{
    // SCAI_REGION( "Communicator.GPI.shiftAsyncImpl" )

    SCAI_LOG_DEBUG( logger,
                    *this << ": recv from " << source << ", send to " << dest << ", both " << size << " values." )
    SCAI_ASSERT_ERROR( source != mRank, "source must not be this partition" )
    SCAI_ASSERT_ERROR( dest != mRank, "dest must not be this partition" )

    shared_ptr<SegmentData> srcSegment( new SegmentData( stype, this, size ) );
    shared_ptr<SegmentData> dstSegment( new SegmentData( stype, this, size ) );

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
    pSyncToken->pushRoutine( common::bind( &SegmentData::copyTo, dstSegment.get(), recvVals, size ) );
    // Keep SegmentData until synchronization of token
    pSyncToken->pushToken( srcSegment );
    pSyncToken->pushToken( dstSegment );
    return pSyncToken.release();
}

/* ---------------------------------------------------------------------------------- */
/*              reduce                                                                */
/* ---------------------------------------------------------------------------------- */

/** Common method for a call of gaspi_allreduce_user */

void GPICommunicator::reduce( 
    void* outValues, 
    const void* inValues, 
    const IndexType n, 
    gaspi_reduce_operation_t op, 
    gaspi_size_t elem_size ) const
{
    gaspi_state_t state = NULL;

    SCAI_GASPI_CALL( gaspi_barrier( group, timeout ) )   // really needed here, otherwise wrong results

    if ( inValues == outValues )
    {
        // unclear whether IN_PLACE is supported, for safety we allocate a copy of inValues 

        size_t tmpSize = n * elem_size;

        common::scoped_array<char> tmp ( new char[ tmpSize ] );

        memcpy( tmp.get(), inValues,  tmpSize );

        SCAI_GASPI_CALL( gaspi_allreduce_user( tmp.get(), outValues, n, elem_size, op, state, group, timeout ) )
    }
    else
    {
        SCAI_GASPI_CALL( gaspi_allreduce_user( const_cast<void*>( inValues ), outValues, n, elem_size, op, state, group, timeout ) )
    }

    SCAI_GASPI_CALL( gaspi_barrier( group, timeout ) )
}

/* ---------------------------------------------------------------------------------- */
/*              sum                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
static gaspi_return_t
reduce_sum( gaspi_pointer_t op1,
           gaspi_pointer_t op2,
           gaspi_pointer_t result,
           gaspi_state_t,
           const gaspi_number_t numValues,
           const gaspi_size_t elemSize,
           gaspi_timeout_t )
{
   SCAI_ASSERT_EQ_ERROR( elemSize, sizeof( ValueType ), "serious type mismatch" )

   ValueType* t_op1 = reinterpret_cast<ValueType*>( op1 );
   ValueType* t_op2 = reinterpret_cast<ValueType*>( op2 );
   ValueType* t_res = reinterpret_cast<ValueType*>( result );

   for ( gaspi_number_t i = 0; i< numValues; ++i )
   {
       t_res[i] = t_op1[i] + t_op2[i];
   }

   return GASPI_SUCCESS;
}

void GPICommunicator::sumImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.GPI.sum" )

    gaspi_size_t elem_size = common::typeSize( stype );

    SCAI_LOG_INFO( logger, *this << ": sumImpl: #values = " << n << ", elemSize = " << elem_size )

    gaspi_reduce_operation_t op;

    switch( stype )
    {
        case common::scalar::INT           : op = reduce_sum<int>; break;
        case common::scalar::UNSIGNED_INT  : op = reduce_sum<unsigned int>; break;
        case common::scalar::LONG          : op = reduce_sum<long>; break;
        case common::scalar::UNSIGNED_LONG : op = reduce_sum<unsigned long>; break;
        case common::scalar::FLOAT         : op = reduce_sum<float>; break;
        case common::scalar::DOUBLE        : op = reduce_sum<double>; break;
        case common::scalar::LONG_DOUBLE   : op = reduce_sum<long double>; break;
#ifdef SCAI_COMPLEX_SUPPORTED
        case common::scalar::COMPLEX             : op = reduce_sum<ComplexFloat>; break;
        case common::scalar::DOUBLE_COMPLEX      : op = reduce_sum<ComplexDouble>; break;
        case common::scalar::LONG_DOUBLE_COMPLEX : op = reduce_sum<ComplexLongDouble>; break;
#endif
        default:
           COMMON_THROWEXCEPTION( "Unsupported reduction type for sum: " << stype );
    }

    reduce( outValues, inValues, n, op, elem_size );
}

/* ---------------------------------------------------------------------------------- */
/*              min                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
static gaspi_return_t
reduce_min( gaspi_pointer_t op1,
            gaspi_pointer_t op2,
            gaspi_pointer_t result,
            gaspi_state_t,
            const gaspi_number_t numValues,
            const gaspi_size_t elemSize,
            gaspi_timeout_t )
{
   SCAI_ASSERT_EQ_ERROR( elemSize, sizeof( ValueType ), "serious type mismatch" )

   ValueType* t_op1 = reinterpret_cast<ValueType*>( op1 );
   ValueType* t_op2 = reinterpret_cast<ValueType*>( op2 );
   ValueType* t_res = reinterpret_cast<ValueType*>( result );

   for ( gaspi_number_t i = 0; i< numValues; ++i )
   {
       t_res[i] = t_op2[i] < t_op1[i] ? t_op2[i] : t_op1[i];
   }

   return GASPI_SUCCESS;
}

void GPICommunicator::minImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.GPI.sum" )

    gaspi_size_t elem_size = common::typeSize( stype );

    SCAI_LOG_INFO( logger, *this << ": minImpl: #values = " << n << ", elemSize = " << elem_size )

    gaspi_reduce_operation_t op;

    switch( stype )
    {
        case common::scalar::INT           : op = reduce_min<int>; break;
        case common::scalar::UNSIGNED_INT  : op = reduce_min<unsigned int>; break;
        case common::scalar::LONG          : op = reduce_min<long>; break;
        case common::scalar::UNSIGNED_LONG : op = reduce_min<unsigned long>; break;
        case common::scalar::FLOAT         : op = reduce_min<float>; break;
        case common::scalar::DOUBLE        : op = reduce_min<double>; break;
        case common::scalar::LONG_DOUBLE   : op = reduce_min<long double>; break;
#ifdef SCAI_COMPLEX_SUPPORTED
        case common::scalar::COMPLEX             : op = reduce_min<ComplexFloat>; break;
        case common::scalar::DOUBLE_COMPLEX      : op = reduce_min<ComplexDouble>; break;
        case common::scalar::LONG_DOUBLE_COMPLEX : op = reduce_min<ComplexLongDouble>; break;
#endif
        default:
           COMMON_THROWEXCEPTION( "Unsupported reduction type for min: " << stype );
    }

    reduce( outValues, inValues, n, op, elem_size );
}

/* ---------------------------------------------------------------------------------- */
/*              max                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
static gaspi_return_t
reduce_max( gaspi_pointer_t op1,
            gaspi_pointer_t op2,
            gaspi_pointer_t result,
            gaspi_state_t,
            const gaspi_number_t numValues,
            const gaspi_size_t elemSize,
            gaspi_timeout_t )
{
   SCAI_ASSERT_EQ_ERROR( elemSize, sizeof( ValueType ), "serious type mismatch" )

   ValueType* t_op1 = reinterpret_cast<ValueType*>( op1 );
   ValueType* t_op2 = reinterpret_cast<ValueType*>( op2 );
   ValueType* t_res = reinterpret_cast<ValueType*>( result );

   for ( gaspi_number_t i = 0; i< numValues; ++i )
   {
       t_res[i] = t_op2[i] > t_op1[i] ? t_op2[i] : t_op1[i];
   }

   return GASPI_SUCCESS;
}

void GPICommunicator::maxImpl( void* outValues, const void* inValues, const IndexType n, common::scalar::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.GPI.sum" )

    gaspi_size_t elem_size = common::typeSize( stype );

    SCAI_LOG_INFO( logger, *this << ": maxImpl: #values = " << n << ", elemSize = " << elem_size )

    gaspi_reduce_operation_t op;

    switch( stype )
    {
        case common::scalar::INT           : op = reduce_max<int>; break;
        case common::scalar::UNSIGNED_INT  : op = reduce_max<unsigned int>; break;
        case common::scalar::LONG          : op = reduce_max<long>; break;
        case common::scalar::UNSIGNED_LONG : op = reduce_max<unsigned long>; break;
        case common::scalar::FLOAT         : op = reduce_max<float>; break;
        case common::scalar::DOUBLE        : op = reduce_max<double>; break;
        case common::scalar::LONG_DOUBLE   : op = reduce_max<long double>; break;
#ifdef SCAI_COMPLEX_SUPPORTED
        case common::scalar::COMPLEX             : op = reduce_max<ComplexFloat>; break;
        case common::scalar::DOUBLE_COMPLEX      : op = reduce_max<ComplexDouble>; break;
        case common::scalar::LONG_DOUBLE_COMPLEX : op = reduce_max<ComplexLongDouble>; break;
#endif
        default:
           COMMON_THROWEXCEPTION( "Unsupported reduction type for max: " << stype );
    }

    reduce( outValues, inValues, n, op, elem_size );
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

void GPICommunicator::bcastImpl( void* val, const IndexType n, const PartitionId root, common::scalar::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.GPI.bcast" )

    if ( n <= 0 )
    {
        synchronize();
        return;
    }

    SegmentData segment( stype, this, n );

    SCAI_LOG_DEBUG( logger, *this << ": bcast<" << stype << ">, root = " << root << ", n = " << n )

    int distance = 1;

    if ( mRank == root )
    {
        segment.assign( val, n );
    }

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
            SCAI_ASSERT_EQ_ERROR( n, val, "notify mismatch" );
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
/*      scatter( myVals, n, root, allVals )                                           */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::scatterImpl(
    void* myVals,
    const IndexType n,
    const PartitionId root,
    const void* allVals,
    common::scalar::ScalarType stype ) const
{
    SCAI_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
    SCAI_LOG_DEBUG( logger, *this << ": scatter of " << n << " elements, root = " << root )

    if ( n < 1 )
    {
        return;
    }

    SegmentData srcSegment( stype, this, n * mSize );  // for allVals
    SegmentData dstSegment( stype, this, n );          // for myVals

    if ( mRank == root )
    {
        srcSegment.assign( reinterpret_cast<const double*>( allVals ), n * mSize ); // copy from ptr to segment
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
    dstSegment.copyTo( myVals, n );  // each proc copies its values from dst

    for ( int i = 0; i < n; ++i )
    {
        // SCAI_LOG_DEBUG( logger, "dstSegment[" << i << "] = " << dstSegment[i] )
        // SCAI_LOG_DEBUG( logger, "myVals[" << i << "] = " << myVals[i] )
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
            SCAI_ASSERT_EQ_ERROR( n, writtenSize, "notify mismatch" )
        }
    }

    SCAI_LOG_DEBUG( logger, *this << ": done scatter of " << n << " elements, root = " << root )
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myVals, n, root, allVals, sizes )                                    */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::scatterVImpl(
    void* myVals,
    const IndexType n,
    const PartitionId root,
    const void* allVals,
    const IndexType sizes[],
    common::scalar::ScalarType stype ) const
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

    SCAI_LOG_DEBUG( logger, *this << ": scatterV<" << stype << ">, n = " << n << ", total = " << totalSize )

    SegmentData srcSegment( stype, this, totalSize );   // for allVals
    SegmentData dstSegment( stype, this, n );           // for myVals

    IndexType rootOffset = 0;  // only used if mRank == root and n > 0

    if ( mRank == root )
    {
        srcSegment.assign( allVals, totalSize ); // copy from ptr to segment
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
        dstSegment.copyTo( myVals, n );  // each proc copies its values from dst
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
                SCAI_ASSERT_EQ_ERROR( size, writtenSize, "notify mismatch" )
            }
        }
    }

    SCAI_LOG_DEBUG( logger, *this << ": scatterV<" << stype << ">, ready" )
}

/* ---------------------------------------------------------------------------------- */
/*      remoteWrite                                                                   */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::remoteWrite( const SegmentData& localSegment, const IndexType localOffset, const PartitionId remP,
                                   SegmentData& remSegment, const IndexType remoteOffset, const IndexType size ) const
{
    SCAI_LOG_DEBUG( logger, *this << " remoteWrite " << remSegment  << ":" << remoteOffset
                    << " @ p = " << remP << ", size = " << size
                    << " = local " << localSegment << ":" << localOffset
                    << ", size = " << size )

    SCAI_ASSERT_EQ_ERROR( localSegment.scalarType(), remSegment.scalarType(), "type mismatch" );

    size_t typeSize = localSegment.typeSize();

    // SegmentData might have an additonal offset to be considered
    const gaspi_offset_t localSegOffset = ( localSegment.getOffset() + localOffset ) * typeSize;
    const gaspi_offset_t remSegOffset   = remoteOffset * typeSize;
    // gaspi_write needs segment ids
    const gaspi_segment_id_t localID = localSegment.getID();
    const gaspi_segment_id_t remID   = remSegment.getID();
    gaspi_pointer_t ptr;
    SCAI_GASPI_CALL ( gaspi_segment_ptr ( localID, &ptr ) )
    gaspi_pointer_t s = static_cast<char*>( ptr ) + localSegOffset;
    SCAI_GASPI_CALL ( gaspi_segment_ptr ( remID, &ptr ) )
    gaspi_pointer_t t = static_cast<char*>( ptr ) + remSegOffset;
    SCAI_LOG_DEBUG( logger, *this << ": localAddr = " << s << ", remAddr = " << t )
    SCAI_GASPI_CALL( gaspi_write( localID, localSegOffset, remP, remID, remSegOffset, size * typeSize, mQueueID, timeout ) )
}

/* ---------------------------------------------------------------------------------- */
/*      localWrite                                                                    */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::localWrite( const SegmentData& srcSegment, const IndexType srcOffset,
                                  SegmentData& dstSegment, const IndexType dstOffset, const IndexType size ) const
{
    SCAI_LOG_DEBUG( logger, *this << " dest " << dstSegment  << ":" << dstOffset
                    << ", size = " << size
                    << " = local " << srcSegment << ":" << srcOffset
                    << ", size = " << size )

    SCAI_ASSERT_EQ_ERROR( srcSegment.scalarType(), dstSegment.scalarType(), "type mismatch" );

    const void* srcPtr = srcSegment.get( srcOffset );
    void* dstPtr = dstSegment.get( dstOffset );

    size_t typeSize = srcSegment.typeSize();

    memcpy( dstPtr, srcPtr, size * typeSize );
}

/* ---------------------------------------------------------------------------------- */
/*      remoteRead                                                                    */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::remoteRead( SegmentData& localSegment, const IndexType localOffset, const PartitionId remP,
                                  const SegmentData& remSegment, const IndexType remoteOffset, const IndexType size ) const
{
    SCAI_LOG_DEBUG( logger, *this << ": remoteRead<" << remSegment.scalarType() 
                    << ">, local " << localSegment << ":" << localOffset
                    << " = remote " << remSegment  << ":" << remoteOffset
                    << " @ p = " << remP << ", size = " << size )

    SCAI_ASSERT_EQ_ERROR( localSegment.scalarType(), remSegment.scalarType(), "type mismatch" );

    size_t typeSize = localSegment.typeSize();

    const gaspi_offset_t localSegOffset = ( localSegment.getOffset() + localOffset ) * typeSize;
    const gaspi_offset_t remSegOffset = remoteOffset * typeSize;
    // Important: remoteOffset must have added remoteSegment.getOffset() of the remote processor
    const gaspi_segment_id_t localID = localSegment.getID();
    const gaspi_segment_id_t remID = remSegment.getID();
    SCAI_GASPI_CALL( gaspi_read( localID, localSegOffset, remP, remID, remSegOffset, size * typeSize, mQueueID, timeout ) )
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allVals, n, root, myVals )                                            */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::gatherImpl(
    void* allVals,
    const IndexType n,
    const PartitionId root,
    const void* myVals,
    const common::scalar::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.GPI.gather" )
    SCAI_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )

    SCAI_LOG_DEBUG( logger, *this << ": gather<" << stype << ">" <<
                    " of " << n << " elements, root = " << root )

    if ( n < 1 )
    {
        return;
    }

    SegmentData srcSegment( stype, this, n );          // keeps myVals to be accessed by root
    SegmentData dstSegment( stype, this, n * mSize );  // keeps allVals on root

    srcSegment.assign( myVals, n );  // copy from ptr to segment

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
            SCAI_ASSERT_EQ_ERROR( n, writtenSize, "notify mismatch" )
        }

        dstSegment.copyTo( allVals, n * mSize );
    }
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allVals, n, root, myVals, sizes )                                     */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::gatherVImpl(
    void* allVals,
    const IndexType n,
    const PartitionId root,
    const void* myVals,
    const IndexType sizes[],
    const common::scalar::ScalarType stype ) const
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

    SCAI_LOG_DEBUG( logger, *this << ": gatherV<" << stype << ">, n = "
                    << n << ", total = " << totalSize )

    // Note: segment data is now safe against different values of n, totalSize

    SegmentData srcSegment( stype, this, n );
    SegmentData dstSegment( stype, this, totalSize );   // allValues on root

    srcSegment.assign( myVals, n ); // copy from ptr to segment

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
                SCAI_ASSERT_EQ_ERROR( size, writtenSize, "notify mismatch" )
            }
        }

        dstSegment.copyTo( allVals, totalSize );
    }

    wait();
}

/* ---------------------------------------------------------------------------------- */
/*           maxloc                                                                   */
/* ---------------------------------------------------------------------------------- */

bool GPICommunicator::supportsLocReduction( common::scalar::ScalarType, common::scalar::ScalarType ) const
{
    // GPI does not support minloc/maxloc, so use default implementation via gather 

    return false;
}

void GPICommunicator::maxlocImpl( void*, IndexType*, PartitionId, common::scalar::ScalarType ) const
{
    COMMON_THROWEXCEPTION( "maxlocRedcution unsupported, should use default" )
}

void GPICommunicator::minlocImpl( void*, IndexType*, PartitionId, common::scalar::ScalarType ) const
{
    COMMON_THROWEXCEPTION( "minlocRedcution unsupported, should use default" )
}

/* ---------------------------------------------------------------------------------- */
/*           swap                                                                     */
/* ---------------------------------------------------------------------------------- */

void GPICommunicator::swapImpl( void* val, const IndexType n, PartitionId partner, common::scalar::ScalarType stype ) const
{
    // this is now safe as there is no global synchronization

    if ( partner == mRank )
    {
        return;
    }

    SegmentData srcSegment( stype, this, n );   // used remote, will be read
    SegmentData dstSegment( stype, this, n );   // used locally

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
    SCAI_ASSERT_EQ_ERROR( n, nn, "notify mismatch" );
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
