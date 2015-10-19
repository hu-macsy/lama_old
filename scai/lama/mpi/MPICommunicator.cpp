/**
 * @file MPICommunicator.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief MPICommunicator.cpp
 * @author Jiri Kraus
 * @date 23.02.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/mpi/MPICommunicator.hpp>

// local library
#include <scai/lama/mpi/MPISyncToken.hpp>
#include <scai/lama/mpi/MPIUtils.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/Assert.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/unique_ptr.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <iostream>
#include <algorithm>

using namespace std;

namespace scai
{

using common::shared_ptr;
using common::unique_ptr;
using common::scoped_array;

namespace lama
{

const int MPICommunicator::defaultTag = 1;
MPI_Op MPICommunicator::mSumComplexLongDouble = 0;

SCAI_LOG_DEF_LOGGER( MPICommunicator::logger, "Communicator.MPICommunicator" )

MPICommunicator::MPICommunicator( int& argc, char** & argv, const std::string& type )
    : CRTPCommunicator<MPICommunicator>( type ), 
      mMainThread( common::Thread::getSelf() ),
      mThreadSafetyLevel( Communicator::Funneled )
{
    SCAI_LOG_DEBUG( logger, "Communicator constructed, type = " << type )
    initialize( argc, argv );
}

MPICommunicator::MPICommunicator()
    : CRTPCommunicator<MPICommunicator>( "MPI" ), 
      mMainThread( common::Thread::getSelf() ),
      mThreadSafetyLevel( Communicator::Funneled )
{
    int argc = 0;
    char** argv = NULL;

    SCAI_LOG_DEBUG( logger, "MPICommunicator constructed, no args" )
    initialize( argc, argv );
}

MPICommunicator::MPICommunicator( int& argc, char** & argv )
    : CRTPCommunicator<MPICommunicator>( "MPI" ), 
      mMainThread( common::Thread::getSelf() ),
      mThreadSafetyLevel( Communicator::Funneled )
{
    SCAI_TRACE_SCOPE( false ) // switch off tracing in this scope as it might call this constructor again

    initialize( argc, argv );
}

void MPICommunicator::initialize( int& argc, char** & argv )
{
    int initialized = 0;

    LAMA_MPICALL( logger, MPI_Initialized( &initialized ), "MPI_Initialized" )

    if( initialized )
    {
        mExternInitialization = true;
        SCAI_LOG_WARN( logger, "MPI_Init: MPI has already been initialized." )
    }
    else
    {
        mExternInitialization = false;

        std::string threadSafetyEnvironment;

        const char* envConfig = getenv( "LAMA_MPICOMM_THREAD_SAFETY" );

        if( envConfig != NULL )
        {
            threadSafetyEnvironment = envConfig;
        }

        int requiredThreadSafety = MPI_THREAD_FUNNELED;

        if( threadSafetyEnvironment == "Multiple" )
        {
            requiredThreadSafety = MPI_THREAD_MULTIPLE;
        }
        else if( threadSafetyEnvironment == "Serialized" )
        {
            requiredThreadSafety = MPI_THREAD_SERIALIZED;
        }
        else if( threadSafetyEnvironment != "" )
        {
            SCAI_LOG_ERROR( logger,
                            "LAMA_MPICOMM_THREAD_SAFETY = " << threadSafetyEnvironment << ", unknown value (try Multiple or Serialized)" )
        }

        int providedThreadSafety = MPI_THREAD_SINGLE;

        LAMA_MPICALL( logger, MPI_Init_thread( &argc, &argv, requiredThreadSafety, &providedThreadSafety ), "MPI_Init" );

        switch( providedThreadSafety )
        {
            case MPI_THREAD_MULTIPLE:
                mThreadSafetyLevel = Communicator::Multiple;
                SCAI_LOG_INFO( logger, "MPI Thread Safety Level: MPI_THREAD_MULTIPLE" )
                break;

            case MPI_THREAD_SERIALIZED:
                mThreadSafetyLevel = Communicator::Serialized;
                SCAI_LOG_INFO( logger, "MPI Thread Safety Level: MPI_THREAD_SERIALIZED" )
                break;

            case MPI_THREAD_FUNNELED:
                mThreadSafetyLevel = Communicator::Funneled;
                SCAI_LOG_INFO( logger, "MPI Thread Safety Level: MPI_THREAD_FUNNELED" )
                break;

            case MPI_THREAD_SINGLE:
                mThreadSafetyLevel = Communicator::Funneled;
                SCAI_LOG_INFO( logger, "MPI Thread Safety Level: MPI_THREAD_SINGLE" )
                break;

            default:
                MPI_Finalize();
                COMMON_THROWEXCEPTION( "MPI which supports at leas thread level MPI_THREAD_FUNNELED is required." )
        }

        isCUDAAware = false;

        bool setCUDA = common::Settings::getEnvironment( isCUDAAware, "LAMA_MPI_CUDA" );

        if( setCUDA )
        {
            SCAI_LOG_ERROR( logger, "MPI isCUDAAware = " << isCUDAAware )
        }
    }

	if( mSumComplexLongDouble == 0)
	{
		MPI_Op_create( &sum_complex_long_double, true, &mSumComplexLongDouble );
		SCAI_LOG_DEBUG( logger, "MPI_Op_create for sum complex long double")
	}

    mCommWorld = MPI_COMM_WORLD;
    MPI_Comm_dup( mCommWorld, &mComm );
    MPI_Comm_dup( mCommWorld, &mCommTask );
    SCAI_LOG_INFO( logger, "MPI_Init" )
    LAMA_MPICALL( logger, MPI_Comm_size( mComm, &mSize ), "MPI_Comm_size" )
    LAMA_MPICALL( logger, MPI_Comm_rank( mComm, &mRank ), "MPI_Comm_rank" )

    setNodeData(); // determine mNodeRank, mNodeSize

    // set rank, output string in an environment variable 
    // so it might be used by logging, tracing, etc.

    std::ostringstream commVal;

    commVal << *this;

    common::Settings::putEnvironment( "SCAI_COMM", commVal.str().c_str() );
    common::Settings::putEnvironment( "SCAI_RANK", mRank );
}

void MPICommunicator::sum_complex_long_double(void *in, void *out, int *count,
                                 MPI_Datatype * UNUSED(dtype) )
{
  ComplexLongDouble *a = reinterpret_cast<ComplexLongDouble*>( in );
  ComplexLongDouble *b = reinterpret_cast<ComplexLongDouble*>( out );
  for(int i = 0; i < *count; ++i) {
      b[i] += a[i];
  }
}


/* ---------------------------------------------------------------------------------- */

void MPICommunicator::setNodeData()
{
    // routine set mNodeRank and mNodeSize

    // processors with same processor_name are assumed to be on the same node

    int nodeNameLength; // lenght of node name for this processor

    char nodeName[MPI_MAX_PROCESSOR_NAME]; // name of node for this processor
    memset( nodeName, '\0', MPI_MAX_PROCESSOR_NAME );

    LAMA_MPICALL( logger, MPI_Get_processor_name( nodeName, &nodeNameLength ), "MPI_Get_processor_name" )

    SCAI_LOG_INFO( logger, "Processor " << mRank << " runs on node " << nodeName )

    char* allNodeNames = (char*) malloc( MPI_MAX_PROCESSOR_NAME * mSize * sizeof(char) );

    if( allNodeNames == NULL )
    {
        COMMON_THROWEXCEPTION( "Can't alloc enough memory for all node names." )
    }

    memset( allNodeNames, '\0', MPI_MAX_PROCESSOR_NAME * mSize );

    LAMA_MPICALL( logger,
                  MPI_Allgather( &nodeName[0], MPI_MAX_PROCESSOR_NAME, MPI_CHAR, &allNodeNames[0],
                                 MPI_MAX_PROCESSOR_NAME, MPI_CHAR, mComm ),
                  "MPI_Allgather( <node_names> )" )

    mNodeSize = 0;
    mNodeRank = mSize; // illegal value to verify that it will be set

    for( int i = 0; i < mSize; ++i )
    {
        if( strcmp( &allNodeNames[i * MPI_MAX_PROCESSOR_NAME], nodeName ) )
        {
            continue; // processor i is not on same node
        }

        // Processor i is on same node as this processor

        if( i == mRank )
        {
            mNodeRank = mNodeSize;
        }

        ++mNodeSize;
    }

    free( allNodeNames );

    SCAI_ASSERT_ERROR( mNodeSize > 0, "Serious problem encountered to get node size" )

    SCAI_ASSERT_ERROR( mNodeRank < mNodeSize, "Serious problem encountered to get node size" )

    SCAI_LOG_INFO( logger, "Processor " << mRank << ": node rank " << mNodeRank << " of " << mNodeSize )
}

/* ---------------------------------------------------------------------------------- */

MPICommunicator::~MPICommunicator()
{
    SCAI_LOG_INFO( logger, *this << ": ~MPICommunicator" )
    int finalized = 0;
    LAMA_MPICALL( logger, MPI_Finalized( &finalized ), "MPI_Finalized" )

    if( !finalized )
    {
        if( !mExternInitialization )
        {
            SCAI_LOG_INFO( logger, "call MPI_Finalize" )
            LAMA_MPICALL( logger, MPI_Finalize(), "MPI_Finalize" )
        }
        else
        {
            SCAI_LOG_INFO( logger, "ATTENTION: no call MPI_Finalize, was externally initialized" )
        }
    }
    else
    {
        SCAI_LOG_WARN( logger, *this << ": tried to finalize MPI, but MPI has already been finalized." )
    }
}

/* ---------------------------------------------------------------------------------- */

bool MPICommunicator::isEqual( const Communicator& other ) const
{
    bool equal = false;
    const MPICommunicator* otherMPI = dynamic_cast<const MPICommunicator*>( &other );

    if( otherMPI )
    {
        equal = mComm == otherMPI->mComm;
    }

    return equal;
}

/* ---------------------------------------------------------------------------------- */

Communicator::ThreadSafetyLevel MPICommunicator::getThreadSafetyLevel() const
{
    return mThreadSafetyLevel;
}

PartitionId MPICommunicator::getSize() const
{
    return mSize;
}

PartitionId MPICommunicator::getRank() const
{
    return mRank;
}

PartitionId MPICommunicator::getNodeSize() const
{
    return mNodeSize;
}

PartitionId MPICommunicator::getNodeRank() const
{
    return mNodeRank;
}

MPI_Comm MPICommunicator::getMPIComm() const
{
    return mComm;
}

template<typename ValueType>
MPI_Request MPICommunicator::startrecv( ValueType* buffer, int count, int source ) const
{
    MPI_Request request;
    MPI_Datatype commType = getMPIType<ValueType>();
    LAMA_MPICALL( logger, MPI_Irecv( buffer, count, commType, source, defaultTag, selectMPIComm(), &request ),
                  "MPI_Irecv" )
    return request;
}

template<typename ValueType>
MPI_Request MPICommunicator::startsend( const ValueType* buffer, int count, int target ) const
{
    MPI_Request request;
    MPI_Datatype commType = getMPIType<ValueType>();
    LAMA_MPICALL( logger,
                  MPI_Isend( const_cast<ValueType*>( buffer ), count, commType, target, defaultTag, selectMPIComm(),
                             &request ),
                  "MPI_Isend" )
    return request;
}

template<typename ValueType>
int MPICommunicator::getCount( MPI_Status& mpiStatus ) const
{
    int size = 0;
    MPI_Datatype commType = getMPIType<ValueType>();
    LAMA_MPICALL( logger, MPI_Get_count( &mpiStatus, commType, &size ), "MPI_Get_count" )
    return size;
}

template<typename ValueType>
void MPICommunicator::send( const ValueType buffer[], int count, int target ) const
{
    MPI_Datatype commType = getMPIType<ValueType>();
    LAMA_MPICALL( logger,
                  MPI_Send( const_cast<ValueType*>( buffer ), count, commType, target, defaultTag, selectMPIComm() ),
                  "MPI_Send" )
}

/* ---------------------------------------------------------------------------------- */

void MPICommunicator::all2all( IndexType recvSizes[], const IndexType sendSizes[] ) const
{
    SCAI_REGION( "Communicator.MPI.all2all" )

    SCAI_ASSERT_ERROR( sendSizes != 0, " invalid sendSizes " )
    SCAI_ASSERT_ERROR( recvSizes != 0, " invalid recvSizes " )

    // MPI is not const-aware so we have to use a const_cast on sendSizes

    MPI_Datatype commType = getMPIType<IndexType>();

    LAMA_MPICALL( logger,
                  MPI_Alltoall( const_cast<IndexType*>( sendSizes ), 1, commType, recvSizes,
                                1, commType, selectMPIComm() ),
                  "MPI_Alltoall" )
}

/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void MPICommunicator::exchangeByPlanImpl(
    ValueType recvData[],
    const CommunicationPlan& recvPlan,
    const ValueType sendData[],
    const CommunicationPlan& sendPlan ) const
{
    SCAI_REGION( "Communicator.MPI.exchangeByPlan" )

    SCAI_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    SCAI_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )

    SCAI_LOG_INFO( logger,
                   *this << ": exchange for values of type " << common::getScalarType<ValueType>() 
                   << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )

    int maxReceives = recvPlan.size();
    int noReceives = 0; // will be incremented
    ValueType* recvDataForMe = NULL;
    IndexType recvDataForMeSize = 0;
    scoped_array<MPI_Request> commRequest( new MPI_Request[maxReceives] );

    // setup receives for each entry in receive plan

    for( PartitionId i = 0; i < maxReceives; ++i )
    {
        IndexType quantity = recvPlan[i].quantity;
        IndexType offset = recvPlan[i].offset;
        ValueType* recvDataForI = recvData + offset;
        PartitionId p = recvPlan[i].partitionId;
        SCAI_LOG_DEBUG( logger,
                        *this << ": receive " << quantity << " elements" << " from processor " << p << " at offset " << offset )

        if( p != mRank )
        {
            commRequest[noReceives] = startrecv( recvDataForI, quantity, p );
            noReceives++;
        }
        else
        {
            recvDataForMe = recvDataForI;
            recvDataForMeSize = quantity;
        }
    }

    // send the data via sendPlan

    for( PartitionId i = 0; i < sendPlan.size(); ++i )
    {
        IndexType quantity = sendPlan[i].quantity;
        IndexType offset = sendPlan[i].offset;
        const ValueType* sendDataForI = sendData + offset;
        PartitionId p = sendPlan[i].partitionId;
        SCAI_LOG_DEBUG( logger,
                        *this << ": send " << quantity << " elements" << " to processor " << p << " at offset " << offset )

        if( p != mRank )
        {
            send( sendDataForI, quantity, p );
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "self-exchange of " << quantity << " elements" )
            SCAI_ASSERT_DEBUG( quantity == recvDataForMeSize, "size mismatch for self exchange" )

            for( IndexType k = 0; k < recvDataForMeSize; k++ )
            {
                recvDataForMe[k] = sendDataForI[k];
            }
        }
    }

    // wait for completion of receives
    scoped_array<MPI_Status> statuses( new MPI_Status[noReceives] );
    LAMA_MPICALL( logger, MPI_Waitall( noReceives, commRequest.get(), statuses.get() ), "MPI_Waitall" )
    // ToDo: check for correct sizes, was done in earlier version, but is now redundant
}

/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* MPICommunicator::exchangeByPlanAsyncImpl(
    ValueType* const recvData,
    const CommunicationPlan& recvPlan,
    const ValueType* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    SCAI_REGION( "Communicator.MPI.exchangeByPlanAsync" )

    SCAI_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    SCAI_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )
    SCAI_LOG_INFO( logger,
                   *this << ": exchange for values of type " << common::getScalarType<ValueType>() 
                    << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )
    int noRequests = sendPlan.size() + recvPlan.size();

    // create MPIToken as auto_ptr, so it will be freed in case of exception

    scai::common::unique_ptr<MPISyncToken> pSyncToken( new MPISyncToken( noRequests ) );

    MPISyncToken& syncToken = *pSyncToken;

    ValueType* recvDataForMe = NULL;
    IndexType recvDataForMeSize = 0;

    // setup receives for each entry in receive plan

    for( PartitionId i = 0; i < recvPlan.size(); ++i )
    {
        IndexType quantity = recvPlan[i].quantity;
        ValueType* recvDataForI = recvData + recvPlan[i].offset;
        PartitionId p = recvPlan[i].partitionId;
        SCAI_LOG_DEBUG( logger, *this << ": receive " << quantity << " elements" << " from processor " << p )

        if( p != mRank )
        {
            syncToken.pushRequest( startrecv( recvDataForI, quantity, p ) );
        }
        else
        {
            recvDataForMe = recvDataForI;
            recvDataForMeSize = quantity;
        }
    }

    // send the data via sendPlan

    for( PartitionId i = 0; i < sendPlan.size(); ++i )
    {
        IndexType quantity = sendPlan[i].quantity;
        const ValueType* sendDataForI = sendData + sendPlan[i].offset;
        PartitionId p = sendPlan[i].partitionId;
        SCAI_LOG_DEBUG( logger, *this << ": send " << quantity << " elements" << " to processor " << p )

        if( p != mRank )
        {
            syncToken.pushRequest( startsend( sendDataForI, quantity, p ) );
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "self-exchange of " << quantity << " elements" )
            SCAI_ASSERT_DEBUG( quantity == recvDataForMeSize, "size mismatch for self exchange" )

            for( IndexType k = 0; k < recvDataForMeSize; k++ )
            {
                recvDataForMe[k] = sendDataForI[k];
            }
        }
    }

    return pSyncToken.release();
}

/* ---------------------------------------------------------------------------------- */

inline MPI_Comm MPICommunicator::selectMPIComm() const
{
    if( mThreadSafetyLevel != Multiple )
    {
        return mComm;
    }

    const common::Thread::Id thisThread = common::Thread::getSelf();

    if( thisThread == mMainThread )
    {
        return mComm;
    }
    else
    {
        return mCommTask;
    }
}

/* ---------------------------------------------------------------------------------- */
/*              bcast                                                                 */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void MPICommunicator::bcastImpl( ValueType val[], const IndexType n, const PartitionId root ) const
{
    SCAI_REGION( "Communicator.MPI.bcast" )

    MPI_Datatype commType = getMPIType<ValueType>();
    LAMA_MPICALL( logger, MPI_Bcast( val, n, commType, root, selectMPIComm() ), "MPI_Bcast<ValueType>" )
}

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
IndexType MPICommunicator::shiftImpl(
    ValueType recvVals[],
    const IndexType recvSize,
    const PartitionId source,
    const ValueType sendVals[],
    const IndexType sendSize,
    const PartitionId dest ) const
{
    SCAI_REGION( "Communicator.MPI.shift" )

    SCAI_ASSERT_ERROR( source != getRank(), "source must not be this partition" )
    SCAI_ASSERT_ERROR( dest != getRank(), "dest must not be this partition" )

    SCAI_LOG_DEBUG( logger,
                    *this << ": recv from " << source << " max " << recvSize << " values " << ", send to " << dest << " " << sendSize << " values." )

    MPI_Datatype commType = getMPIType<ValueType>();
    MPI_Status mpiStatus;

    LAMA_MPICALL( logger,
                  MPI_Sendrecv( const_cast<ValueType*>( sendVals ), sendSize, commType, dest, 4711, recvVals, recvSize,
                                commType, source, 4711, selectMPIComm(), &mpiStatus ),
                  "MPI_Sendrecv" )

    // extract number of read values from status

    int count = 0;

    LAMA_MPICALL( logger, MPI_Get_count( &mpiStatus, commType, &count ), "MPI_Get_count(ValueType)" )

    SCAI_LOG_DEBUG( logger, "received from " << source << " #values = " << count << ", max was " << recvSize )

    return count;
}

/* ---------------------------------------------------------------------------------- */
/*              shiftAsync                                                            */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* MPICommunicator::shiftAsyncImpl(
    ValueType recvVals[],
    const PartitionId source,
    const ValueType sendVals[],
    const PartitionId dest,
    const IndexType size ) const
{
    SCAI_LOG_DEBUG( logger,
                    *this << ": recv from " << source << ", send to " << dest << ", both " << size << " values." )

    SCAI_ASSERT_ERROR( source != getRank(), "source must not be this partition" )
    SCAI_ASSERT_ERROR( dest != getRank(), "dest must not be this partition" )

    // need an MPI communicator with 2 requests, no clean up needed

    unique_ptr<MPISyncToken> pSyncToken( new MPISyncToken( 2 ) );

    pSyncToken->pushRequest( startrecv( recvVals, size, source ) );
    pSyncToken->pushRequest( startsend( sendVals, size, dest ) );

    return pSyncToken.release();
}

/* ---------------------------------------------------------------------------------- */
/*              sum                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MPICommunicator::sumImpl( const ValueType value ) const
{
    SCAI_REGION( "Communicator.MPI.sum" )

    ValueType sum;
    MPI_Datatype commType = getMPIType<ValueType>();
    MPI_Op opType = getMPISum<ValueType>();
    LAMA_MPICALL( logger, MPI_Allreduce( (void* ) &value, (void* ) &sum, 1, commType, opType,
                  selectMPIComm() ), "MPI_Allreduce(MPI_SUM)" )
    SCAI_LOG_DEBUG( logger, "sum: my value = " << value << ", sum = " << sum )
    return sum;
}

/* ---------------------------------------------------------------------------------- */
/*              min / max reduction                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MPICommunicator::minImpl( const ValueType value ) const
{
    SCAI_REGION( "Communicator.MPI.min" )

    MPI_Datatype commType = getMPIType<ValueType>();

    ValueType globalMin; // no initialization needed, done in MPI call

    LAMA_MPICALL( logger, MPI_Allreduce( (void* ) &value, (void* ) &globalMin, 1, commType,
                                         MPI_MIN, selectMPIComm() ), "MPI_Allreduce( MPI_MIN )" )
    return globalMin;
}

template<typename ValueType>
ValueType MPICommunicator::maxImpl( const ValueType value ) const
{
    SCAI_REGION( "Communicator.MPI.max" )

    MPI_Datatype commType = getMPIType<ValueType>();

    ValueType globalMax; // no initialization needed, done in MPI call

    SCAI_LOG_DEBUG( logger, "maxImpl: local value = " << value )

    LAMA_MPICALL( logger, MPI_Allreduce( (void* ) &value, (void* ) &globalMax, 1, commType, MPI_MAX,
                                         selectMPIComm() ), "MPI_Allreduce( MPI_MAX )" )

    SCAI_LOG_DEBUG( logger, "maxImpl: global value = " << globalMax )

    return globalMax;
}

void MPICommunicator::synchronize() const
{
    LAMA_MPICALL( logger, MPI_Barrier( selectMPIComm() ), "MPI_Barrier()" )
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myvals, n, root, allvals )                                           */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void MPICommunicator::scatterImpl(
    ValueType myvals[],
    const IndexType n,
    const PartitionId root,
    const ValueType allvals[] ) const
{
    SCAI_REGION( "Communicator.MPI.scatter" )

    SCAI_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
    SCAI_LOG_DEBUG( logger, *this << ": scatter of " << n << " elements, root = " << root )
    MPI_Datatype commType = getMPIType<ValueType>();
    // MPI interface is not aware of const, so const_cast is required
    LAMA_MPICALL( logger,
                  MPI_Scatter( const_cast<ValueType*>( allvals ), n, commType, myvals, n, commType, root,
                               selectMPIComm() ),
                  "MPI_Scatter" )
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myvals, n, root, allvals, sizes )                                    */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void MPICommunicator::scatterVImpl(
    ValueType myvals[],
    const IndexType n,
    const PartitionId root,
    const ValueType allvals[],
    const IndexType sizes[] ) const
{
    SCAI_REGION( "Communicator.MPI.scatterV" )

    SCAI_ASSERT_ERROR( root < getSize(), "illegal root, root = " << root )
    MPI_Datatype commType = getMPIType<ValueType>();

    if( root == getRank() )
    {
        void* sendbuf = const_cast<ValueType*>( allvals );
        PartitionId np = getSize();
        scoped_array<int> counts( new int[np] );
        scoped_array<int> displs( new int[np] );
        int displacement = 0;

        for( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = static_cast<int>( sizes[i] );
            displs[i] = displacement;
            displacement += counts[i];
        }

        SCAI_LOG_DEBUG( logger,
                        *this << ": scatter of " << displacement << " elements, I receive " << n << " elements" )
        LAMA_MPICALL( logger,
                      MPI_Scatterv( sendbuf, counts.get(), displs.get(), commType, myvals, n, commType, root,
                                    selectMPIComm() ),
                      "MPI_Scatterv" )
    }
    else
    {
        // VampirTrace: requires valid counts array, even if values will be ignored

        PartitionId np = getSize();
        scoped_array<int> counts( new int[np] );

        for( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = 0;
        }

        SCAI_LOG_DEBUG( logger, *this << ": root = " << root << " scatters " << n << " elements to me" )
        LAMA_MPICALL( logger,
                      MPI_Scatterv( NULL, counts.get(), NULL, commType, myvals, n, commType, root, selectMPIComm() ),
                      "MPI_Scatterv" )
    }
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allvals, n, root, myvals )                                            */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void MPICommunicator::gatherImpl(
    ValueType allvals[],
    const IndexType n,
    const PartitionId root,
    const ValueType myvals[] ) const
{
    SCAI_REGION( "Communicator.MPI.gather" )

    SCAI_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
    SCAI_LOG_DEBUG( logger, *this << ": gather of " << n << " elements, root = " << root )
    MPI_Datatype commType = getMPIType<ValueType>();
    // MPI interface is not aware of const, so const_cast is required
    void* sendbuf = const_cast<ValueType*>( myvals );
    LAMA_MPICALL( logger, MPI_Gather( sendbuf, n, commType, allvals, n, commType, root, selectMPIComm() ),
                  "MPI_Gather<ValueType>" )
}

/* ---------------------------------------------------------------------------------- */
/*      gatherV( allvals, n, root, myvals, sizes )                                    */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void MPICommunicator::gatherVImpl(
    ValueType allvals[],
    const IndexType n,
    const PartitionId root,
    const ValueType myvals[],
    const IndexType sizes[] ) const
{
    SCAI_REGION( "Communicator.MPI.gatherV" )

    SCAI_ASSERT_ERROR( root < getSize(), "illegal root, root = " << root )
    void* sendbuf = const_cast<ValueType*>( myvals );
    MPI_Datatype commType = getMPIType<ValueType>();

    if( root == getRank() )
    {
        PartitionId np = getSize();
        scoped_array<int> counts( new int[np] );
        scoped_array<int> displs( new int[np] );
        int displacement = 0;

        for( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = static_cast<int>( sizes[i] );
            displs[i] = displacement;
            displacement += counts[i];
        }

        SCAI_LOG_DEBUG( logger,
                        *this << ": scatter of " << displacement << " elements, I receive " << n << " elements" )
        LAMA_MPICALL( logger,
                      MPI_Gatherv( sendbuf, n, commType, allvals, counts.get(), displs.get(), commType, root,
                                   selectMPIComm() ),
                      "MPI_Gatherv<ValueType>" )
    }
    else
    {
        // VampirTrace: requires valid counts array, even if values will be ignored

        PartitionId np = getSize();
        scoped_array<int> counts( new int[np] );

        for( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = 0;
        }

        SCAI_LOG_DEBUG( logger, *this << ": root = " << root << " scatters " << n << " elements to me" )
        LAMA_MPICALL( logger,
                      MPI_Gatherv( sendbuf, n, commType, NULL, counts.get(), NULL, commType, root, selectMPIComm() ),
                      "MPI_Gatherv<ValueType>" )
    }
}

/* ---------------------------------------------------------------------------------- */
/*           maxloc                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void MPICommunicator::maxlocImpl( ValueType& val, IndexType& location, PartitionId root ) const
{
    SCAI_REGION( "Communicator.MPI.maxloc" )

    struct ValAndLoc
    {
        ValueType val;
        int location;
    };
    ValAndLoc in;
    in.val = val;
    in.location = location;
    ValAndLoc out;
    MPI_Datatype commType = getMPI2Type<ValueType,int>();
    MPI_Reduce( &in, &out, 1, commType, MPI_MAXLOC, root, selectMPIComm() );

    if( mRank == root )
    {
        val = out.val;
        location = out.location;
    }
}

/* ---------------------------------------------------------------------------------- */
/*           swap                                                                     */
/* ---------------------------------------------------------------------------------- */

template<typename ValueType>
void MPICommunicator::swapImpl( ValueType val[], const IndexType n, PartitionId partner ) const
{
    SCAI_REGION( "Communicator.MPI.swap" )

    if( partner == mRank )
    {
        return;
    }

    scoped_array<ValueType> tmp( new ValueType[n] );

    for( IndexType i = 0; i < n; i++ )
    {
        tmp[i] = val[i];
    }

    MPI_Status mpiStatus;
    MPI_Datatype commType = getMPIType<ValueType>();
    LAMA_MPICALL( logger,
                  MPI_Sendrecv( tmp.get(), n, commType, partner, defaultTag, val, n, commType, partner, defaultTag,
                                selectMPIComm(), &mpiStatus ),
                  "MPI_Sendrecv" )
    SCAI_ASSERT_ERROR( getCount<ValueType>( mpiStatus ) == n, "size mismatch for swap" )
}

hmemo::ContextPtr MPICommunicator::getCommunicationContext( const hmemo::ContextArray& array ) const
{
    // get a valid context, i.e. a context that contains valid data

    hmemo::ContextPtr validContext = array.getValidContext( hmemo::context::Host );

    SCAI_LOG_DEBUG( logger, "CommunicationContext: valid context for " << array << ": " << *validContext )

    if ( validContext->getType() == hmemo::context::Host )
    {
        return validContext;
    }

    // This can only be used for CUDAaware MPI

    if( isCUDAAware && ( validContext->getType() == hmemo::context::CUDA ) )
    {
        return validContext;
    }

    return hmemo::Context::getContextPtr( hmemo::context::Host );
}

/* ---------------------------------------------------------------------------------- */
/*           writeAt                                                                  */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::writeAt( std::ostream& stream ) const
{
    // Info about rank and size of the communicator is very useful
    stream << "MPI(" << mRank << ":" << mSize << ")";
}

// template instantiation for the supported array types

#define LAMA_MPI_METHODS_INSTANTIATE(z, I, _)                         \
    template COMMON_DLL_IMPORTEXPORT                                    \
    void MPICommunicator::maxlocImpl(                                 \
            ARRAY_TYPE##I &, IndexType&, PartitionId) const;

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_MPI_METHODS_INSTANTIATE, _ )

#undef LAMA_MPI_METHODS_INSTANTIATE

/* --------------------------------------------------------------- */

static shared_ptr<const MPICommunicator> theMPICommunicator;

MPICommunicator::MPIGuard::MPIGuard()
{
}

MPICommunicator::MPIGuard::~MPIGuard()
{
    if ( theMPICommunicator.use_count() > 1 )
    {
        SCAI_LOG_WARN( logger,
                   "MPICommunicator has " << theMPICommunicator.use_count() - 1 << " remaining references, seems that not all LAMA data structures have been freed" )
    }
}

// create a guard whose destructor at program exit takes care of MPI exit call

// MPICommunicator::MPIGuard MPICommunicator::guard = 

CommunicatorPtr MPICommunicator::create()
{
    if( !theMPICommunicator )
    {
        SCAI_LOG_INFO( logger, "create new MPICommunicator" )

        // create a new instance of MPICommunicator, will call MPI_Init

        theMPICommunicator = shared_ptr<const MPICommunicator>( new MPICommunicator() );
    }

    return theMPICommunicator;
}

/* --------------------------------------------------------------- */

std::string MPICommunicator::createValue()
{
    return "MPI";
}

} /* end namespace lama */

} /* end namespace scai */
