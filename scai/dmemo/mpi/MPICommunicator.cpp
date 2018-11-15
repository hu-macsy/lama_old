/**
 * @file MPICommunicator.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief MPICommunicator.cpp
 * @author Jiri Kraus
 * @date 23.02.2011
 */

// hpp
#include <scai/dmemo/mpi/MPICommunicator.hpp>

// local library
#include <scai/dmemo/mpi/MPISyncToken.hpp>
#include <scai/dmemo/mpi/MPIUtils.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/loop.hpp>
#include <scai/common/BinaryOp.hpp>

// std
#include <iostream>
#include <algorithm>
#include <memory>

using namespace std;

using std::unique_ptr;

namespace scai
{

namespace dmemo
{

const int MPICommunicator::defaultTag = 1;

#ifdef SCAI_COMPLEX_SUPPORTED
MPI_Datatype MPICommunicator::mComplexLongDoubleType = 0;

MPI_Op MPICommunicator::mSumComplexLongDouble = 0;

#endif

SCAI_LOG_DEF_LOGGER( MPICommunicator::logger, "Communicator.MPICommunicator" )

MPICommunicator::MPICommunicator( int& argc, char** & argv, const CommunicatorKind& type )
    : Communicator( type ),
      mMainThread( std::this_thread::get_id() ),
      mThreadSafetyLevel( Communicator::Funneled )
{
    SCAI_LOG_DEBUG( logger, "Communicator constructed, type = " << type )
    initialize( argc, argv );
}

MPICommunicator::MPICommunicator()
    : Communicator( MPI ),
      mMainThread( std::this_thread::get_id() ),
      mThreadSafetyLevel( Communicator::Funneled )
{
    int argc = 0;
    char** argv = NULL;
    SCAI_LOG_DEBUG( logger, "MPICommunicator constructed, no args" )
    initialize( argc, argv );
}

MPICommunicator::MPICommunicator( int& argc, char** & argv )
    : Communicator( MPI ),
      mMainThread( std::this_thread::get_id() ),
      mThreadSafetyLevel( Communicator::Funneled )
{
    SCAI_TRACE_SCOPE( false ) // switch off tracing in this scope as it might call this constructor again
    initialize( argc, argv );
}

void MPICommunicator::initialize( int& argc, char** & argv )
{
    int initialized = 0;
    SCAI_MPICALL( logger, MPI_Initialized( &initialized ), "MPI_Initialized" )

    if ( initialized )
    {
        mKind = MPICommKind::EXTERNAL;
        SCAI_LOG_WARN( logger, "MPI_Init: MPI has already been initialized." )
    }
    else
    {
        mKind = MPICommKind::INTERNAL;
        std::string threadSafetyEnvironment;
        const char* envConfig = getenv( "LAMA_MPICOMM_THREAD_SAFETY" );

        if ( envConfig != NULL )
        {
            threadSafetyEnvironment = envConfig;
        }

        int requiredThreadSafety = MPI_THREAD_FUNNELED;

        if ( threadSafetyEnvironment == "Multiple" )
        {
            requiredThreadSafety = MPI_THREAD_MULTIPLE;
        }
        else if ( threadSafetyEnvironment == "Serialized" )
        {
            requiredThreadSafety = MPI_THREAD_SERIALIZED;
        }
        else if ( threadSafetyEnvironment != "" )
        {
            SCAI_LOG_ERROR( logger,
                            "LAMA_MPICOMM_THREAD_SAFETY = " << threadSafetyEnvironment << ", unknown value (try Multiple or Serialized)" )
        }

        int providedThreadSafety = MPI_THREAD_SINGLE;
        SCAI_MPICALL( logger, MPI_Init_thread( &argc, &argv, requiredThreadSafety, &providedThreadSafety ), "MPI_Init" );

        switch ( providedThreadSafety )
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
        bool setCUDA = common::Settings::getEnvironment( isCUDAAware, "SCAI_MPI_CUDA" );

        if ( setCUDA )
        {
            SCAI_LOG_INFO( logger, "MPI isCUDAAware = " << isCUDAAware )
        }
    }

#ifdef SCAI_COMPLEX_SUPPORTED

    if ( mSumComplexLongDouble == 0 )
    {
        SCAI_MPICALL( logger, MPI_Op_create( &sum_complex_long_double, true, &mSumComplexLongDouble ), "Create MPI_Op sum for ComplexLongDouble" )
        SCAI_LOG_DEBUG( logger, "MPI_Op_create for sum complex long double" )
    }

    if ( mComplexLongDoubleType == 0 )
    {
        SCAI_MPICALL( logger, MPI_Type_contiguous( 2, MPI_LONG_DOUBLE, &mComplexLongDoubleType ), "Call init of MPI_Datatype for ComplexLongDouble" )
        SCAI_MPICALL( logger, MPI_Type_commit( &mComplexLongDoubleType ), "Register ComplexLongDouble MPI_Datatype" )
        SCAI_LOG_DEBUG( logger, "register MPI_Datatype for ComplexLongDouble" )
    }

#endif

    MPI_Comm_dup( MPI_COMM_WORLD, &mComm );
    MPI_Comm_set_errhandler( mComm, MPI_ERRORS_RETURN );
    SCAI_LOG_INFO( logger, "MPI_Init" )

    {
        // temporary variables required as PartitionId can be different from int

        int mpiRank;
        int mpiSize;

        SCAI_MPICALL( logger, MPI_Comm_size( mComm, &mpiSize ), "MPI_Comm_size" )
        SCAI_MPICALL( logger, MPI_Comm_rank( mComm, &mpiRank ), "MPI_Comm_rank" )

        setSizeAndRank( static_cast<PartitionId>( mpiSize ), static_cast<PartitionId>( mpiRank ) );
    }

    // set rank, output string in an environment variable
    // so it might be used by logging, tracing, etc.

    std::ostringstream commVal;
    commVal << *this;

    common::Settings::putEnvironment( "SCAI_COMM", commVal.str().c_str() );
    common::Settings::putEnvironment( "SCAI_RANK", getRank() );

    // tracing of MPI calls for getting node data can already be traced

    setNodeData(); // determine mNodeRank, mNodeSize
}

#ifdef SCAI_COMPLEX_SUPPORTED

void MPICommunicator::sum_complex_long_double( void* in, void* out, int* count,
        MPI_Datatype* SCAI_UNUSED( dtype ) )
{
    ComplexLongDouble* a = reinterpret_cast<ComplexLongDouble*>( in );
    ComplexLongDouble* b = reinterpret_cast<ComplexLongDouble*>( out );

    for ( int i = 0; i < *count; ++i )
    {
        b[i] += a[i];
    }
}

#endif

/* ---------------------------------------------------------------------------------- */


void MPICommunicator::getProcessorName( char* name ) const
{
    size_t len = maxProcessorName();

    memset( name, '\0', len * sizeof( char ) );

    int nodeNameLength;   // not really neded as terminated with \0

    SCAI_MPICALL( logger, MPI_Get_processor_name( name, &nodeNameLength ), "MPI_Get_processor_name" )

    SCAI_LOG_INFO( logger, "Processor " << getRank() << " runs on node " << name )
}

/* --------------------------------------------------------------- */

size_t MPICommunicator::maxProcessorName() const
{
    return MPI_MAX_PROCESSOR_NAME;
}

/* ---------------------------------------------------------------------------------- */

MPICommunicator::~MPICommunicator()
{
    SCAI_LOG_INFO( logger, *this << ": ~MPICommunicator" )
    int finalized = 0;

    MPI_Finalized( &finalized );

    if ( !finalized )
    {
        if ( mKind == MPICommKind::INTERNAL )
        {
#ifdef SCAI_COMPLEX_SUPPORTED
            MPI_Type_free( &mComplexLongDoubleType );
            mComplexLongDoubleType = 0;
#endif

            SCAI_LOG_INFO( logger, "call MPI_Finalize" )

            SCAI_MPICALL_NOTHROW( logger, MPI_Finalize(), "~MPICommunicator" )

        }
        else if ( mKind == MPICommKind::EXTERNAL )
        {
            SCAI_LOG_INFO( logger, "ATTENTION: no call MPI_Finalize, was externally initialized" )
        }
        else
        {
            // free the communicator

            SCAI_MPICALL_NOTHROW( logger, MPI_Comm_free( &mComm ), "~MPICommunicator" )
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

    if ( otherMPI )
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

MPI_Comm MPICommunicator::getMPIComm() const
{
    return mComm;
}

MPI_Request MPICommunicator::startrecv( void* buffer, int count, int source, common::ScalarType stype ) const
{
    MPI_Request request;

    MPI_Datatype commType = getMPIType( stype );

    SCAI_MPICALL( logger, MPI_Irecv( buffer, count, commType, source, defaultTag, mComm, &request ),
                  "MPI_Irecv<" << stype << ">" )

    return request;
}

MPI_Request MPICommunicator::startsend( const void* buffer, int count, int target, common::ScalarType stype ) const
{
    MPI_Request request;

    MPI_Datatype commType = getMPIType( stype );

    void* sBuffer = const_cast<void*>( buffer );  // MPI is not const aware

    SCAI_MPICALL( logger,
                  MPI_Isend( sBuffer, count, commType, target, defaultTag, mComm, &request ),
                  "MPI_Isend<" << stype << ">" )

    return request;
}

int MPICommunicator::getCount( MPI_Status& mpiStatus, common::ScalarType stype ) const
{
    int size = 0;
    MPI_Datatype commType = getMPIType( stype );
    SCAI_MPICALL( logger, MPI_Get_count( &mpiStatus, commType, &size ), "MPI_Get_count<" << stype << ">" )
    return size;
}

void MPICommunicator::send( const void* buffer, int count, int target, common::ScalarType stype ) const
{
    MPI_Datatype commType = getMPIType( stype );
    SCAI_MPICALL( logger,
                  MPI_Send( const_cast<void*>( buffer ), count, commType, target, defaultTag, mComm ),
                  "MPI_Send<" << stype << ">" )
}

/* ---------------------------------------------------------------------------------- */

void MPICommunicator::all2allImpl( void* recvValues, const void* sendValues, common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.all2all" )

    SCAI_ASSERT_ERROR( sendValues != NULL, "invalid sendValues" )
    SCAI_ASSERT_ERROR( recvValues != NULL, "invalid recvValues" )

    MPI_Datatype commType = getMPIType( stype );

    // MPI is not const-aware so we have to use a const_cast on sendValues

    SCAI_MPICALL( logger,
                  MPI_Alltoall( const_cast<void*>( sendValues ), 1, commType, recvValues, 1, commType, mComm ),
                  "MPI_Alltoall" )
}

/* ---------------------------------------------------------------------------------- */

/** The following method returns the memory of the actual context.
 *
 *  This memory is required for send/recv communication on the same processor
 *  where memcpy is called to do this instead of send/recv. As memcpy is not
 *  CUDA-aware we must use the routine of the corresponding device memory.
 *
 *  This is safe, as source and target are always in the same memory and not
 *  in different memories; otherwise memory transfer between different memories would
 *  be required.
 */

static inline hmemo::Memory& getActualMemory()
{
    const hmemo::Context* ctx = hmemo::Context::getCurrentContext();

    if ( ctx == NULL )
    {
        ctx = hmemo::Context::getHostPtr().get();
    }

    return *ctx->getMemoryPtr();
}

/* ---------------------------------------------------------------------------------- */

void MPICommunicator::exchangeByPlanImpl(
    void* recvData,
    const CommunicationPlan& recvPlan,
    const void* sendData,
    const CommunicationPlan& sendPlan,
    common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.exchangeByPlan" )

    SCAI_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    SCAI_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )

    SCAI_LOG_INFO( logger,
                   *this << ": exchange for values of type " << stype
                   << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )

    PartitionId maxReceives = recvPlan.size();
    PartitionId noReceives = 0; // will be incremented
    void* recvDataForMe = NULL;
    IndexType recvDataForMeSize = 0;
    unique_ptr<MPI_Request[]> commRequest( new MPI_Request[maxReceives] );

    size_t typeSize = common::typeSize( stype );

    // setup receives for each entry in receive plan

    for ( PartitionId i = 0; i < maxReceives; ++i )
    {
        IndexType quantity = recvPlan[i].quantity;
        IndexType offset = recvPlan[i].offset;
        char* recvDataForI = reinterpret_cast<char*>( recvData ) + offset * typeSize;
        PartitionId p = recvPlan[i].partitionId;
        SCAI_LOG_DEBUG( logger,
                        *this << ": receive " << quantity << " elements" << " from processor " << p << " at offset " << offset )

        if ( p != getRank() )
        {
            commRequest[noReceives] = startrecv( recvDataForI, quantity, p, stype );
            noReceives++;
        }
        else
        {
            recvDataForMe = recvDataForI;
            recvDataForMeSize = quantity;
        }
    }

    // send the data via sendPlan

    for ( PartitionId i = 0; i < sendPlan.size(); ++i )
    {
        IndexType quantity = sendPlan[i].quantity;
        IndexType offset = sendPlan[i].offset;
        const char* sendDataForI = reinterpret_cast<const char*>( sendData ) + offset * typeSize;
        PartitionId p = sendPlan[i].partitionId;
        SCAI_LOG_DEBUG( logger,
                        *this << ": send " << quantity << " elements" << " to processor " << p << " at offset " << offset )

        if ( p != getRank() )
        {
            send( sendDataForI, quantity, p, stype );
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "self-exchange of " << quantity << " elements" )
            SCAI_ASSERT_DEBUG( quantity == recvDataForMeSize, "size mismatch for self exchange" )

            getActualMemory().memcpy( recvDataForMe, sendDataForI, recvDataForMeSize * typeSize );
        }
    }

    // wait for completion of receives
    unique_ptr<MPI_Status[]> statuses( new MPI_Status[noReceives] );
    SCAI_MPICALL( logger, MPI_Waitall( noReceives, commRequest.get(), statuses.get() ), "MPI_Waitall" )
    // ToDo: check for correct sizes, was done in earlier version, but is now redundant
}

/* ---------------------------------------------------------------------------------- */

tasking::SyncToken* MPICommunicator::exchangeByPlanAsyncImpl(
    void* const recvData,
    const CommunicationPlan& recvPlan,
    const void* const sendData,
    const CommunicationPlan& sendPlan,
    const common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.exchangeByPlanAsync" )
    SCAI_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    SCAI_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )
    SCAI_LOG_INFO( logger,
                   *this << ": exchange for values of type " << stype
                   << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )
    int noRequests = sendPlan.size() + recvPlan.size();
    // create MPIToken as unique_ptr, so it will be freed in case of exception
    std::unique_ptr<MPISyncToken> pSyncToken( new MPISyncToken( noRequests ) );
    MPISyncToken& syncToken = *pSyncToken;
    void* recvDataForMe = NULL;
    IndexType recvDataForMeSize = 0;

    size_t typeSize = common::typeSize( stype );

    // setup receives for each entry in receive plan

    for ( PartitionId i = 0; i < recvPlan.size(); ++i )
    {
        IndexType quantity = recvPlan[i].quantity;
        char* recvDataForI = reinterpret_cast<char*>( recvData ) + recvPlan[i].offset * typeSize;
        PartitionId p = recvPlan[i].partitionId;
        SCAI_LOG_DEBUG( logger, *this << ": receive " << quantity << " elements" << " from processor " << p )

        if ( p != getRank() )
        {
            syncToken.pushRequest( startrecv( recvDataForI, quantity, p, stype ) );
        }
        else
        {
            recvDataForMe = recvDataForI;
            recvDataForMeSize = quantity;
        }
    }

    // send the data via sendPlan

    for ( PartitionId i = 0; i < sendPlan.size(); ++i )
    {
        IndexType quantity = sendPlan[i].quantity;
        const char* sendDataForI = reinterpret_cast<const char*>( sendData ) + sendPlan[i].offset * typeSize;
        PartitionId p = sendPlan[i].partitionId;
        SCAI_LOG_DEBUG( logger, *this << ": send " << quantity << " elements" << " to processor " << p )

        if ( p != getRank() )
        {
            syncToken.pushRequest( startsend( sendDataForI, quantity, p, stype ) );
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "self-exchange of " << quantity << " elements" )
            SCAI_ASSERT_DEBUG( quantity == recvDataForMeSize, "size mismatch for self exchange" )

            getActualMemory().memcpy( recvDataForMe, sendDataForI, recvDataForMeSize * typeSize );
        }
    }

    return pSyncToken.release();
}

/* ---------------------------------------------------------------------------------- */
/*              bcast                                                                 */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::bcastImpl( void* val, const IndexType n, const PartitionId root, common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.bcast" )
    MPI_Datatype commType = getMPIType( stype );
    SCAI_MPICALL( logger, MPI_Bcast( val, n, commType, root, mComm ), "MPI_Bcast<" << stype << ">" )
}

/* ---------------------------------------------------------------------------------- */
/*           all2allv                                                                 */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::all2allvImpl( void* recvBuffer[], const IndexType recvCount[],
                                    const void* sendBuffer[], const IndexType sendCount[],
                                    common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.all2allv" )

    int noReceives = 0;

    unique_ptr<MPI_Request[]> commRequest( new MPI_Request[getSize()] );

    for ( PartitionId i = 0; i < getSize(); ++i )
    {
        commRequest[noReceives] = startrecv( recvBuffer[i], recvCount[i], i, stype );
        noReceives++;
    }

    for ( PartitionId i = 0; i < getSize(); ++i )
    {
        send( sendBuffer[i], sendCount[i], i, stype );
    }

    // wait for completion of receives

    unique_ptr<MPI_Status[]> statuses( new MPI_Status[noReceives] );

    SCAI_MPICALL( logger, MPI_Waitall( noReceives, commRequest.get(), statuses.get() ), "MPI_Waitall" )
}

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

IndexType MPICommunicator::shiftImpl(
    void* recvVals,
    const IndexType recvSize,
    const PartitionId source,
    const void* sendVals,
    const IndexType sendSize,
    const PartitionId dest,
    common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.shift" )

    SCAI_ASSERT_NE_ERROR( source, getRank(), "source must not be this partition" )
    SCAI_ASSERT_NE_ERROR( dest, getRank(), "dest must not be this partition" )

    SCAI_LOG_DEBUG( logger,
                    *this << ": recv from " << source << " max " << recvSize << " values " << ", send to " << dest << " " << sendSize << " values." )

    MPI_Datatype commType = getMPIType( stype );

    MPI_Status mpiStatus;

    SCAI_MPICALL( logger,
                  MPI_Sendrecv( const_cast<void*>( sendVals ), sendSize, commType, dest, 4711, recvVals, recvSize,
                                commType, source, 4711, mComm, &mpiStatus ),
                  "MPI_Sendrecv" )

    // extract number of read values from status

    int count = 0;

    SCAI_MPICALL( logger, MPI_Get_count( &mpiStatus, commType, &count ),
                  "MPI_Get_count<" << stype << ">" )

    SCAI_LOG_DEBUG( logger, "received from " << source << " #values = " << count << ", max was " << recvSize )

    return count;
}

/* ---------------------------------------------------------------------------------- */
/*              shiftAsync                                                            */
/* ---------------------------------------------------------------------------------- */

tasking::SyncToken* MPICommunicator::shiftAsyncImpl(
    void* recvVals,
    const PartitionId source,
    const void* sendVals,
    const PartitionId dest,
    const IndexType size,
    common::ScalarType stype ) const
{
    SCAI_LOG_DEBUG( logger,
                    *this << ": recv from " << source << ", send to " << dest << ", both " << size << " values." )

    SCAI_ASSERT_NE_ERROR( source, getRank(), "source must not be this partition" )
    SCAI_ASSERT_NE_ERROR( dest, getRank(), "dest must not be this partition" )

    // need an MPI communicator with 2 requests, no clean up needed

    unique_ptr<MPISyncToken> pSyncToken( new MPISyncToken( 2 ) );

    pSyncToken->pushRequest( startrecv( recvVals, size, source, stype ) );
    pSyncToken->pushRequest( startsend( sendVals, size, dest, stype ) );

    return pSyncToken.release();
}

/* ---------------------------------------------------------------------------------- */
/*              scan                                                                  */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::scanImpl( void* outValues, const void* inValues, const IndexType n, common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.scan" )

    SCAI_LOG_DEBUG( logger, "MPI_Scan, n = " << n << ", type = " << stype )

    MPI_Datatype commType = getMPIType( stype );
    MPI_Op opType = getMPISum( stype );

    if ( inValues == outValues )
    {
        SCAI_MPICALL( logger, MPI_Scan( MPI_IN_PLACE, outValues, 1, commType, opType,
                                        mComm ), "MPI_Scan" )
    }
    else
    {
        SCAI_MPICALL( logger, MPI_Scan( const_cast<void*>( inValues ), outValues, n, commType, opType,
                                        mComm ), "MPI_Scan" )
    }
}

/* ---------------------------------------------------------------------------------- */
/*              sum                                                                   */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::sumImpl( void* outValues, const void* inValues, const IndexType n, common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.sum" )

    MPI_Datatype commType = getMPIType( stype );
    MPI_Op opType = getMPISum( stype );

    if ( inValues == outValues )
    {
        SCAI_MPICALL( logger, MPI_Allreduce( MPI_IN_PLACE, outValues, n, commType, opType,
                                             mComm ), "MPI_Allreduce(MPI_SUM)" )
    }
    else
    {
        SCAI_MPICALL( logger, MPI_Allreduce( const_cast<void*>( inValues ), outValues, n, commType, opType,
                                             mComm ), "MPI_Allreduce(MPI_SUM)" )
    }
}

/* ---------------------------------------------------------------------------------- */
/*              min                                                                   */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::minImpl( void* outValues, const void* inValues, const IndexType n, common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.min" )

    MPI_Datatype commType = getMPIType( stype );
    MPI_Op opType = getMPIMin( stype );

    if ( inValues == outValues )
    {
        SCAI_MPICALL( logger, MPI_Allreduce( MPI_IN_PLACE, outValues, n, commType, opType,
                                             mComm ), "MPI_Allreduce(MPI_SUM)" )
    }
    else
    {
        SCAI_MPICALL( logger, MPI_Allreduce( const_cast<void*>( inValues ), outValues, n, commType, opType,
                                             mComm ), "MPI_Allreduce(MPI_SUM)" )
    }
}

/* ---------------------------------------------------------------------------------- */
/*              min                                                                   */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::maxImpl( void* outValues, const void* inValues, const IndexType n, common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.max" )

    MPI_Datatype commType = getMPIType( stype );
    MPI_Op opType = getMPIMax( stype );

    if ( inValues == outValues )
    {
        SCAI_MPICALL( logger, MPI_Allreduce( MPI_IN_PLACE, outValues, n, commType, opType,
                                             mComm ), "MPI_Allreduce(MPI_SUM)" )
    }
    else
    {
        SCAI_MPICALL( logger, MPI_Allreduce( const_cast<void*>( inValues ), outValues, n, commType, opType,
                                             mComm ), "MPI_Allreduce(MPI_SUM)" )
    }
}

void MPICommunicator::synchronize() const
{
    SCAI_REGION( "Communicator.MPI.sync" )
    SCAI_MPICALL( logger, MPI_Barrier( mComm ), "MPI_Barrier()" )
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myVals, n, root, allVals )                                           */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::scatterImpl(
    void* myVals,
    const IndexType n,
    const PartitionId root,
    const void* allVals,
    common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.scatter" )

    SCAI_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )

    SCAI_LOG_DEBUG( logger, *this << ": scatter of " << n << " elements, root = " << root )

    MPI_Datatype commType = getMPIType( stype );

    // MPI interface is not aware of const, so const_cast is required

    SCAI_MPICALL( logger,
                  MPI_Scatter( const_cast<void*>( allVals ), n, commType, myVals, n, commType, root,
                               mComm ),
                  "MPI_Scatter<" << stype << ">" )
}

/* ---------------------------------------------------------------------------------- */
/*      scatterV( myVals, n, root, allVals, sizes )                                   */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::scatterVImpl(
    void* myVals,
    const IndexType n,
    const PartitionId root,
    const void* allVals,
    const IndexType sizes[],
    common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.scatterV" )
    SCAI_ASSERT_ERROR( root < getSize(), "illegal root, root = " << root )
    MPI_Datatype commType = getMPIType( stype );

    if ( root == getRank() )
    {
        void* sendbuf = const_cast<void*>( allVals );
        PartitionId np = getSize();
        unique_ptr<int[]> counts( new int[np] );
        unique_ptr<int[]> displs( new int[np] );
        int displacement = 0;

        for ( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = static_cast<int>( sizes[i] );
            displs[i] = displacement;
            displacement += counts[i];
        }

        SCAI_LOG_DEBUG( logger,
                        *this << ": scatter of " << displacement << " elements, I receive " << n << " elements" )
        SCAI_MPICALL( logger,
                      MPI_Scatterv( sendbuf, counts.get(), displs.get(), commType, myVals, n, commType, root, mComm ),
                      "MPI_Scatterv" )
    }
    else
    {
        // VampirTrace: requires valid counts array, even if values will be ignored
        PartitionId np = getSize();
        unique_ptr<int[]> counts( new int[np] );

        for ( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = 0;
        }

        SCAI_LOG_DEBUG( logger, *this << ": root = " << root << " scatters " << n << " elements to me" )
        SCAI_MPICALL( logger,
                      MPI_Scatterv( NULL, counts.get(), NULL, commType, myVals, n, commType, root, mComm ),
                      "MPI_Scatterv" )
    }
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allVals, n, root, myVals )                                            */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::gatherImpl(
    void* allVals,
    const IndexType n,
    const PartitionId root,
    const void* myVals,
    const common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.gather" )

    SCAI_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )

    SCAI_LOG_DEBUG( logger, *this << ": gather of " << n << " elements, root = " << root )

    MPI_Datatype commType = getMPIType( stype );

    void* sendbuf = const_cast<void*>( myVals );  // MPI interface is not const aware

    SCAI_MPICALL( logger,
                  MPI_Gather( sendbuf, n, commType, allVals, n, commType, root, mComm ),
                  "MPI_Gather<" << stype << ">" )
}

/* ---------------------------------------------------------------------------------- */
/*      gatherV( allVals, n, root, myVals, sizes )                                    */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::gatherVImpl(
    void* allVals,
    const IndexType n,
    const PartitionId root,
    const void* myVals,
    const IndexType sizes[],
    const common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.gatherV" )
    SCAI_ASSERT_ERROR( root < getSize(), "illegal root, root = " << root )
    void* sendbuf = const_cast<void*>( myVals );
    MPI_Datatype commType = getMPIType( stype );

    if ( root == getRank() )
    {
        PartitionId np = getSize();
        unique_ptr<int[]> counts( new int[np] );
        unique_ptr<int[]> displs( new int[np] );
        int displacement = 0;

        for ( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = static_cast<int>( sizes[i] );
            displs[i] = displacement;
            displacement += counts[i];
        }

        SCAI_LOG_DEBUG( logger,
                        *this << ": scatter of " << displacement << " elements, I receive " << n << " elements" )
        SCAI_MPICALL( logger,
                      MPI_Gatherv( sendbuf, n, commType, allVals, counts.get(), displs.get(), commType, root, mComm ),
                      "MPI_Gatherv<" << stype << ">" )
    }
    else
    {
        // VampirTrace: requires valid counts array, even if values will be ignored

        PartitionId np = getSize();

        unique_ptr<int[]> counts( new int[np] );

        for ( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = 0;
        }

        SCAI_LOG_DEBUG( logger, *this << ": root = " << root << " scatters " << n << " elements to me" )

        SCAI_MPICALL( logger,
                      MPI_Gatherv( sendbuf, n, commType, NULL, counts.get(), NULL, commType, root, mComm ),
                      "MPI_Gatherv<" << stype << ">" )
    }
}

/* ---------------------------------------------------------------------------------- */
/*           maxloc                                                                   */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::maxlocImpl( void* val, IndexType* location, PartitionId root, common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.maxloc" )

    MPI_Datatype commType = getMPI2Type( stype, common::TypeTraits<IndexType>::stype );

    // verify that the value for location follows directly the location of val

    size_t typeSize = common::typeSize( stype );
    size_t tmpSize  = typeSize + sizeof( IndexType );

    void* expectedLoc = reinterpret_cast<char*>( val ) + typeSize;

    SCAI_ASSERT_EQ_ERROR( expectedLoc, location,
                          "val and loc not contiguously in memory"
                          << ": val type = " << stype << ", size = " << typeSize
                          << ", index type = " << common::TypeTraits<IndexType>::id() << ", size = " << sizeof( IndexType ) )

    // only on Host

    std::unique_ptr<char[]> tmp ( new char[ tmpSize ] );
    memcpy( tmp.get(), val,  tmpSize );

    MPI_Reduce( tmp.get(), val, 1, commType, MPI_MAXLOC, root, mComm );
}

/* ---------------------------------------------------------------------------------- */
/*           minloc                                                                   */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::minlocImpl( void* val, IndexType* location, PartitionId root, common::ScalarType stype ) const
{
    SCAI_REGION( "Communicator.MPI.minloc" )

    MPI_Datatype commType = getMPI2Type( stype, common::TypeTraits<IndexType>::stype );

    // verify that the value for location follows directly the location of val

    size_t typeSize = common::typeSize( stype );
    size_t tmpSize  = typeSize + sizeof( IndexType );

    void* expectedLoc = reinterpret_cast<char*>( val ) + typeSize;

    SCAI_ASSERT_EQ_ERROR( expectedLoc, location,
                          "val and loc not contiguously in memory"
                          << ": val type = " << stype << ", size = " << typeSize
                          << ", index type = " << common::TypeTraits<IndexType>::id() << ", size = " << sizeof( IndexType ) )

    // only on host

    std::unique_ptr<char[]> tmp ( new char[ tmpSize ] );
    memcpy( tmp.get(), val,  tmpSize );

    MPI_Reduce( tmp.get(), val, 1, commType, MPI_MINLOC, root, mComm );
}

/* ---------------------------------------------------------------------------------- */
/*          supportsLocReduction                                                      */
/* ---------------------------------------------------------------------------------- */

bool MPICommunicator::supportsLocReduction( common::ScalarType vType, common::ScalarType iType ) const
{
    if ( common::TypeTraits<int>::stype != iType )
    {
        // IndexType != int is not supported here
        return false;
    }

    // min/maxloc reduction not supported for all data types

    switch ( vType )
    {
        case common::ScalarType::INT       :
            return true;
        case common::ScalarType::FLOAT     :
            return true;
        case common::ScalarType::DOUBLE    :
            return true;
        default                        :
            return false;
    }
}

/* ---------------------------------------------------------------------------------- */
/*           swap                                                                     */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::swapImpl( void* val, const IndexType n, PartitionId partner, common::ScalarType stype ) const
{
    if ( partner == getRank() )
    {
        return;   // swap with same processor is redundant
    }

    SCAI_REGION( "Communicator.MPI.swap" )

    // only on host

    std::unique_ptr<char[]> tmp ( new char[ n * common::typeSize( stype ) ] );

    memcpy( tmp.get(), val, n * common::typeSize( stype ) );

    MPI_Status mpiStatus;

    MPI_Datatype commType = getMPIType( stype );

    SCAI_MPICALL( logger,
                  MPI_Sendrecv( tmp.get(), n, commType, partner, defaultTag, val, n, commType, partner, defaultTag,
                                mComm, &mpiStatus ),
                  "MPI_Sendrecv<" << stype << "> for swap" )

    IndexType nPartner = getCount( mpiStatus, stype );

    SCAI_ASSERT_EQ_ERROR( n, nPartner, "size mismatch for swap" )
}

hmemo::ContextPtr MPICommunicator::getCommunicationContext( const hmemo::_HArray& array ) const
{
    // get a valid context, i.e. a context that contains valid data
    hmemo::ContextPtr validContext = array.getValidContext();
    SCAI_LOG_DEBUG( logger, "CommunicationContext: valid context for " << array << ": " << *validContext )

    if ( validContext->getType() == common::ContextType::Host )
    {
        return validContext;
    }

    // This can only be used for CUDAaware MPI

    if ( isCUDAAware && ( validContext->getType() == common::ContextType::CUDA ) )
    {
        return validContext;
    }

    return hmemo::Context::getHostPtr();
}

/* ---------------------------------------------------------------------------------- */
/*           writeAt                                                                  */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::writeAt( std::ostream& stream ) const
{
    // Info about rank and size of the communicator is very useful

    if ( mKind == MPICommKind::CREATED )
    {
        stream << "MPI_created(" << getRank() << ":" << getSize() << ")";
    }
    else if ( mKind == MPICommKind::INTERNAL )
    {
        stream << "MPI(" << getRank() << ":" << getSize() << ")";
    }
    else if ( mKind == MPICommKind::EXTERNAL )
    {
        stream << "_MPI(" << getRank() << ":" << getSize() << ")";
    }
}

/* --------------------------------------------------------------- */

MPICommunicator::MPICommunicator( const MPICommunicator& comm, int color, int key ) :

    Communicator( MPI ),
    mKind( MPICommKind::CREATED ),
    mThreadSafetyLevel( comm.mThreadSafetyLevel )

{
    SCAI_LOG_INFO( logger, *this << ": split ( color = " << color << ", key = " << key << " )" )
    SCAI_MPICALL( logger, MPI_Comm_split( comm.mComm, color, key, &mComm ), "comm split" );
    
    {
        // temporary variables required as PartitionId can be different from int

        int mpiRank;
        int mpiSize;

        SCAI_MPICALL( logger, MPI_Comm_size( mComm, &mpiSize ), "MPI_Comm_size" )
        SCAI_MPICALL( logger, MPI_Comm_rank( mComm, &mpiRank ), "MPI_Comm_rank" )

        setSizeAndRank( static_cast<PartitionId>( mpiSize ), static_cast<PartitionId>( mpiRank ) );
    }
}

/* --------------------------------------------------------------- */

MPICommunicator* MPICommunicator::splitIt( PartitionId color, PartitionId key ) const
{
    return new MPICommunicator( *this, color, key );
}

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

MPICommunicator::MPIGuard MPICommunicator::guard;

CommunicatorPtr MPICommunicator::create()
{
    if ( !theMPICommunicator )
    {
        SCAI_LOG_INFO( logger, "create new MPICommunicator" )
        // create a new instance of MPICommunicator, will call MPI_Init
        theMPICommunicator = shared_ptr<const MPICommunicator>( new MPICommunicator() );
    }

    return theMPICommunicator;
}

/* --------------------------------------------------------------- */

Communicator::CommunicatorKind MPICommunicator::createValue()
{
    return MPI;
}

/* --------------------------------------------------------------- */

} /* end namespace dmemo */

} /* end namespace scai */
