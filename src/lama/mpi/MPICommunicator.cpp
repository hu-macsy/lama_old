/**
 * @file MPICommunicator.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

// hpp
#include <lama/mpi/MPICommunicator.hpp>

// others
#include <lama/mpi/MPISyncToken.hpp>
#include <lama/mpi/MPIUtils.hpp>

#include <lama/exception/LAMAAssert.hpp>

// tracing
#include <lama/tracing.hpp>

// boost
#include <boost/scoped_array.hpp>
#include <boost/bind.hpp>

#include <iostream>
#include <algorithm>

#include <mpi.h>

using namespace std;

namespace lama
{

const int MPICommunicator::defaultTag = 1;

LAMA_LOG_DEF_LOGGER( MPICommunicator::logger, "Communicator.MPICommunicator" )

MPICommunicator::MPICommunicator( int& argc, char** & argv )
    : Communicator( "MPI" ), mThreadSafetyLevel( Communicator::Funneled )
{
    int initialized = 0;

    LAMA_MPICALL( logger, MPI_Initialized( &initialized ), "MPI_Initialized" );

    if ( initialized )
    {
        mExternInitialization = true;
        LAMA_LOG_WARN( logger, "MPI_Init: MPI has already been initialized." )
    }
    else
    {
        mExternInitialization = false;

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
            LAMA_LOG_ERROR( logger,
                            "LAMA_MPICOMM_THREAD_SAFETY = " << threadSafetyEnvironment << ", unknown value (try Multiple or Serialized)" )
        }

        int providedThreadSafety = MPI_THREAD_SINGLE;

        LAMA_MPICALL( logger, MPI_Init_thread( &argc, &argv, requiredThreadSafety, &providedThreadSafety ),
                      "MPI_Init" );

        switch ( providedThreadSafety )
        {
        case MPI_THREAD_MULTIPLE:
            mThreadSafetyLevel = Communicator::Multiple;
            LAMA_LOG_INFO( logger, "MPI Thread Safety Level: MPI_THREAD_MULTIPLE" )
            break;
        case MPI_THREAD_SERIALIZED:
            mThreadSafetyLevel = Communicator::Serialized;
            LAMA_LOG_INFO( logger, "MPI Thread Safety Level: MPI_THREAD_SERIALIZED" )
            break;
        case MPI_THREAD_FUNNELED:
            mThreadSafetyLevel = Communicator::Funneled;
            LAMA_LOG_INFO( logger, "MPI Thread Safety Level: MPI_THREAD_FUNNELED" )
            break;
        case MPI_THREAD_SINGLE:
            mThreadSafetyLevel = Communicator::Funneled;
            LAMA_LOG_INFO( logger, "MPI Thread Safety Level: MPI_THREAD_SINGLE" )
            break;
        default:
            MPI_Finalize();
            LAMA_THROWEXCEPTION( "MPI which supports at leas thread level MPI_THREAD_FUNNELED is required." )
        }

    }

    mCommWorld = MPI_COMM_WORLD;
    MPI_Comm_dup( mCommWorld, &mComm );
    MPI_Comm_dup( mCommWorld, &mCommTask );
    LAMA_LOG_INFO( logger, "MPI_Init" )
    LAMA_MPICALL( logger, MPI_Comm_size( mComm, &mSize ), "MPI_Comm_size" );
    LAMA_MPICALL( logger, MPI_Comm_rank( mComm, &mRank ), "MPI_Comm_rank" );
}

MPICommunicator::~MPICommunicator()
{
    LAMA_LOG_INFO( logger, *this << ": ~MPICommunicator" )
    int finalized = 0;
    LAMA_MPICALL( logger, MPI_Finalized( &finalized ), "MPI_Finalized" );

    if ( !finalized )
    {
        if ( !mExternInitialization )
        {
            LAMA_LOG_INFO( logger, "call MPI_Finalize" )
            LAMA_MPICALL( logger, MPI_Finalize(), "MPI_Finalize" );
        }
        else
        {
            LAMA_LOG_INFO( logger, "ATTENTION: no call MPI_Finalize, was externally initialized" )
        }
    }
    else
    {
        LAMA_LOG_WARN( logger, *this << ": tried to finalize MPI, but MPI has already been finalized." )
    }
}

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

MPI_Comm MPICommunicator::getMPIComm() const
{
    return mComm;
}

template<typename T>
MPI_Request MPICommunicator::startrecv( T* buffer, int count, int source ) const
{
    MPI_Request request;
    MPI_Datatype commType = getMPIType<T>();
    LAMA_MPICALL( logger, MPI_Irecv( buffer, count, commType, source, defaultTag, selectMPIComm(), &request ),
                  "MPI_Irecv" );
    return request;
}

template<typename T>
MPI_Request MPICommunicator::startsend( const T* buffer, int count, int target ) const
{
    MPI_Request request;
    MPI_Datatype commType = getMPIType<T>();
    LAMA_MPICALL( logger,
                  MPI_Isend( const_cast<T*>( buffer ), count, commType, target, defaultTag, selectMPIComm(), &request ),
                  "MPI_Isend" );
    return request;
}

template<typename T>
int MPICommunicator::getCount( MPI_Status& mpiStatus ) const
{
    int size = 0;
    MPI_Datatype commType = getMPIType<T>();
    LAMA_MPICALL( logger, MPI_Get_count( &mpiStatus, commType, &size ), "MPI_Get_count" );
    return size;
}

template<typename T>
void MPICommunicator::send( const T* buffer, int count, int target ) const
{
    MPI_Datatype commType = getMPIType<T>();
    LAMA_MPICALL( logger, MPI_Send( const_cast<T*>( buffer ), count, commType, target, defaultTag, selectMPIComm() ),
                  "MPI_Send" );
}

/* ---------------------------------------------------------------------------------- */

void MPICommunicator::all2all( int* recvSizes, const int* sendSizes ) const
{
    LAMA_ASSERT_ERROR( sendSizes != 0, " invalid sendSizes " )
    LAMA_ASSERT_ERROR( recvSizes != 0, " invalid recvSizes " )
    LAMA_MPICALL( logger,
                  MPI_Alltoall( const_cast<int*>( sendSizes ), 1, MPI_INT, recvSizes, 1, MPI_INT, selectMPIComm() ),
                  "MPI_Alltoall" );
}

/* ---------------------------------------------------------------------------------- */

void MPICommunicator::exchangeByPlan(
    int* const recvData,
    const CommunicationPlan& recvPlan,
    const int* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan );
}

void MPICommunicator::exchangeByPlan(
    float* const recvData,
    const CommunicationPlan& recvPlan,
    const float* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan );
}

void MPICommunicator::exchangeByPlan(
    double* const recvData,
    const CommunicationPlan& recvPlan,
    const double* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    exchangeByPlanImpl( recvData, recvPlan, sendData, sendPlan );
}

std::auto_ptr<SyncToken> MPICommunicator::exchangeByPlanAsync(
    int* const recvData,
    const CommunicationPlan& recvPlan,
    const int* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    return exchangeByPlanAsyncImpl( recvData, recvPlan, sendData, sendPlan );
}

std::auto_ptr<SyncToken> MPICommunicator::exchangeByPlanAsync(
    float* const recvData,
    const CommunicationPlan& recvPlan,
    const float* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    return exchangeByPlanAsyncImpl( recvData, recvPlan, sendData, sendPlan );
}

std::auto_ptr<SyncToken> MPICommunicator::exchangeByPlanAsync(
    double* const recvData,
    const CommunicationPlan& recvPlan,
    const double* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    return exchangeByPlanAsyncImpl( recvData, recvPlan, sendData, sendPlan );
}

template<typename T>
void MPICommunicator::exchangeByPlanImpl(
    T* const recvData,
    const CommunicationPlan& recvPlan,
    const T* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    LAMA_REGION( "Communicator.MPI.exchangeByPlan" )
    LAMA_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    LAMA_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )
    LAMA_LOG_INFO( logger,
                   *this << ": exchange for values of type " << typeid( T ).name() << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )
    int maxReceives = recvPlan.size();
    int noReceives = 0; // will be incremented
    T* recvDataForMe = NULL;
    IndexType recvDataForMeSize = 0;
    boost::scoped_array<MPI_Request> commRequest( new MPI_Request[maxReceives] );

    // setup receives for each entry in receive plan

    for ( PartitionId i = 0; i < maxReceives; ++i )
    {
        IndexType quantity = recvPlan[i].quantity;
        IndexType offset = recvPlan[i].offset;
        T* recvDataForI = recvData + offset;
        PartitionId p = recvPlan[i].partitionId;
        LAMA_LOG_DEBUG( logger,
                        *this << ": receive " << quantity << " elements" << " from processor " << p << " at offset " << offset )

        if ( p != mRank )
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

    for ( PartitionId i = 0; i < sendPlan.size(); ++i )
    {
        IndexType quantity = sendPlan[i].quantity;
        IndexType offset = sendPlan[i].offset;
        const T* sendDataForI = sendData + offset;
        PartitionId p = sendPlan[i].partitionId;
        LAMA_LOG_DEBUG( logger,
                        *this << ": send " << quantity << " elements" << " to processor " << p << " at offset " << offset )

        if ( p != mRank )
        {
            send( sendDataForI, quantity, p );
        }
        else
        {
            LAMA_LOG_DEBUG( logger, "self-exchange of " << quantity << " elements" )
            LAMA_ASSERT_DEBUG( quantity == recvDataForMeSize, "size mismatch for self exchange" )

            for ( IndexType k = 0; k < recvDataForMeSize; k++ )
            {
                recvDataForMe[k] = sendDataForI[k];
            }
        }
    }

    // wait for completion of receives
    boost::scoped_array<MPI_Status> statuses( new MPI_Status[noReceives] );
    LAMA_MPICALL( logger, MPI_Waitall( noReceives, commRequest.get(), statuses.get() ), "MPI_Waitall" );
    // ToDo: check for correct sizes, was done in earlier version, but is now redundant
}

/* ---------------------------------------------------------------------------------- */

template<typename T>
auto_ptr<SyncToken> MPICommunicator::exchangeByPlanAsyncImpl(
    T* const recvData,
    const CommunicationPlan& recvPlan,
    const T* const sendData,
    const CommunicationPlan& sendPlan ) const
{
    LAMA_ASSERT_ERROR( sendPlan.allocated(), "sendPlan not allocated" )
    LAMA_ASSERT_ERROR( recvPlan.allocated(), "recvPlan not allocated" )
    LAMA_LOG_INFO( logger,
                   *this << ": exchange for values of type " << typeid( T ).name() << ", send to " << sendPlan.size() << " processors, recv from " << recvPlan.size() )
    int noRequests = sendPlan.size() + recvPlan.size();

    // create MPIToken as auto_ptr, so it will be freed in case of exception

    auto_ptr<MPISyncToken> pSyncToken( new MPISyncToken( noRequests ) );

    MPISyncToken& syncToken = *pSyncToken;

    T* recvDataForMe = NULL;
    IndexType recvDataForMeSize = 0;

    // setup receives for each entry in receive plan

    for ( PartitionId i = 0; i < recvPlan.size(); ++i )
    {
        IndexType quantity = recvPlan[i].quantity;
        T* recvDataForI = recvData + recvPlan[i].offset;
        PartitionId p = recvPlan[i].partitionId;
        LAMA_LOG_DEBUG( logger, *this << ": receive " << quantity << " elements" << " from processor " << p )

        if ( p != mRank )
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

    for ( PartitionId i = 0; i < sendPlan.size(); ++i )
    {
        IndexType quantity = sendPlan[i].quantity;
        const T* sendDataForI = sendData + sendPlan[i].offset;
        PartitionId p = sendPlan[i].partitionId;
        LAMA_LOG_DEBUG( logger, *this << ": send " << quantity << " elements" << " to processor " << p )

        if ( p != mRank )
        {
            syncToken.pushRequest( startsend( sendDataForI, quantity, p ) );
        }
        else
        {
            LAMA_LOG_DEBUG( logger, "self-exchange of " << quantity << " elements" )
            LAMA_ASSERT_DEBUG( quantity == recvDataForMeSize, "size mismatch for self exchange" )

            for ( IndexType k = 0; k < recvDataForMeSize; k++ )
            {
                recvDataForMe[k] = sendDataForI[k];
            }
        }
    }

    return auto_ptr<SyncToken>( pSyncToken.release() );
}

/* ---------------------------------------------------------------------------------- */

inline MPI_Comm MPICommunicator::selectMPIComm() const
{
    if ( mThreadSafetyLevel != Multiple )
    {
        return mComm;
    }

    const boost::thread thisThread;

    if ( thisThread == mMainThread )
    {
        return mComm;
    }
    else
    {
        return mCommTask;
    }
}

/* ---------------------------------------------------------------------------------- */
/*              getMPIType                                                            */
/* ---------------------------------------------------------------------------------- */

template<>
inline MPI_Datatype MPICommunicator::getMPIType<float>()
{
    return MPI_FLOAT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<double>()
{
    return MPI_DOUBLE;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<int>()
{
    return MPI_INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPIType<unsigned long>()
{
    return MPI_UNSIGNED_LONG;
}

/* ---------------------------------------------------------------------------------- */
/*              getMPI2Type                                                           */
/* ---------------------------------------------------------------------------------- */

template<typename T1,typename T2>
inline MPI_Datatype MPICommunicator::getMPI2Type()
{
    LAMA_THROWEXCEPTION( "unsupported type for MPI communication" )
    return MPI_2INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<float,int>()
{
    return MPI_FLOAT_INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<double,int>()
{
    return MPI_DOUBLE_INT;
}

template<>
inline MPI_Datatype MPICommunicator::getMPI2Type<int,int>()
{
    return MPI_2INT;
}

/* ---------------------------------------------------------------------------------- */
/*              shift                                                                 */
/* ---------------------------------------------------------------------------------- */

template<typename T>
IndexType MPICommunicator::shiftMPI(
    T recvVals[],
    const IndexType recvSize,
    const PartitionId source,
    const T sendVals[],
    const IndexType sendSize,
    const PartitionId dest ) const
{
    LAMA_ASSERT_ERROR( source != getRank(), "source must not be this partition" )
    LAMA_ASSERT_ERROR( dest != getRank(), "dest must not be this partition" )

    LAMA_LOG_DEBUG( logger,
                    *this << ": recv from " << source << " max " << recvSize << " values " << ", send to " << dest << " " << sendSize << " values." )

    MPI_Datatype commType = getMPIType<T>();
    MPI_Status mpiStatus;

    LAMA_MPICALL( logger,
                  MPI_Sendrecv( const_cast<T*>( sendVals ), sendSize, commType, dest, 4711, recvVals, recvSize, commType, source, 4711, selectMPIComm(), &mpiStatus ),
                  "MPI_Sendrecv" );

    // extract number of read values from status

    int count = 0;

    LAMA_MPICALL( logger, MPI_Get_count( &mpiStatus, commType, &count ), "MPI_Get_count(T)" );

    LAMA_LOG_DEBUG( logger, "received from " << source << " #values = " << count << ", max was " << recvSize )

    return count;
}

IndexType MPICommunicator::shiftImpl(
    double recvData[],
    const IndexType recvSize,
    const double sendVals[],
    const IndexType sendSize,
    const int direction ) const
{
    LAMA_LOG_DEBUG( logger,
                    *this << ": shift, direction = " << direction << ", sendsize = " << sendSize << ", recvsize = " << recvSize )

    if ( direction % getSize() == 0 )
    {
        return shift0( recvData, recvSize, sendVals, sendSize );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftMPI( recvData, recvSize, source, sendVals, sendSize, dest );
}

IndexType MPICommunicator::shiftImpl(
    float recvData[],
    const IndexType recvSize,
    const float sendVals[],
    const IndexType sendSize,
    const int direction ) const
{
    LAMA_LOG_DEBUG( logger,
                    *this << ": shift, direction = " << direction << ", sendsize = " << sendSize << ", recvsize = " << recvSize )

    if ( direction % getSize() == 0 )
    {
        return shift0( recvData, recvSize, sendVals, sendSize );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftMPI( recvData, recvSize, source, sendVals, sendSize, dest );
}

IndexType MPICommunicator::shiftImpl(
    int recvData[],
    const IndexType recvSize,
    const int sendVals[],
    const IndexType sendSize,
    const int direction ) const
{
    LAMA_LOG_DEBUG( logger,
                    *this << ": shift, direction = " << direction << ", sendsize = " << sendSize << ", recvsize = " << recvSize )

    if ( direction % getSize() == 0 )
    {
        return shift0( recvData, recvSize, sendVals, sendSize );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftMPI( recvData, recvSize, source, sendVals, sendSize, dest );
}

/* ---------------------------------------------------------------------------------- */
/*              shiftAsync                                                            */
/* ---------------------------------------------------------------------------------- */

template<typename T>
auto_ptr<SyncToken> MPICommunicator::shiftAsyncMPI(
    T recvVals[],
    const PartitionId source,
    const T sendVals[],
    const PartitionId dest,
    const IndexType size ) const
{
    LAMA_LOG_DEBUG( logger,
                    *this << ": recv from " << source << ", send to " << dest << ", both " << size << " values." )

    LAMA_ASSERT_ERROR( source != getRank(), "source must not be this partition" )
    LAMA_ASSERT_ERROR( dest != getRank(), "dest must not be this partition" )

    // need an MPI communicator with 2 requests, no clean up needed

    auto_ptr<MPISyncToken> pSyncToken( new MPISyncToken( 2 ) );

    pSyncToken->pushRequest( startrecv( recvVals, size, source ) );
    pSyncToken->pushRequest( startsend( sendVals, size, dest ) );

    return auto_ptr<SyncToken>( pSyncToken.release() );
}

auto_ptr<SyncToken> MPICommunicator::shiftAsyncImpl(
    double recvVals[],
    const double sendVals[],
    const IndexType size,
    const int direction ) const
{
    LAMA_LOG_DEBUG( logger, *this << ": shiftAsync size = " << size << ", direction = " << direction )

    if ( direction % getSize() == 0 )
    {
        return defaultShiftAsync( recvVals, sendVals, size, 0 );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftAsyncMPI( recvVals, source, sendVals, dest, size );
}

auto_ptr<SyncToken> MPICommunicator::shiftAsyncImpl(
    float recvVals[],
    const float sendVals[],
    const IndexType size,
    const int direction ) const
{
    if ( direction % getSize() == 0 )
    {
        return defaultShiftAsync( recvVals, sendVals, size, 0 );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftAsyncMPI( recvVals, source, sendVals, dest, size );
}

auto_ptr<SyncToken> MPICommunicator::shiftAsyncImpl(
    int recvVals[],
    const int sendVals[],
    const IndexType size,
    const int direction ) const
{
    if ( direction % getSize() == 0 )
    {
        return defaultShiftAsync( recvVals, sendVals, size, 0 );
    }

    PartitionId dest = getNeighbor( direction );
    PartitionId source = getNeighbor( -direction );
    return shiftAsyncMPI( recvVals, source, sendVals, dest, size );
}

/* ---------------------------------------------------------------------------------- */
/*              sum                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename T>
T MPICommunicator::sumImpl( const T value ) const
{
    T sum;
    MPI_Datatype commType = getMPIType<T>();
    LAMA_MPICALL( logger, MPI_Allreduce( ( void* ) &value, ( void* ) &sum, 1, commType, MPI_SUM, selectMPIComm() ),
                  "MPI_Allreduce(MPI_SUM)" );
    return sum;
}

float MPICommunicator::sum( const float value ) const
{
    return sumImpl( value );
}

double MPICommunicator::sum( const double value ) const
{
    return sumImpl( value );
}

size_t MPICommunicator::sum( const size_t value ) const
{
    // conversion: size_t <-> unsigned long should always be okay

    return sumImpl<unsigned long>( static_cast<unsigned long>( value ) );
}

int MPICommunicator::sum( const int value ) const
{
    return sumImpl( value );
}

/* ---------------------------------------------------------------------------------- */
/*              min / max reduction                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename T>
T MPICommunicator::minval( const T value ) const
{
    MPI_Datatype commType = getMPIType<T>();

    T globalMin; // no initialization needed, done in MPI call

    LAMA_MPICALL( logger,
                  MPI_Allreduce( ( void* ) &value, ( void* ) &globalMin, 1, commType, MPI_MIN, selectMPIComm() ),
                  "MPI_Allreduce( MPI_MIN )" );
    return globalMin;
}

template<typename T>
T MPICommunicator::maxval( const T value ) const
{
    MPI_Datatype commType = getMPIType<T>();

    T globalMax; // no initialization needed, done in MPI call

    LAMA_MPICALL( logger,
                  MPI_Allreduce( ( void* ) &value, ( void* ) &globalMax, 1, commType, MPI_MAX, selectMPIComm() ),
                  "MPI_Allreduce( MPI_MAX )" );
    return globalMax;
}

float MPICommunicator::min( const float value ) const
{
    return minval( value );
}

float MPICommunicator::max( const float value ) const
{
    return maxval( value );
}

double MPICommunicator::min( const double value ) const
{
    return minval( value );
}

double MPICommunicator::max( const double value ) const
{
    return maxval( value );
}

int MPICommunicator::min( const int value ) const
{
    return minval( value );
}

int MPICommunicator::max( const int value ) const
{
    return maxval( value );
}

void MPICommunicator::gather( vector<float>& values, float value ) const
{
    // build a vector of just a single value
    values.clear();
    values.resize( mSize, 0.0 );
    LAMA_MPICALL( logger, MPI_Allgather( &value, 1, MPI_FLOAT, &values[0], 1, MPI_FLOAT, selectMPIComm() ),
                  "MPI_Allgather(MPI_FLOAT,MPI_FLOAT)" );
}

void MPICommunicator::synchronize() const
{
    LAMA_MPICALL( logger, MPI_Barrier( selectMPIComm() ), "MPI_Barrier()" );
}

/* ---------------------------------------------------------------------------------- */
/*      bcast                                                                         */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::bcast( int val[], const IndexType n, const PartitionId root ) const
{
    LAMA_MPICALL( logger, MPI_Bcast( val, n, MPI_INT, root, selectMPIComm() ), "MPI_Bcast<int>" );
}

void MPICommunicator::bcast( double val[], const IndexType n, const PartitionId root ) const
{
    LAMA_MPICALL( logger, MPI_Bcast( val, n, MPI_DOUBLE, root, selectMPIComm() ), "MPI_Bcast<double>" );
}

void MPICommunicator::bcast( float val[], const IndexType n, const PartitionId root ) const
{
    LAMA_MPICALL( logger, MPI_Bcast( val, n, MPI_FLOAT, root, selectMPIComm() ), "MPI_Bcast<float>" );
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myvals, n, root, allvals )                                           */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void MPICommunicator::scatterImpl( T myvals[], const IndexType n, const PartitionId root, const T allvals[] ) const
{
    LAMA_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
    LAMA_LOG_DEBUG( logger, *this << ": scatter of " << n << " elements, root = " << root )
    MPI_Datatype commType = getMPIType<T>();
    // MPI interface is not aware of const, so const_cast is required
    LAMA_MPICALL( logger,
                  MPI_Scatter( const_cast<T*>( allvals ), n, commType, myvals, n, commType, root, selectMPIComm() ),
                  "MPI_Scatter" );
}

void MPICommunicator::scatter( float myvals[], const IndexType n, const PartitionId root, const float allvals[] ) const
{
    scatterImpl( myvals, n, root, allvals );
}

void MPICommunicator::scatter(
    double myvals[],
    const IndexType n,
    const PartitionId root,
    const double allvals[] ) const
{
    scatterImpl( myvals, n, root, allvals );
}

void MPICommunicator::scatter( int myvals[], const IndexType n, const PartitionId root, const int allvals[] ) const
{
    scatterImpl( myvals, n, root, allvals );
}

/* ---------------------------------------------------------------------------------- */
/*      scatter( myvals, n, root, allvals, sizes )                                    */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void MPICommunicator::scatterImpl(
    T myvals[],
    const IndexType n,
    const PartitionId root,
    const T allvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_ERROR( root < getSize(), "illegal root, root = " << root )
    MPI_Datatype commType = getMPIType<T>();

    if ( root == getRank() )
    {
        void* sendbuf = const_cast<T*>( allvals );
        PartitionId np = getSize();
        boost::scoped_array<int> counts( new int[np] );
        boost::scoped_array<int> displs( new int[np] );
        int displacement = 0;

        for ( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = static_cast<int>( sizes[i] );
            displs[i] = displacement;
            displacement += counts[i];
        }

        LAMA_LOG_DEBUG( logger,
                        *this << ": scatter of " << displacement << " elements, I receive " << n << " elements" )
        LAMA_MPICALL( logger,
                      MPI_Scatterv( sendbuf, counts.get(), displs.get(), commType, myvals, n, commType, root, selectMPIComm() ),
                      "MPI_Scatterv" );
    }
    else
    {
        // VampirTrace: requires valid counts array, even if values will be ignored

        PartitionId np = getSize();
        boost::scoped_array<int> counts( new int[np] );
        for ( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = 0;
        }

        LAMA_LOG_DEBUG( logger, *this << ": root = " << root << " scatters " << n << " elements to me" )
        LAMA_MPICALL( logger,
                      MPI_Scatterv( NULL, counts.get(), NULL, commType, myvals, n, commType, root, selectMPIComm() ),
                      "MPI_Scatterv" );
    }
}

void MPICommunicator::scatter(
    float myvals[],
    const IndexType n,
    const PartitionId root,
    const float allvals[],
    const IndexType sizes[] ) const
{
    scatterImpl( myvals, n, root, allvals, sizes );
}

void MPICommunicator::scatter(
    double myvals[],
    const IndexType n,
    const PartitionId root,
    const double allvals[],
    const IndexType sizes[] ) const
{
    scatterImpl( myvals, n, root, allvals, sizes );
}

void MPICommunicator::scatter(
    int myvals[],
    const IndexType n,
    const PartitionId root,
    const int allvals[],
    const IndexType sizes[] ) const
{
    scatterImpl( myvals, n, root, allvals, sizes );
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allvals, n, root, myvals )                                            */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void MPICommunicator::gatherImpl( T allvals[], const IndexType n, const PartitionId root, const T myvals[] ) const
{
    LAMA_ASSERT_DEBUG( root < getSize(), "illegal root, root = " << root )
    LAMA_LOG_DEBUG( logger, *this << ": gather of " << n << " elements, root = " << root )
    MPI_Datatype commType = getMPIType<T>();
    // MPI interface is not aware of const, so const_cast is required
    void* sendbuf = const_cast<T*>( myvals );
    LAMA_MPICALL( logger, MPI_Gather ( sendbuf, n, commType, allvals, n, commType, root, selectMPIComm() ),
                  "MPI_Gather<T>" );
}

void MPICommunicator::gather( double allvals[], const IndexType n, const PartitionId root, const double myvals[] ) const
{
    gatherImpl( allvals, n, root, myvals );
}

void MPICommunicator::gather( float allvals[], const IndexType n, const PartitionId root, const float myvals[] ) const
{
    gatherImpl( allvals, n, root, myvals );
}

void MPICommunicator::gather( int allvals[], const IndexType n, const PartitionId root, const int myvals[] ) const
{
    gatherImpl( allvals, n, root, myvals );
}

/* ---------------------------------------------------------------------------------- */
/*      gather( allvals, n, root, myvals, sizes )                                     */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void MPICommunicator::gatherImpl(
    T allvals[],
    const IndexType n,
    const PartitionId root,
    const T myvals[],
    const IndexType sizes[] ) const
{
    LAMA_ASSERT_ERROR( root < getSize(), "illegal root, root = " << root )
    void* sendbuf = const_cast<T*>( myvals );
    MPI_Datatype commType = getMPIType<T>();

    if ( root == getRank() )
    {
        PartitionId np = getSize();
        boost::scoped_array<int> counts( new int[np] );
        boost::scoped_array<int> displs( new int[np] );
        int displacement = 0;

        for ( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = static_cast<int>( sizes[i] );
            displs[i] = displacement;
            displacement += counts[i];
        }

        LAMA_LOG_DEBUG( logger,
                        *this << ": scatter of " << displacement << " elements, I receive " << n << " elements" )
        LAMA_MPICALL( logger,
                      MPI_Gatherv( sendbuf, n, commType, allvals, counts.get(), displs.get(), commType, root, selectMPIComm() ),
                      "MPI_Gatherv<T>" );
    }
    else
    {
        // VampirTrace: requires valid counts array, even if values will be ignored

        PartitionId np = getSize();
        boost::scoped_array<int> counts( new int[np] );
        for ( PartitionId i = 0; i < np; i++ )
        {
            counts[i] = 0;
        }

        LAMA_LOG_DEBUG( logger, *this << ": root = " << root << " scatters " << n << " elements to me" )
        LAMA_MPICALL( logger,
                      MPI_Gatherv( sendbuf, n, commType, NULL, counts.get(), NULL, commType, root, selectMPIComm() ),
                      "MPI_Gatherv<T>" );
    }
}

void MPICommunicator::gather(
    float allvals[],
    const IndexType n,
    const PartitionId root,
    const float myvals[],
    const IndexType sizes[] ) const
{
    gatherImpl( allvals, n, root, myvals, sizes );
}

void MPICommunicator::gather(
    double allvals[],
    const IndexType n,
    const PartitionId root,
    const double myvals[],
    const IndexType sizes[] ) const
{
    gatherImpl( allvals, n, root, myvals, sizes );
}

void MPICommunicator::gather(
    int allvals[],
    const IndexType n,
    const PartitionId root,
    const int myvals[],
    const IndexType sizes[] ) const
{
    gatherImpl( allvals, n, root, myvals, sizes );
}

/* ---------------------------------------------------------------------------------- */
/*           maxloc                                                                   */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void MPICommunicator::maxlocImpl( T& val, int& location, PartitionId root ) const
{
    struct ValAndLoc
    {
        T val;
        int location;
    };
    ValAndLoc in;
    in.val = val;
    in.location = location;
    ValAndLoc out;
    MPI_Datatype commType = getMPI2Type<T,int>();
    MPI_Reduce( &in, &out, 1, commType, MPI_MAXLOC, root, selectMPIComm() );

    if ( mRank == root )
    {
        val = out.val;
        location = out.location;
    }
}

void MPICommunicator::maxloc( float& val, int& location, PartitionId root ) const
{
    maxlocImpl( val, location, root );
}

void MPICommunicator::maxloc( double& val, int& location, PartitionId root ) const
{
    maxlocImpl( val, location, root );
}

void MPICommunicator::maxloc( int& val, int& location, PartitionId root ) const
{
    maxlocImpl( val, location, root );
}

/* ---------------------------------------------------------------------------------- */
/*           swap                                                                     */
/* ---------------------------------------------------------------------------------- */

template<typename T>
void MPICommunicator::swapImpl( T val[], const IndexType n, PartitionId partner ) const
{
    if ( partner == mRank )
    {
        return;
    }

    boost::scoped_array<T> tmp( new T[n] );

    for ( IndexType i = 0; i < n; i++ )
    {
        tmp[i] = val[i];
    }

    MPI_Status mpiStatus;
    MPI_Datatype commType = getMPIType<T>();
    LAMA_MPICALL( logger,
                  MPI_Sendrecv( tmp.get(), n, commType, partner, defaultTag, val, n, commType, partner, defaultTag, selectMPIComm(), &mpiStatus ),
                  "MPI_Sendrecv" );
    LAMA_ASSERT_ERROR( getCount<T>( mpiStatus ) == n, "size mismatch for swap" )
}

void MPICommunicator::swap( double val[], const IndexType n, PartitionId partner ) const
{
    swapImpl( val, n, partner );
}

void MPICommunicator::swap( float val[], const IndexType n, PartitionId partner ) const
{
    swapImpl( val, n, partner );
}

void MPICommunicator::swap( int val[], const IndexType n, PartitionId partner ) const
{
    swapImpl( val, n, partner );
}

ContextPtr MPICommunicator::getCommunicationContext() const
{
    return ContextFactory::getContext( Context::Host );
}

/* ---------------------------------------------------------------------------------- */
/*           writeAt                                                                  */
/* ---------------------------------------------------------------------------------- */

void MPICommunicator::writeAt( std::ostream& stream ) const
{
    // Info about rank and size of the communicator is very useful
    stream << "MPI(" << mRank << ":" << mSize << ")";
}

} // namespace lama
