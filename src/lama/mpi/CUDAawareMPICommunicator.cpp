/**
 * @file CUDAawareMPICommunicator.cpp
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
 * @brief CUDAawareMPICommunicator.cpp
 * @author Lauretta Schubert
 * @date 11.03.2014
 * @since 1.0.1
 */

// hpp
#include <lama/mpi/CUDAawareMPICommunicator.hpp>

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

const int CUDAawareMPICommunicator::defaultTag = 1;

LAMA_LOG_DEF_LOGGER( CUDAawareMPICommunicator::logger, "Communicator.CUDAawareMPICommunicator" )

CUDAawareMPICommunicator::CUDAawareMPICommunicator( int& argc, char** & argv )
    : MPICommunicator( argc, argv, "CUDAawareMPI" )
{
    int initialized = 0;

    LAMA_MPICALL( logger, MPI_Initialized( &initialized ), "MPI_Initialized" )

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
    LAMA_MPICALL( logger, MPI_Comm_size( mComm, &mSize ), "MPI_Comm_size" )
    LAMA_MPICALL( logger, MPI_Comm_rank( mComm, &mRank ), "MPI_Comm_rank" )

    setNodeData();   // determine mNodeRank, mNodeSize
}

/* ---------------------------------------------------------------------------------- */

CUDAawareMPICommunicator::~CUDAawareMPICommunicator()
{
    LAMA_LOG_INFO( logger, *this << ": ~CUDAawareMPICommunicator" )
    int finalized = 0;
    LAMA_MPICALL( logger, MPI_Finalized( &finalized ), "MPI_Finalized" )

    if ( !finalized )
    {
        if ( !mExternInitialization )
        {
            LAMA_LOG_INFO( logger, "call MPI_Finalize" )
            LAMA_MPICALL( logger, MPI_Finalize(), "MPI_Finalize" )
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

/* -------------------------------------------------------------------------- */

template<typename T>
void CUDAawareMPICommunicator::shiftArray( LAMAArray<T>& recvArray, const LAMAArray<T>& sendArray, const int direction ) const
{
    LAMA_ASSERT_ERROR( &recvArray != &sendArray, "send and receive array are same, not allowed for shift" )

    if ( direction % getSize() == 0 )
    {
        // self assignment

        recvArray = sendArray;
        return;
    }

    ReadAccess<T> sendData( sendArray, sendArray.getValidContext() );
    IndexType numSendElems = sendData.size();

    // make recv array large enough to fit for send data

    WriteOnlyAccess<T> recvData( recvArray, recvArray.getValidContext(), numSendElems );

    // but we are able to receive even more data if array is large enough

    IndexType maxNumRecvElems = recvData.capacity();

    // For shifting of data we use the pure virtual methods implemened by each communicator

    IndexType numRecvElems = shiftImpl( recvData.get(), maxNumRecvElems, sendData.get(), numSendElems, direction );

    LAMA_LOG_INFO( logger,
                   "shift, direction = " << direction << ", sent " << numSendElems << ", recvd " << numRecvElems << "( max was " << maxNumRecvElems << ")" )

    recvData.resize( numRecvElems ); // take over the size
}

template LAMA_DLL_IMPORTEXPORT
void CUDAawareMPICommunicator::shiftArray( LAMAArray<float>& recvArray, const LAMAArray<float>& sendArray, const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
void CUDAawareMPICommunicator::shiftArray( LAMAArray<double>& recvArray, const LAMAArray<double>& sendArray, const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
void CUDAawareMPICommunicator::shiftArray( LAMAArray<int>& recvArray, const LAMAArray<int>& sendArray, const int direction ) const;

/* -------------------------------------------------------------------------- */

template<typename T>
SyncToken* CUDAawareMPICommunicator::shiftAsync(
    LAMAArray<T>& recvArray,
    const LAMAArray<T>& sendArray,
    const int direction ) const
{
    LAMA_ASSERT_ERROR( &recvArray != &sendArray, "send and receive array are same, not allowed for shift" )

    recvArray.clear(); // do not keep any old data, keep capacities

    boost::shared_ptr<WriteAccess<T> > recvData( new WriteAccess<T>( recvArray, recvArray.getValidContext() ) );
    boost::shared_ptr<ReadAccess<T> > sendData( new ReadAccess<T>( sendArray, sendArray.getValidContext() ) );

    IndexType numElems = sendData->size();

    recvData->resize( numElems ); // size should fit at least to keep own data

    // For shifting of data we use the pure virtual methods implemened by each communicator
    // Note: get is the method of the accesses and not of the auto_ptr

    SyncToken* syncToken = shiftAsyncImpl( recvData->get(), sendData->get(), numElems, direction );

    LAMA_ASSERT_DEBUG( syncToken, "NULL pointer for sync token" )

    // accesses are pushed in the sync token so they are freed after synchronization

    syncToken->pushAccess( sendData );
    syncToken->pushAccess( recvData );

    return syncToken;
}

template LAMA_DLL_IMPORTEXPORT
SyncToken* CUDAawareMPICommunicator::shiftAsync(
    LAMAArray<float>& recvArray,
    const LAMAArray<float>& sendArray,
    const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
SyncToken* CUDAawareMPICommunicator::shiftAsync(
    LAMAArray<double>& recvArray,
    const LAMAArray<double>& sendArray,
    const int direction ) const;

template LAMA_DLL_IMPORTEXPORT
SyncToken* CUDAawareMPICommunicator::shiftAsync(
    LAMAArray<int>& recvArray,
    const LAMAArray<int>& sendArray,
    const int direction ) const;

/* ---------------------------------------------------------------------------------- */
/*           writeAt                                                                  */
/* ---------------------------------------------------------------------------------- */

void CUDAawareMPICommunicator::writeAt( std::ostream& stream ) const
{
    // Info about rank and size of the communicator is very useful
    stream << "MPI(" << mRank << ":" << mSize << ")";
}

} // namespace lama
