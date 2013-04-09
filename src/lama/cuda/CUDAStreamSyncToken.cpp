/**
 * @file CUDAStreamSyncToken.cpp
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
 * @brief CUDAStreamSyncToken.cpp
 * @author Jiri Kraus
 * @date 20.05.2011
 * $Id$
 */

// hpp
#include <lama/cuda/CUDAStreamSyncToken.hpp>

// others
#include <lama/cuda/CUDAContext.hpp>
#include <lama/cuda/CUDAError.hpp>

#include <lama/ContextAccess.hpp>
#include <lama/tracing.hpp>

namespace lama
{

CUDAStreamSyncToken::CUDAStreamSyncToken( CUDAContextPtr cudaContext, CUstream stream )
    : mCUDAContext( cudaContext ), mStream( stream ), mEvent( 0 )
{
    LAMA_LOG_DEBUG( logger, "StreamSyncToken for " << *cudaContext << " generated, stream = " << stream )
}

CUDAStreamSyncToken::~CUDAStreamSyncToken()
{
    wait();
}

void CUDAStreamSyncToken::wait()
{
    if ( isSynchronized() )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "wait on CUDA stream synchronization" )

    {
        LAMA_REGION( "CUDA_synchronize" )

        LAMA_CONTEXT_ACCESS( mCUDAContext )

        if ( mEvent != 0 )
        {
            LAMA_CUDA_DRV_CALL( cuEventSynchronize( mEvent ), "cuEventSynchronize( " << mEvent << " ) failed." )
            LAMA_CUDA_DRV_CALL( cuEventDestroy( mEvent ), "cuEventDestroy( " << mEvent << " ) failed." );
            mEvent = 0;
        }
        else
        {
            LAMA_CUDA_DRV_CALL( cuStreamSynchronize( mStream ), "cuStreamSynchronize( " << mStream <<" ) failed." );
        }
    }

    setSynchronized();
}

bool CUDAStreamSyncToken::probe() const
{
    if ( isSynchronized() )
    {
        return true;
    }
    return probeEvent( mEvent );
}

cudaStream_t CUDAStreamSyncToken::getCUDAStream() const
{
    return (cudaStream_t) mStream;
}

void CUDAStreamSyncToken::createEvent( CUevent& event ) const
{
    LAMA_CONTEXT_ACCESS( mCUDAContext )
    LAMA_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ),
                        "Could not create event " )
}

void CUDAStreamSyncToken::createTimingEvent( CUevent& event ) const
{
    LAMA_CONTEXT_ACCESS( mCUDAContext )
    LAMA_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT ), "Could not create event " )
}

void CUDAStreamSyncToken::getTime( float* time, CUevent& startEvent, CUevent& stopEvent ) const
{
    LAMA_CONTEXT_ACCESS( mCUDAContext )
    LAMA_CUDA_DRV_CALL( cuEventElapsedTime( time, startEvent, stopEvent ),
                        "cuEventElapsedTime failed for CUevent " << stopEvent << '.' )
}

void CUDAStreamSyncToken::recordEvent( const CUevent event )
{
    LAMA_CONTEXT_ACCESS( mCUDAContext )
    LAMA_CUDA_DRV_CALL( cuEventRecord ( event, mStream ), "cuEventRecord failed for CUevent "<<event<<'.' )
}

bool CUDAStreamSyncToken::probeEvent( const CUevent& stopEvent ) const
{
    LAMA_ASSERT_ERROR( stopEvent != 0, "probe on invalid event" )

    LAMA_CONTEXT_ACCESS( mCUDAContext )

    CUresult result = cuEventQuery( stopEvent );

    if ( result != CUDA_SUCCESS && result != CUDA_ERROR_NOT_READY )
    {
        LAMA_CUDA_DRV_CALL( result, "cuEventQuery failed for CUevent " << stopEvent << '.' );
    }

    return result == CUDA_SUCCESS;
}

bool CUDAStreamSyncToken::queryEvent( const CUevent event ) const
{
    LAMA_CONTEXT_ACCESS( mCUDAContext )

    CUresult result = cuEventQuery( event );

    if ( result != CUDA_SUCCESS || result != CUDA_ERROR_NOT_READY )
    {
        LAMA_CUDA_DRV_CALL( result, "cuEventQuery failed for CUevent "<<event<<'.' );
    }

    return result == CUDA_SUCCESS;
}

void CUDAStreamSyncToken::synchronizeEvent( const CUevent event ) const
{
    LAMA_CONTEXT_ACCESS( mCUDAContext )
    LAMA_CUDA_DRV_CALL( cuEventSynchronize( event ), "cuEventSynchronize failed for CUevent "<<event<<'.' )
}

CUDAStreamSyncTokenPtr::CUDAStreamSyncTokenPtr() throw ()
{
}

CUDAStreamSyncTokenPtr::CUDAStreamSyncTokenPtr( CUDAStreamSyncTokenPtr& other ) throw ()
    : mCUDAStreamSyncToken( other.mCUDAStreamSyncToken )
{
}

CUDAStreamSyncTokenPtr::CUDAStreamSyncTokenPtr( std::auto_ptr<CUDAStreamSyncToken> other ) throw ()
    : mCUDAStreamSyncToken( other )
{
}

CUDAStreamSyncTokenPtr::CUDAStreamSyncTokenPtr( CUDAStreamSyncToken* pointer ) throw ()
    : mCUDAStreamSyncToken( pointer )
{
}

CUDAStreamSyncTokenPtr::~CUDAStreamSyncTokenPtr()
{
}

CUDAStreamSyncTokenPtr& CUDAStreamSyncTokenPtr::operator=( std::auto_ptr<CUDAStreamSyncToken> other ) throw ()
{
    mCUDAStreamSyncToken = other;
    return *this;
}

CUDAStreamSyncTokenPtr& CUDAStreamSyncTokenPtr::operator=( CUDAStreamSyncTokenPtr& other ) throw ()
{
    mCUDAStreamSyncToken = other.mCUDAStreamSyncToken;
    return *this;
}

CUDAStreamSyncTokenPtr::operator std::auto_ptr<SyncToken>&()
{
    LAMA_ASSERT_ERROR( mCUDAStreamSyncToken->mEvent == 0, "Tried to create a already created event." )
    mCUDAStreamSyncToken->createEvent( mCUDAStreamSyncToken->mEvent );
    mCUDAStreamSyncToken->recordEvent( mCUDAStreamSyncToken->mEvent );
    mSyncToken = mCUDAStreamSyncToken;
    return mSyncToken;
}

CUDAStreamSyncToken* CUDAStreamSyncTokenPtr::operator->() const throw ()
{
    LAMA_ASSERT_ERROR( mCUDAStreamSyncToken.get() != 0, "Tried to dereference a NULL Pointer." )
    return mCUDAStreamSyncToken.get();
}

CUDAStreamSyncToken& CUDAStreamSyncTokenPtr::operator*() const throw ()
{
    LAMA_ASSERT_ERROR( mCUDAStreamSyncToken.get() != 0, "Tried to dereference a NULL Pointer." )
    return *mCUDAStreamSyncToken;
}

CUDAStreamSyncToken* CUDAStreamSyncTokenPtr::release() throw ()
{
    return mCUDAStreamSyncToken.release();
}

}

