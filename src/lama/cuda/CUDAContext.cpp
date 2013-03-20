/**
 * @file CUDAContext.cpp
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
 * @file CUDAContext.cpp
 * @brief Contains the implementation of the class CUDAContext.
 * @author Thomas Brandes
 * Created on: 15.07.2011
 */

// hpp
#include <lama/cuda/CUDAContext.hpp>

// others
#include <lama/cuda/CUDAStreamSyncToken.hpp>
#include <lama/cuda/CUDAError.hpp>

#include <lama/task/TaskSyncToken.hpp>
#include <lama/task/Thread.hpp>

#include <lama/ContextAccess.hpp>

#include <lama/exception/LAMAAssert.hpp>

// tracing
#include <lama/tracing.hpp>

// boost
#include <boost/bind.hpp>

#include <cublas.h>
#include <cuda.h>
#include <cusparse.h>

#include <memory>

namespace lama
{

/**  static variables *****************************************************/

LAMA_LOG_DEF_LOGGER( CUDAContext::logger, "Context.CUDAContext" );

int CUDAContext::currentDeviceNr = -1;

int CUDAContext::numUsedDevices = 0;

//TODO: Determine good general default value
//TODO: Add setter routine for this
size_t CUDAContext::minPinnedSize = std::numeric_limits<size_t>::max();

cusparseHandle_t CUDAContext_cusparseHandle = 0;

/**  constructor  *********************************************************/

CUDAContext::CUDAContext( int deviceNr )
    : Context( CUDA ), mDeviceNr( deviceNr )
{
    LAMA_LOG_DEBUG( logger, "construct CUDAContext, device nr = = " << deviceNr );

    mNumberOfAllocatedBytes = 0;
    mNumberOfAllocates = 0;
    mMaxNumberOfAllocatedBytes = 0;

    // Note: logging is safe as CUDA context is always created after static initializations

    if( numUsedDevices == 0 )
    {
        unsigned int flags = 0; // must be set to zero
        LAMA_CUDA_DRV_CALL( cuInit( flags ), "cuInit failed for first used CUDA device" );
    }

    currentDeviceNr = deviceNr;

    numUsedDevices++;

    Context::enable( __FILE__, __LINE__ );

    LAMA_CUDA_DRV_CALL( cuDeviceGet( &mCUdevice, mDeviceNr ), "cuDeviceGet device " << mDeviceNr );

    char deviceName[256];

    LAMA_CUDA_DRV_CALL( cuDeviceGetName( deviceName, 256, mCUdevice ), "cuDeviceGetName" );

    mDeviceName = deviceName; // save it as string member variable for output

    LAMA_LOG_DEBUG( logger, "got device " << mDeviceName );

    LAMA_CUDA_DRV_CALL( cuCtxCreate( &mCUcontext, CU_CTX_SCHED_SPIN | CU_CTX_MAP_HOST, mCUdevice ),
                        "cuCtxCreate for " << *this );

    LAMA_CUBLAS_CALL( cublasInit(), "Initialization of CUBLAS library" );

    LAMA_CUSPARSE_CALL( cusparseCreate( &CUDAContext_cusparseHandle ),
                        "Initialization of CUSparse library: cusparseCreate" );

    LAMA_LOG_INFO( logger, "Initialized: CUBLAS, CuSparse" );

    int flags = 0; // must be 0 by specification of CUDA driver API

    LAMA_CUDA_DRV_CALL( cuStreamCreate( &mTransferStream, flags ), "cuStreamCreate for transfer failed" );
    LAMA_CUDA_DRV_CALL( cuStreamCreate( &mComputeStream, flags ), "cuStreamCreate for compute failed" );

    mOwnerThread = Thread::getSelf(); // thread that can use the context

    disable( __FILE__, __LINE__ );

    // count used devices to shutdown CUDA if no more device is used

    LAMA_LOG_INFO( logger, *this << " constructed by thread " << mOwnerThread << ", now disabled" );
}

/**  destructor   *********************************************************/

CUDAContext::~CUDAContext()
{
    LAMA_LOG_INFO( logger, "~CUDAContext: " << *this );

    if( mNumberOfAllocates > 0 )
    {
        LAMA_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " allocate without free" );
    }
    else if( mNumberOfAllocates < 0 )
    {
        LAMA_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " free without allocate" );
    }

    if( mNumberOfAllocatedBytes != 0 )
    {
        LAMA_LOG_ERROR( logger,
                        *this << ": number of allocated bytes = " << mNumberOfAllocatedBytes << ", mismatch of free/allocate sizes" );
    }

    numUsedDevices--;

    // context must be valid before streams are destroyed

    CUresult res = cuCtxPushCurrent( mCUcontext );

    if( res == CUDA_ERROR_DEINITIALIZED )
    {
        // this might happen with other software, e.g. VampirTrace

        LAMA_LOG_WARN( logger, "CUDA driver already deinitialized" );
        return;
    }
    else if( res != CUDA_SUCCESS )
    {
        LAMA_LOG_ERROR( logger, "Could not push any more context for " << *this );
        return;
    }

    LAMA_LOG_DEBUG( logger, "pushed context: synchronize/destroy compute stream" );

    LAMA_CUDA_DRV_CALL( cuStreamSynchronize( mComputeStream ), "cuStreamSynchronize for compute failed" );
    LAMA_CUDA_DRV_CALL( cuStreamDestroy( mComputeStream ), "cuStreamDestroy for compute failed" );

    LAMA_LOG_DEBUG( logger, "synchronize/destroy transfer stream" );

    LAMA_CUDA_DRV_CALL( cuStreamSynchronize( mTransferStream ), "cuStreamSynchronize for transfer failed" );
    LAMA_CUDA_DRV_CALL( cuStreamDestroy( mTransferStream ), "cuStreamDestroy for transfer failed" );

    if( !numUsedDevices )
    {
        LAMA_LOG_DEBUG( logger, "no more devices in use -> shutdown cuda" );

        cublasShutdown();
        LAMA_CHECK_CUBLAS_ERROR
        ;

//        lama_shutdown_cuda();
        if( CUDAContext_cusparseHandle )
        {
            cusparseStatus_t error = cusparseDestroy( CUDAContext_cusparseHandle );

            if( error != CUSPARSE_STATUS_SUCCESS )
            {
                LAMA_LOG_ERROR( logger, "Could not destroy cusparse handle, status = " << error );
            }

            LAMA_LOG_INFO( logger, "cusparse handle successfully destroyed" );

            CUDAContext_cusparseHandle = 0;
        }
    }

    LAMA_LOG_DEBUG( logger, "destroy cuda context" );

    // do not do it if CUDA tracing is enabled
    LAMA_CUDA_DRV_CALL( cuCtxDestroy( mCUcontext ), "cuCtxDestroy failed" );

    // if we are current device, set it back to -1

    if( currentDeviceNr == mDeviceNr )
    {
        currentDeviceNr = -1;
    }

    LAMA_LOG_INFO( logger, "Max allocated device memory " << mMaxNumberOfAllocatedBytes << " bytes." );
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::writeAt( std::ostream& stream ) const
{
    stream << "CUDAContext(" << mDeviceNr << ": " << mDeviceName << ")";
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::disable( const char* file, int line ) const
{
    /* not serious
     LAMA_ASSERT_ERROR( currentDeviceNr != -1, "disable context, but no device in use" );
     */

    Context::disable( file, line ); // call routine of base class

    LAMA_LOG_DEBUG( logger, *this << ": disable by thread = " << Thread::getSelf() );

    CUcontext tmp; // temporary for last context, not necessary to save it
    LAMA_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" );
    currentDeviceNr = -1;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::enable( const char* file, int line ) const
{
    LAMA_LOG_DEBUG( logger, *this << ": enable by thread = " << Thread::getSelf() );

    Context::enable( file, line ); // call routine of base class

    /* This is not serious:

     LAMA_ASSERT_ERROR( currentDeviceNr == -1, "enable context, but other device "
     << currentDeviceNr << " in use" );
     */

    LAMA_CUDA_DRV_CALL( cuCtxPushCurrent( mCUcontext ), "could not push context for " << *this );

//    lama_init0_cuda ();  // otherwise some properties of device might be wrong

    currentDeviceNr = mDeviceNr;
}

/* ----------------------------------------------------------------------------- */

bool CUDAContext::canUseData( const Context& other ) const
{
    // same object by pointer can always use same data.

    if( this == &other )
    {
        return true;
    }

    // CUDA device can use only data on same CUDA device

    if( other.getType() == CUDA )
    {
        const CUDAContext& otherCUDA = static_cast<const CUDAContext&>( other );
        return otherCUDA.mDeviceNr == mDeviceNr;
    }

    return false;
}

/* ----------------------------------------------------------------------------- */

void* CUDAContext::allocate( const size_t size ) const
{
    LAMA_REGION( "CUDAContext::allocate" );

    LAMA_CONTEXT_ACCESS( shared_from_this() );

    LAMA_ASSERT_ERROR( size > 0, "should not call allocate for size = " << size );

    LAMA_LOG_TRACE( logger, *this << ": allocate " << size << " bytes" );

    CUdeviceptr pointer = 0;

    LAMA_CUDA_DRV_CALL(
        cuMemAlloc( &pointer, size),
        "cuMemAlloc( " << pointer << ", " << size << " ) failed. This allocation would require a total of " << mMaxNumberOfAllocatedBytes + size << " bytes global memory." );

    LAMA_LOG_DEBUG( logger, *this << ": allocated " << size << " bytes, ptr = " << pointer );

    mNumberOfAllocatedBytes += size;
    mNumberOfAllocates++;

    mMaxNumberOfAllocatedBytes = std::max( mNumberOfAllocatedBytes, mMaxNumberOfAllocatedBytes );

    return (void *) pointer;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::allocate( ContextData& contextData, const size_t size ) const
{
    contextData.pointer = allocate( size );
    contextData.setPinned();
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::free( void* pointer, const size_t size ) const
{
    LAMA_REGION( "CUDAContext::free" );

    LAMA_CONTEXT_ACCESS( shared_from_this() );

    LAMA_LOG_DEBUG( logger, *this << ": free " << size << " bytes, ptr = " << pointer );

    LAMA_CUDA_DRV_CALL( cuMemFree( (CUdeviceptr) pointer ), "cuMemFree( " << pointer << " ) failed" );

    mNumberOfAllocatedBytes -= size;
    mNumberOfAllocates--;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::free( ContextData& contextData ) const
{
    LAMA_ASSERT_EQUAL_ERROR( contextData.context->getType(), getType() );
    free( contextData.pointer, contextData.size );
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpy( void* dst, const void* src, const size_t size ) const
{
    LAMA_CONTEXT_ACCESS( shared_from_this() );

    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (device) " );

    LAMA_CUDA_DRV_CALL( cuMemcpyDtoD( (CUdeviceptr) dst, (CUdeviceptr) src, size ),
                        "cuMemcpyDtoD( " << dst << ", " << src << ", " << size << " ) failed" );
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<SyncToken> CUDAContext::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    LAMA_CONTEXT_ACCESS( shared_from_this() );

    CUDAStreamSyncTokenPtr syncToken( new CUDAStreamSyncToken( shared_from_this(), mTransferStream ) );

    LAMA_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (device) to " << dst << " (device) " );

    LAMA_CUDA_DRV_CALL( cuMemcpyDtoDAsync( (CUdeviceptr) dst,(CUdeviceptr) src, size, mTransferStream ),
                        "cuMemcpyDtoDAsync( " << dst << ", " << src << ", " << size << ") failed " );

    return syncToken;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpyFromHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_CONTEXT_ACCESS( shared_from_this() );

    LAMA_REGION( "memcpyHost2Cuda" );

    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (host) to " << dst << " (device) " );

    LAMA_CUDA_DRV_CALL( cuMemcpyHtoD( (CUdeviceptr) dst, src, size),
                        "cuMemcpyHToD( " << dst << ", " << src << ", " << size << ") failed " );
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<SyncToken> CUDAContext::memcpyAsyncFromHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_LOG_INFO( logger, "async copy " << size << " bytes from " << src << " (host) to " << dst << " (device) " );

    // as current thread has disabled the context, another thread might use it
    return std::auto_ptr<SyncToken>(
               new TaskSyncToken( boost::bind( &CUDAContext::memcpyFromHost, this, dst, src, size ) ) );
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpyToHost( void* dst, const void* src, const size_t size ) const
{
    //LAMA_THROWEXCEPTION("memcpyToHost is temporarily forbidden for tracking");
    LAMA_CONTEXT_ACCESS( shared_from_this() );

    LAMA_REGION("memcpyCuda2Host");

    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (host) " );

    LAMA_CUDA_DRV_CALL( cuMemcpyDtoH( dst, (CUdeviceptr) src, size),
                        "cuMemcpyDToH( " << dst << ", " << src << ", " << size << ") failed " );
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<SyncToken> CUDAContext::memcpyAsyncToHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_LOG_INFO( logger, "async copy " << size << " bytes from " << src << " (device) to " << dst << " (host) " );

    // as current thread has disabled the context, another thread might use it

    return std::auto_ptr<SyncToken>(
               new TaskSyncToken( boost::bind( &CUDAContext::memcpyToHost, this, dst, src, size ) ) );
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpyFromCUDAHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_REGION("memcpyCudaHost2CudaDev");

    LAMA_CONTEXT_ACCESS( shared_from_this() );

    LAMA_LOG_INFO( logger, "copy " << size <<" bytes from " << src << " (host) to " << dst << " (device) " );

    LAMA_CUDA_DRV_CALL( cuMemcpyHtoD( (CUdeviceptr) dst, src, size),
                        "cuMemcpyHToD( " << dst << ", " << src << ", " << size << ") failed " );
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<SyncToken> CUDAContext::memcpyAsyncFromCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // Alternative solution for comparison:
    //   return std::auto_ptr<SyncToken>(new TaskSyncToken( boost::bind( &CUDAContext::memcpyFromCUDAHost, this, dst, src, size ) ) );

    LAMA_CONTEXT_ACCESS( shared_from_this() );

    CUDAStreamSyncTokenPtr syncToken( new CUDAStreamSyncToken( shared_from_this(), mTransferStream ) );

    LAMA_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (host) to " << dst << " (device) " );

    LAMA_CUDA_DRV_CALL( cuMemcpyHtoDAsync( (CUdeviceptr) dst, src, size, mTransferStream ),
                        "cuMemcpyHtoDAsync( " << dst << ", " << src << ", " << size << ") failed " );

    return syncToken;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpyToCUDAHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_REGION("memcpyCudaDev2CudaHost");

    LAMA_CONTEXT_ACCESS( shared_from_this() );

    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (cuda host) " );

    LAMA_CUDA_DRV_CALL( cuMemcpyDtoH( dst, (CUdeviceptr) src, size),
                        "cuMemcpyDToH( " << dst << ", " << src << ", " << size << ") failed " );
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<SyncToken> CUDAContext::memcpyAsyncToCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // return std::auto_ptr<SyncToken>(new TaskSyncToken( boost::bind( &CUDAContext::memcpyToCUDAHost, this, dst, src, size ) ) );

    LAMA_CONTEXT_ACCESS( shared_from_this() );

    CUDAStreamSyncTokenPtr syncToken( new CUDAStreamSyncToken( shared_from_this(), mTransferStream ) );

    LAMA_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (device) to " << dst << " (host) " );

    LAMA_CUDA_DRV_CALL(
        cuMemcpyDtoHAsync( dst, (CUdeviceptr) src, size, mTransferStream ),
        "cuMemcpyDtoHAsync( " << dst << ", " << src << ", " << size << ", " << mTransferStream << ") failed " );

    return syncToken;
}

/* ----------------------------------------------------------------------------- */

bool CUDAContext::cancpy( const ContextData& dst, const ContextData& src ) const
{
    const Context::ContextType dstType = dst.context->getType();
    const Context::ContextType srcType = src.context->getType();

    return ( srcType == Host && dstType == CUDA ) || ( srcType == CUDA && dstType == Host )
           || ( srcType == CUDA && dstType == CUDA );
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpy( ContextData& dst, const ContextData& src, const size_t size ) const
{
    LAMA_ASSERT_ERROR( cancpy(dst,src), "Can not copy from "<< *(src.context) << " to " << *(dst.context) );
    const Context::ContextType dstType = dst.context->getType();
    const Context::ContextType srcType = src.context->getType();

    if( srcType == Host && dstType == CUDA )
    {
        if( !src.isPinned() && size > minPinnedSize )
        {
            LAMA_REGION("RegisterHostMemory");
            LAMA_CONTEXT_ACCESS( shared_from_this() );
            LAMA_CUDA_DRV_CALL( cuMemHostRegister( src.pointer, src.size,0),
                                "cuMemHostRegister( " << src.pointer << ", " << src.size << ", " << 0 << ") failed " );
            src.setPinned();
            src.setCleanFunction( cuMemHostUnregister );
        }

        if( src.isPinned() )
        {
            memcpyFromCUDAHost( dst.pointer, src.pointer, size );
        }
        else
        {
            memcpyFromHost( dst.pointer, src.pointer, size );
        }
    }
    else if( srcType == CUDA && dstType == Host )
    {
        if( !dst.isPinned() && size > minPinnedSize )
        {
            LAMA_REGION("RegisterHostMemory");
            LAMA_CONTEXT_ACCESS( shared_from_this() );
            LAMA_CUDA_DRV_CALL( cuMemHostRegister( dst.pointer, dst.size,0),
                                "cuMemHostRegister( " << src.pointer << ", " << dst.size << ", " << 0 << ") failed " );
            dst.setPinned();
            dst.setCleanFunction( cuMemHostUnregister );
        }
        if( dst.isPinned() )
        {
            memcpyToCUDAHost( dst.pointer, src.pointer, size );
        }
        else
        {
            memcpyToHost( dst.pointer, src.pointer, size );
        }
    }
    else if( srcType == CUDA && dstType == CUDA )
    {
//Compile it only if sufficient CUDA Version found.
#if CUDA_VERSION >= 4000
        //Both Contexts are cuda. So we can cast safely to CUDAContext.
        const CUDAContext* dstContext = static_cast<const CUDAContext*>( dst.context.get() );
        const CUDAContext* srcContext = static_cast<const CUDAContext*>( src.context.get() );

        if( dstContext->getDeviceNr() != srcContext->getDeviceNr() )
        {

            //Check for the access capability
            int accessCapability = 0;
            LAMA_CUDA_DRV_CALL(
                cuDeviceCanAccessPeer(&accessCapability, srcContext->getDeviceNr(),dstContext->getDeviceNr()),
                "cuDeviceCanAccessPeer failed" );

            //Find the Other Context for enabling peer access
            const CUDAContext* otherContext = srcContext;
            if( this->getDeviceNr() == srcContext->getDeviceNr() )
            {
                otherContext = dstContext;
            }

            //Activate Peer Access if possible
            if( accessCapability >= 1 )
            {
                CUresult resu = cuCtxEnablePeerAccess( otherContext->mCUcontext, 0 );
                if( resu != CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED && resu != CUDA_SUCCESS )
                {
                    LAMA_CUDA_DRV_CALL( resu, "cuCtxEnablePeerAccess failed" );
                }
            }
        }
#endif /*__CUDA_API_VERSION >= 4000*/

        //Copy as normal
        memcpy( dst.pointer, src.pointer, size );
    }
    else
    {
        LAMA_THROWEXCEPTION( "Can not copy from "<< *(src.context) << " to " << *(dst.context) );
    }
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<SyncToken> CUDAContext::memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const
{
    LAMA_LOG_INFO( logger, "memcpyAsync from " << *src.context << " to " << *dst.context << ", size = " << size );

    LAMA_ASSERT_ERROR( cancpy(dst,src), "Can not copy from " << *src.context << " to " << *dst.context );

    const Context::ContextType dstType = dst.context->getType();
    const Context::ContextType srcType = src.context->getType();

    if( srcType == Host && dstType == CUDA )
    {
        if( !src.isPinned() && size > minPinnedSize )
        {
            LAMA_REGION("RegisterHostMemory");
            LAMA_CONTEXT_ACCESS( shared_from_this() );
            LAMA_LOG_DEBUG( logger, "register host memory, size = " << src.size );
            LAMA_CUDA_DRV_CALL( cuMemHostRegister( src.pointer, src.size,0),
                                "cuMemHostRegister( " << src.pointer << ", " << src.size << ", " << 0 << ") failed " );
            src.setPinned();
            src.setCleanFunction( cuMemHostUnregister );
        }

        if( src.isPinned() )
        {
            return memcpyAsyncFromCUDAHost( dst.pointer, src.pointer, size );
        }
        else
        {
            return memcpyAsyncFromHost( dst.pointer, src.pointer, size );
        }
    }
    else if( srcType == CUDA && dstType == Host )
    {
        if( !dst.isPinned() && size > minPinnedSize )
        {
            LAMA_REGION("RegisterHostMemory");
            LAMA_CONTEXT_ACCESS( shared_from_this() );
            LAMA_CUDA_DRV_CALL( cuMemHostRegister( dst.pointer, dst.size,0),
                                "cuMemHostRegister( " << src.pointer << ", " << dst.size << ", " << 0 << ") failed " );
            dst.setPinned();
            dst.setCleanFunction( cuMemHostUnregister );
        }
        if( dst.isPinned() )
        {
            return memcpyAsyncToCUDAHost( dst.pointer, src.pointer, size );
        }
        else
        {
            return memcpyAsyncToHost( dst.pointer, src.pointer, size );
        }
    }
    else if( srcType == CUDA && dstType == CUDA )
    {
        return memcpyAsync( dst.pointer, src.pointer, size );
    }
    else
    {
        LAMA_THROWEXCEPTION( "Can not copy from "<< *(src.context) << " to " << *(dst.context) );
    }
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<CUDAStreamSyncToken> CUDAContext::getComputeSyncToken() const
{
    return std::auto_ptr<CUDAStreamSyncToken>( new CUDAStreamSyncToken( shared_from_this(), mComputeStream ) );
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<SyncToken> CUDAContext::getSyncToken() const
{
    return std::auto_ptr<SyncToken>( new CUDAStreamSyncToken( shared_from_this(), mComputeStream ) );
}

/* ----------------------------------------------------------------------------- */

std::auto_ptr<CUDAStreamSyncToken> CUDAContext::getTransferSyncToken() const
{
    return std::auto_ptr<CUDAStreamSyncToken>( new CUDAStreamSyncToken( shared_from_this(), mTransferStream ) );
}

/* ----------------------------------------------------------------------------- */

} //namespace lama

