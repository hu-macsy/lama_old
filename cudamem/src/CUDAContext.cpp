/**
 * @file CUDAContext.cpp
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
 * @brief Contains the implementation of the class CUDAContext.
 * @author Thomas Brandes
 * @date 15.07.2011
 * @since 1.0.0
 */

// hpp
#include <cudamem/CUDAContext.hpp>
#include <cudamem/CUDAHostContext.hpp>
#include <cudamem/CUDAError.hpp>

// others
#include <cudamem/CUDAStreamSyncToken.hpp>

#include <tasking/TaskSyncToken.hpp>
#include <memory/ContextAccess.hpp>

// boost
#include <boost/bind.hpp>

#include <cublas_v2.h>
#include <cuda.h>
#include <cusparse.h>

#include <memory>

using common::Thread;
using tasking::SyncToken;

namespace memory
{

/**  static variables *****************************************************/

LAMA_LOG_DEF_LOGGER( CUDAContext::logger, "Context.CUDAContext" )

int CUDAContext::currentDeviceNr = -1;

int CUDAContext::numUsedDevices = 0;

//TODO: Determine good general default value
//TODO: Add setter routine for this
size_t CUDAContext::minPinnedSize = std::numeric_limits<size_t>::max();

cusparseHandle_t CUDAContext_cusparseHandle = 0;
cublasHandle_t CUDAContext_cublasHandle = 0;

/**  constructor  *********************************************************/

CUDAContext::CUDAContext( int deviceNr )
    : Context( context::CUDA ), mDeviceNr( deviceNr )
{
    LAMA_LOG_DEBUG( logger, "construct CUDAContext, device nr = = " << deviceNr )
    mNumberOfAllocatedBytes = 0;
    mNumberOfAllocates = 0;
    mMaxNumberOfAllocatedBytes = 0;

    // Note: logging is safe as CUDA context is always created after static initializations

    if ( numUsedDevices == 0 )
    {
        unsigned int flags = 0; // must be set to zero
        LAMA_CUDA_DRV_CALL( cuInit( flags ), "cuInit failed for first used CUDA device" )
    }

    currentDeviceNr = deviceNr;
    numUsedDevices++;
    Context::enable( __FILE__, __LINE__ );
    LAMA_CUDA_DRV_CALL( cuDeviceGet( &mCUdevice, mDeviceNr ), "cuDeviceGet device " << mDeviceNr );

    cudaDeviceProp properties;
    cudaGetDeviceProperties( &properties, mCUdevice );

    LAMA_LOG_ERROR( logger, "canMapHostMemory = " << properties.canMapHostMemory );


    char deviceName[256];
    LAMA_CUDA_DRV_CALL( cuDeviceGetName( deviceName, 256, mCUdevice ), "cuDeviceGetName" );
    mDeviceName = deviceName; // save it as string member variable for output
    LAMA_LOG_DEBUG( logger, "got device " << mDeviceName )
    LAMA_CUDA_DRV_CALL( cuCtxCreate( &mCUcontext, CU_CTX_SCHED_SPIN | CU_CTX_MAP_HOST, mCUdevice ),
                        "cuCtxCreate for " << *this )
    LAMA_CUSPARSE_CALL( cusparseCreate( &CUDAContext_cusparseHandle ),
                        "Initialization of CUSparse library: cusparseCreate" );
    LAMA_CUBLAS_CALL( cublasCreate( &CUDAContext_cublasHandle ), "Initialization of CUBlas library: cublasCreate" );
    LAMA_LOG_INFO( logger, "Initialized: CUBLAS, CuSparse" )
    int flags = 0; // must be 0 by specification of CUDA driver API
    LAMA_CUDA_DRV_CALL( cuStreamCreate( &mTransferStream, flags ), "cuStreamCreate for transfer failed" )
    LAMA_CUDA_DRV_CALL( cuStreamCreate( &mComputeStream, flags ), "cuStreamCreate for compute failed" );
    mOwnerThread = Thread::getSelf(); // thread that can use the context
    disable( __FILE__, __LINE__ );
    // count used devices to shutdown CUDA if no more device is used
    LAMA_LOG_INFO( logger, *this << " constructed by thread " << mOwnerThread << ", now disabled" )
}

/**  destructor   *********************************************************/

CUDAContext::~CUDAContext()
{
    LAMA_LOG_INFO( logger, "~CUDAContext: " << *this )

    if ( mNumberOfAllocates > 0 )
    {
        LAMA_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " allocate without free" )
    }
    else if ( mNumberOfAllocates < 0 )
    {
        LAMA_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " free without allocate" )
    }

    if ( mNumberOfAllocatedBytes != 0 )
    {
        LAMA_LOG_ERROR( logger,
                        *this << ": number of allocated bytes = " << mNumberOfAllocatedBytes << ", mismatch of free/allocate sizes" )
    }

    numUsedDevices--;
    // context must be valid before streams are destroyed
    CUresult res = cuCtxPushCurrent( mCUcontext );

    if ( res == CUDA_ERROR_DEINITIALIZED )
    {
        // this might happen with other software, e.g. VampirTrace
        LAMA_LOG_WARN( logger, "CUDA driver already deinitialized" )
        return;
    }
    else if ( res != CUDA_SUCCESS )
    {
        LAMA_LOG_ERROR( logger, "Could not push any more context for " << *this )
        return;
    }

    LAMA_LOG_DEBUG( logger, "pushed context: synchronize/destroy compute stream" )
    LAMA_CUDA_DRV_CALL( cuStreamSynchronize( mComputeStream ), "cuStreamSynchronize for compute failed" )
    LAMA_CUDA_DRV_CALL( cuStreamDestroy( mComputeStream ), "cuStreamDestroy for compute failed" );
    LAMA_LOG_DEBUG( logger, "synchronize/destroy transfer stream" )
    LAMA_CUDA_DRV_CALL( cuStreamSynchronize( mTransferStream ), "cuStreamSynchronize for transfer failed" )
    LAMA_CUDA_DRV_CALL( cuStreamDestroy( mTransferStream ), "cuStreamDestroy for transfer failed" );

    if ( !numUsedDevices )
    {
        LAMA_LOG_DEBUG( logger, "no more devices in use -> shutdown cuda" )

        if ( CUDAContext_cublasHandle )
        {
            cublasStatus_t error = cublasDestroy( CUDAContext_cublasHandle );

            if ( error != CUBLAS_STATUS_SUCCESS )
            {
                LAMA_LOG_ERROR( logger, "Could not destroy cublas handle, status = " << error )
            }

            LAMA_LOG_INFO( logger, "cublas handle successfully destroyed" )
            CUDAContext_cublasHandle = 0;
        }

//        cublasShutdown();

//        lama_shutdown_cuda();
        if ( CUDAContext_cusparseHandle )
        {
            cusparseStatus_t error = cusparseDestroy( CUDAContext_cusparseHandle );

            if ( error != CUSPARSE_STATUS_SUCCESS )
            {
                LAMA_LOG_ERROR( logger, "Could not destroy cusparse handle, status = " << error )
            }

            LAMA_LOG_INFO( logger, "cusparse handle successfully destroyed" )
            CUDAContext_cusparseHandle = 0;
        }
    }

    LAMA_LOG_DEBUG( logger, "destroy cuda context" )
    // do not do it if CUDA tracing is enabled
    LAMA_CUDA_DRV_CALL( cuCtxDestroy( mCUcontext ), "cuCtxDestroy failed" )

    // if we are current device, set it back to -1

    if ( currentDeviceNr == mDeviceNr )
    {
        currentDeviceNr = -1;
    }

    LAMA_LOG_INFO( logger, "Max allocated device memory " << mMaxNumberOfAllocatedBytes << " bytes." )
}

/* ----------------------------------------------------------------------------- */

ContextPtr CUDAContext::getHostContext() const
{
    ContextPtr context;

    if ( mHostContext.expired() )
    {
        context = Context::getContext( context::CUDAHost, mDeviceNr );
        mHostContext = context; // save it here as a weak pointer to avoid cycles
    }
    else
    {
        // the last host context instance is still valid, so we return just shared pointer to it
        context = mHostContext.lock();
    }

    return context;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::writeAt( std::ostream& stream ) const
{
    stream << "CUDAContext(" << mDeviceNr << ": " << mDeviceName << ")";
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::disable( const char* file, int line ) const
{
    Context::disable( file, line ); // call routine of base class
    LAMA_LOG_DEBUG( logger, *this << ": disable by thread = " << Thread::getSelf() )
    CUcontext tmp; // temporary for last context, not necessary to save it
    LAMA_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" )
    currentDeviceNr = -1;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::enable( const char* file, int line ) const
{
    LAMA_LOG_DEBUG( logger, *this << ": enable by thread = " << Thread::getSelf() )
    Context::enable( file, line ); // call routine of base class
    LAMA_CUDA_DRV_CALL( cuCtxPushCurrent( mCUcontext ), "could not push context for " << *this )
//    lama_init0_cuda ();  // otherwise some properties of device might be wrong
    currentDeviceNr = mDeviceNr;
}

/* ----------------------------------------------------------------------------- */

bool CUDAContext::canUseData( const Context& other ) const
{
    // same object by pointer can always use same data.
    if ( this == &other )
    {
        return true;
    }

    // CUDA device can use only data on same CUDA device

    if ( other.getType() == context::CUDA )
    {
        const CUDAContext& otherCUDA = static_cast<const CUDAContext&>( other );
        return otherCUDA.mDeviceNr == mDeviceNr;
    }

    return false;
}

/* ----------------------------------------------------------------------------- */

void* CUDAContext::allocate( const size_t size ) const
{
    // LAMA_REGION( "CUDA.allocate" )
    LAMA_CONTEXT_ACCESS( shared_from_this() )
    COMMON_ASSERT( size > 0, "should not call allocate for size = " << size )
    LAMA_LOG_TRACE( logger, *this << ": allocate " << size << " bytes" )
    CUdeviceptr pointer = 0;
    LAMA_CUDA_DRV_CALL(
        cuMemAlloc( &pointer, size ),
        "cuMemAlloc( size = " << size << " ) failed. This allocation would require a total of " << mMaxNumberOfAllocatedBytes + size << " bytes global memory." )
    LAMA_LOG_DEBUG( logger, *this << ": allocated " << size << " bytes, ptr = " << ( ( void* ) pointer ) )
    mNumberOfAllocatedBytes += size;
    mNumberOfAllocates++;
    mMaxNumberOfAllocatedBytes = std::max( mNumberOfAllocatedBytes, mMaxNumberOfAllocatedBytes );
    return ( void* ) pointer;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::free( void* pointer, const size_t size ) const
{
    // LAMA_REGION( "CUDA.free" )
    LAMA_CONTEXT_ACCESS( shared_from_this() )
    LAMA_LOG_DEBUG( logger, *this << ": free " << size << " bytes, ptr = " << pointer )
    LAMA_CUDA_DRV_CALL( cuMemFree( ( CUdeviceptr ) pointer ), "cuMemFree( " << pointer << " ) failed" )
    mNumberOfAllocatedBytes -= size;
    mNumberOfAllocates--;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpy( void* dst, const void* src, const size_t size ) const
{
    LAMA_CONTEXT_ACCESS( shared_from_this() )
    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (device) " )
    LAMA_CUDA_DRV_CALL( cuMemcpyDtoD( ( CUdeviceptr ) dst, ( CUdeviceptr ) src, size ),
                        "cuMemcpyDtoD( " << dst << ", " << src << ", " << size << " ) failed" )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAContext::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    // LAMA_REGION( "CUDA.memcpyDtoDAsync" )
    LAMA_CONTEXT_ACCESS( shared_from_this() )
    // use auto pointer so memory will be freed in case of exceptions
    LAMA_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (device) to " << dst << " (device) " )
    LAMA_CUDA_DRV_CALL( cuMemcpyDtoDAsync( ( CUdeviceptr ) dst, ( CUdeviceptr ) src, size, mTransferStream ),
                        "cuMemcpyDtoDAsync( " << dst << ", " << src << ", " << size << ") failed " )
    // sync token should not synchronize on the full stream but only on the transfer, so add event
    CUevent event;
    LAMA_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ), "Could not create event " )
    LAMA_CUDA_DRV_CALL( cuEventRecord( event, mTransferStream ), "cuEventRecord failed for CUevent " << event << '.' )
    // This way timing might be added
    // CuEvent startEevent;
    // CuEvent stopEevent;
    // startEvent = CUDAStreamSyncToken::getStartEvent( shared_from_this(), mTransferStream )
    // LAMA_CUDA_DRV_CALL( cuEventCreate( &startEvent, CU_EVENT_DEFAULT )
    // LAMA_CUDA_DRV_CALL( cuEventRecord( startEvent, mTransferStream )
    // now cuMemcpyDToDAsync
    // stopEvent = CUDAStreamSyncToken::getStopEvent( shared_from_this(), mTransferStream )
    // LAMA_CUDA_DRV_CALL( cuEventCreate( &stopEvent, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ),
    // LAMA_CUDA_DRV_CALL( cuEventRecord( &stopEvent, mTransferStream )
    // new CUDAStreamSyncToken( shared_from_this(), mTransferStream, startEvent, stopEvent )
    // -> will synchronize on stopEvent, gets time as follows
    // LAMA_CUDA_DRV_CALL( cuEventElapsedTime( time, startEvent, stopEvent ),
    return new CUDAStreamSyncToken( shared_from_this(), mTransferStream, event );
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpyFromHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (host) to " << dst << " (device) " )

    LAMA_CONTEXT_ACCESS( shared_from_this() )

    LAMA_CUDA_DRV_CALL( cuMemcpyHtoD( ( CUdeviceptr ) dst, src, size ),
                        "cuMemcpyHToD( " << dst << ", " << src << ", " << size << ") failed " )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAContext::memcpyAsyncFromHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_LOG_INFO( logger, "async copy " << size << " bytes from " << src << " (host) to " << dst << " (device) " )

    // as current thread has disabled the context, another thread might use it

    memcpyFromHost( dst, src, size );

    return NULL;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpyToHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (host) " )

    LAMA_CONTEXT_ACCESS( shared_from_this() )

    LAMA_CUDA_DRV_CALL( cuMemcpyDtoH( dst, ( CUdeviceptr ) src, size ),
                        "cuMemcpyDToH( " << dst << ", " << src << ", " << size << ") failed " )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAContext::memcpyAsyncToHost( void* dst, const void* src, const size_t size ) const
{
    LAMA_LOG_INFO( logger, "async copy " << size << " bytes from " << src << " (device) to " << dst << " (host) " )

    // as current thread has disabled the context, another thread might use it

    memcpyToHost( dst, src, size );

    return NULL;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpyFromCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // LAMA_REGION( "CUDA.memcpyCUDAHost->Dev")
    LAMA_CONTEXT_ACCESS( shared_from_this() )
    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (host) to " << dst << " (device) " )
    LAMA_CUDA_DRV_CALL( cuMemcpyHtoD( ( CUdeviceptr ) dst, src, size ),
                        "cuMemcpyHToD( " << dst << ", " << src << ", " << size << ") failed " )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAContext::memcpyAsyncFromCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // LAMA_REGION( "CUDA.memcpyHtoDAsync" )
    LAMA_CONTEXT_ACCESS( shared_from_this() )
    LAMA_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (host) to " << dst << " (device) " )
    LAMA_CUDA_DRV_CALL( cuMemcpyHtoDAsync( ( CUdeviceptr ) dst, src, size, mTransferStream ),
                        "cuMemcpyHtoDAsync( " << dst << ", " << src << ", " << size << ") failed " )
    return new CUDAStreamSyncToken( shared_from_this(), mTransferStream );
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::memcpyToCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // LAMA_REGION( "CUDA.memcpyDev->CUDAHost" )
    LAMA_CONTEXT_ACCESS( shared_from_this() )
    LAMA_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (cuda host) " )
    LAMA_CUDA_DRV_CALL( cuMemcpyDtoH( dst, ( CUdeviceptr ) src, size ),
                        "cuMemcpyDToH( " << dst << ", " << src << ", " << size << ") failed " )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAContext::memcpyAsyncToCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // LAMA_REGION( "CUDA.memcpyDtoHAsync" )
    LAMA_CONTEXT_ACCESS( shared_from_this() )
    LAMA_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (device) to " << dst << " (host) " )
    LAMA_CUDA_DRV_CALL(
        cuMemcpyDtoHAsync( dst, ( CUdeviceptr ) src, size, mTransferStream ),
        "cuMemcpyDtoHAsync( " << dst << ", " << src << ", " << size << ", " << mTransferStream << ") failed " )
    // sync token should not synchronize on the full stream but only on the transfer, so add event
    CUevent event;
    LAMA_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ), "Could not create event " )
    LAMA_CUDA_DRV_CALL( cuEventRecord( event, mTransferStream ), "cuEventRecord failed for CUevent " << event << '.' )
    return new CUDAStreamSyncToken( shared_from_this(), mTransferStream, event );
}

/* ----------------------------------------------------------------------------- */

CUDAStreamSyncToken* CUDAContext::getComputeSyncToken() const
{
    return new CUDAStreamSyncToken( shared_from_this(), mComputeStream );
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAContext::getSyncToken() const
{
    return new CUDAStreamSyncToken( shared_from_this(), mComputeStream );
}

/* ----------------------------------------------------------------------------- */

CUDAStreamSyncToken* CUDAContext::getTransferSyncToken() const
{
    return new CUDAStreamSyncToken( shared_from_this(), mTransferStream );
}

/* ----------------------------------------------------------------------------- */

bool CUDAContext::canCopyFrom( const Context& other ) const
{
    // copy from host to this context should always be supported

    ContextType otherType = other.getType();

    return ( otherType == context::Host ) || ( otherType == context::CUDAHost );
}

bool CUDAContext::canCopyTo( const Context& other ) const
{
    // copy from this context to host should always be supported

    ContextType otherType = other.getType();

    return ( otherType == context::Host ) || ( otherType == context::CUDAHost );
}

void CUDAContext::memcpyFrom( void* dst, const Context& srcContext, const void* src, size_t size ) const
{
    if ( srcContext.getType() == context::Host )
    {
        memcpyFromHost( dst, src, size );
    }
    else if ( srcContext.getType() == context::CUDAHost )
    {
        memcpyFromCUDAHost( dst, src, size );
    }
    else
    {
        COMMON_THROWEXCEPTION( "copy from " << srcContext << " to " << *this << " not supported" )
    }
}

void CUDAContext::memcpyTo( const Context& dstContext, void* dst, const void* src, size_t size ) const
{
    if ( dstContext.getType() == context::Host )
    {
        memcpyToHost( dst, src, size );
    }
    else if ( dstContext.getType() == context::CUDAHost )
    {
        memcpyToCUDAHost( dst, src, size );
    }
    else
    {
        COMMON_THROWEXCEPTION( "copy to " << dstContext << " from " << *this << " not supported" )
    }
}

/* ----------------------------------------------------------------------------- */

#define LAMA_DEFAULT_DEVICE_NUMBER -1
#define LAMA_MAX_CUDA_DEVICES 4

static int getDefaultDeviceNr() 
{
    return 0;
}

/* ----------------------------------------------------------------------------- */
/*      Factory::Register - create( int )                                        */
/* ----------------------------------------------------------------------------- */

static boost::weak_ptr<CUDAContext> mCUDAContext[LAMA_MAX_CUDA_DEVICES];

ContextPtr CUDAContext::create( int deviceNr )
{
    int cudaDeviceNr = deviceNr;

    if( cudaDeviceNr == LAMA_DEFAULT_DEVICE_NUMBER )
    {
        cudaDeviceNr = getDefaultDeviceNr();

        // no need here to check for a good value
    }
    else
    {
        COMMON_ASSERT(
            0 <= cudaDeviceNr && cudaDeviceNr < LAMA_MAX_CUDA_DEVICES,
            "device = " << cudaDeviceNr << " out of range" << ", max supported device = " << LAMA_MAX_CUDA_DEVICES )
    }

    boost::shared_ptr<CUDAContext> context = boost::shared_ptr<CUDAContext>();

    if( mCUDAContext[cudaDeviceNr].expired() )
    {
        // create a new context for the device and return the shared pointer

        context = boost::shared_ptr<CUDAContext>( new CUDAContext( cudaDeviceNr ) );

        // we keep a weak pointer so that we can return

        mCUDAContext[cudaDeviceNr] = context;
    }
    else
    {
        // the weak pointer to the device is still okay, so return a shared pointer for it

        context = mCUDAContext[cudaDeviceNr].lock();
    }

    return context;
}

/* ----------------------------------------------------------------------------- */

} //namespace

