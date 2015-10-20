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
 * @author Thomas Brandes, Jiri Kraus
 * @date 15.07.2011
 */

// local library
#include <scai/hmemo/cuda/CUDAContext.hpp>
#include <scai/hmemo/cuda/CUDAHostMemory.hpp>
#include <scai/hmemo/cuda/CUDAMemory.hpp>
#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>

#include <scai/hmemo/ContextAccess.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Assert.hpp>

// CUDA
#include <cublas_v2.h>
#include <cuda.h>
#include <cusparse.h>

// std
#include <memory>

namespace scai
{

cusparseHandle_t CUDAContext_cusparseHandle = 0;
cublasHandle_t CUDAContext_cublasHandle = 0;

using common::Thread;
using tasking::SyncToken;
using tasking::CUDAStreamSyncToken;

namespace hmemo
{

/**  static variables *****************************************************/

SCAI_LOG_DEF_LOGGER( CUDAContext::logger, "Context.CUDAContext" )

int CUDAContext::currentDeviceNr = -1;

int CUDAContext::numUsedDevices = 0;


/**  constructor  *********************************************************/

CUDAContext::CUDAContext( int deviceNr ) : 

    Context( common::context::CUDA ), 
    mDeviceNr( deviceNr )

{
    SCAI_LOG_DEBUG( logger, "construct CUDAContext, device nr = = " << deviceNr )

    // Note: logging is safe as CUDA context is always created after static initializations

    if ( numUsedDevices == 0 )
    {
        unsigned int flags = 0;    // must be set to zero
        SCAI_CUDA_DRV_CALL( cuInit( flags ), "cuInit failed, probably no GPU devices available" )
    }

    currentDeviceNr = deviceNr;
    numUsedDevices++;
    Context::enable( __FILE__, __LINE__ );
    SCAI_CUDA_DRV_CALL( cuDeviceGet( &mCUdevice, mDeviceNr ), "cuDeviceGet device " << mDeviceNr );

    cudaDeviceProp properties;
    cudaGetDeviceProperties( &properties, mCUdevice );

    // feature might be important to indicate that we can use CUDAHostMemory on device
    // SCAI_LOG_ERROR( logger, "canMapHostMemory = " << properties.canMapHostMemory );


    char deviceName[256];
    SCAI_CUDA_DRV_CALL( cuDeviceGetName( deviceName, 256, mCUdevice ), "cuDeviceGetName" );
    mDeviceName = deviceName; // save it as string member variable for output
    SCAI_LOG_DEBUG( logger, "got device " << mDeviceName )
    SCAI_CUDA_DRV_CALL( cuCtxCreate( &mCUcontext, CU_CTX_SCHED_SPIN | CU_CTX_MAP_HOST, mCUdevice ),
                        "cuCtxCreate for " << *this )
    SCAI_CUSPARSE_CALL( cusparseCreate( &CUDAContext_cusparseHandle ),
                        "Initialization of CUSparse library: cusparseCreate" );
    SCAI_CUBLAS_CALL( cublasCreate( &CUDAContext_cublasHandle ), "Initialization of CUBlas library: cublasCreate" );
    SCAI_LOG_INFO( logger, "Initialized: CUBLAS, CuSparse" )
    int flags = 0; // must be 0 by specification of CUDA driver API
    SCAI_CUDA_DRV_CALL( cuStreamCreate( &mTransferStream, flags ), "cuStreamCreate for transfer failed" )
    SCAI_CUDA_DRV_CALL( cuStreamCreate( &mComputeStream, flags ), "cuStreamCreate for compute failed" );
    mOwnerThread = Thread::getSelf(); // thread that can use the context
    disable( __FILE__, __LINE__ );
    // count used devices to shutdown CUDA if no more device is used
    SCAI_LOG_INFO( logger, *this << " constructed by thread " << mOwnerThread << ", now disabled" )
}

/**  destructor   *********************************************************/

CUDAContext::~CUDAContext()
{
    SCAI_LOG_INFO( logger, "~CUDAContext: " << *this )

    numUsedDevices--;

    // context must be valid before streams are destroyed

    CUresult res = cuCtxPushCurrent( mCUcontext );

    if ( res == CUDA_ERROR_DEINITIALIZED )
    {
        // this might happen with other software, e.g. VampirTrace
        SCAI_LOG_WARN( logger, "CUDA driver already deinitialized" )
        return;
    }
    else if ( res != CUDA_SUCCESS )
    {
        SCAI_LOG_ERROR( logger, "Could not push any more context for " << *this )
        return;
    }

    SCAI_LOG_DEBUG( logger, "pushed context: synchronize/destroy compute stream" )
    SCAI_CUDA_DRV_CALL( cuStreamSynchronize( mComputeStream ), "cuStreamSynchronize for compute failed" )
    SCAI_CUDA_DRV_CALL( cuStreamDestroy( mComputeStream ), "cuStreamDestroy for compute failed" );
    SCAI_LOG_DEBUG( logger, "synchronize/destroy transfer stream" )
    SCAI_CUDA_DRV_CALL( cuStreamSynchronize( mTransferStream ), "cuStreamSynchronize for transfer failed" )
    SCAI_CUDA_DRV_CALL( cuStreamDestroy( mTransferStream ), "cuStreamDestroy for transfer failed" );

    if ( !numUsedDevices )
    {
        SCAI_LOG_DEBUG( logger, "no more devices in use -> shutdown cuda" )

        if ( CUDAContext_cublasHandle )
        {
            cublasStatus_t error = cublasDestroy( CUDAContext_cublasHandle );

            if ( error != CUBLAS_STATUS_SUCCESS )
            {
                SCAI_LOG_ERROR( logger, "Could not destroy cublas handle, status = " << error )
            }

            SCAI_LOG_INFO( logger, "cublas handle successfully destroyed" )
            CUDAContext_cublasHandle = 0;
        }

        if ( CUDAContext_cusparseHandle )
        {
            cusparseStatus_t error = cusparseDestroy( CUDAContext_cusparseHandle );

            if ( error != CUSPARSE_STATUS_SUCCESS )
            {
                SCAI_LOG_ERROR( logger, "Could not destroy cusparse handle, status = " << error )
            }

            SCAI_LOG_INFO( logger, "cusparse handle successfully destroyed" )
            CUDAContext_cusparseHandle = 0;
        }
    }

    SCAI_LOG_DEBUG( logger, "destroy cuda context" )
    // do not do it if CUDA tracing is enabled
    SCAI_CUDA_DRV_CALL( cuCtxDestroy( mCUcontext ), "cuCtxDestroy failed" )

    // if we are current device, set it back to -1

    if ( currentDeviceNr == mDeviceNr )
    {
        currentDeviceNr = -1;
    }
}

/* ----------------------------------------------------------------------------- */

MemoryPtr CUDAContext::getMemoryPtr() const
{
    MemoryPtr memory;

    if ( mMemory.expired() )
    {
        memory.reset( new CUDAMemory( shared_from_this() ) );
        mMemory = memory; // save it here as a weak pointer to avoid cycles
    }
    else
    {
        // the last memory instance is still valid, so we return just shared pointer to it
        memory = mMemory.lock();
    }

    return memory;
}

/* ----------------------------------------------------------------------------- */

MemoryPtr CUDAContext::getHostMemoryPtr() const
{
    MemoryPtr memory;

    if ( mHostMemory.expired() )
    {
        memory.reset( new CUDAHostMemory( shared_from_this() ) );
        mHostMemory = memory; // save it here as a weak pointer to avoid cycles
    }
    else
    {
        // the last host context instance is still valid, so we return just shared pointer to it
        memory = mHostMemory.lock();
    }

    return memory;
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
    SCAI_LOG_DEBUG( logger, *this << ": disable by thread = " << Thread::getSelf() )
    CUcontext tmp; // temporary for last context, not necessary to save it
    SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" )
    currentDeviceNr = -1;
}

/* ----------------------------------------------------------------------------- */

void CUDAContext::enable( const char* file, int line ) const
{
    SCAI_LOG_DEBUG( logger, *this << ": enable by thread = " << Thread::getSelf() )
    Context::enable( file, line ); // call routine of base class
    SCAI_CUDA_DRV_CALL( cuCtxPushCurrent( mCUcontext ), "could not push context for " << *this )
//    lama_init0_cuda ();  // otherwise some properties of device might be wrong
    currentDeviceNr = mDeviceNr;
}

/* ----------------------------------------------------------------------------- */

bool CUDAContext::canUseMemory( const Memory& other ) const
{
    bool canUse = false;

    // CUDA device can use only data on same CUDA device

    if ( other.getType() == memtype::CUDAMemory )
    {
        const CUDAMemory* otherCUDAMem = dynamic_cast<const CUDAMemory*>( &other );

        SCAI_ASSERT( otherCUDAMem, "serious type mismatch" )

        canUse = otherCUDAMem->getDeviceNr() == mDeviceNr;
    }

    // Zero-Copy: we can use CUDA Host memory 

    if ( other.getType() == memtype::CUDAHostMemory )
    {
        const CUDAHostMemory* otherCUDAHostMem = dynamic_cast<const CUDAHostMemory*>( &other );

        SCAI_ASSERT( otherCUDAHostMem, "serious type mismatch" )

        canUse = otherCUDAHostMem->getCUDAContext().getDeviceNr() == mDeviceNr;
    }

    SCAI_LOG_DEBUG( logger, *this << ": " << ( canUse ? "can use " : "can't use " )
                            << other )

    return canUse;
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

#define SCAI_DEFAULT_DEVICE_NUMBER -1
#define SCAI_MAX_CUDA_DEVICES 4

static int getDefaultDeviceNr() 
{
    return 0;
}

/* ----------------------------------------------------------------------------- */
/*      Factory::Register - create( int )                                        */
/* ----------------------------------------------------------------------------- */

static common::weak_ptr<CUDAContext> mCUDAContext[SCAI_MAX_CUDA_DEVICES];

ContextPtr CUDAContext::create( int deviceNr )
{
    int cudaDeviceNr = deviceNr;

    if( cudaDeviceNr == SCAI_DEFAULT_DEVICE_NUMBER )
    {
        cudaDeviceNr = getDefaultDeviceNr();

        // no need here to check for a good value
    }
    else
    {
        SCAI_ASSERT(
            0 <= cudaDeviceNr && cudaDeviceNr < SCAI_MAX_CUDA_DEVICES,
            "device = " << cudaDeviceNr << " out of range" << ", max supported device = " << SCAI_MAX_CUDA_DEVICES )
    }

    common::shared_ptr<CUDAContext> context = common::shared_ptr<CUDAContext>();

    if( mCUDAContext[cudaDeviceNr].expired() )
    {
        // create a new context for the device and return the shared pointer

        context = common::shared_ptr<CUDAContext>( new CUDAContext( cudaDeviceNr ) );

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

} /* end namespace hmemo */

} /* end namespace scai */
