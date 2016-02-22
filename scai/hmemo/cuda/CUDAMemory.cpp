/**
 * @file CUDAMemory.cpp
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
 * @brief Contains the implementation of methods for the class CUDAMemory.
 * @author Thomas Brandes
 * @date 15.07.2011
 */

// hpp
#include <scai/hmemo/cuda/CUDAMemory.hpp>

// local library
#include <scai/hmemo/cuda/CUDAHostMemory.hpp>
#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>

#include <scai/hmemo/ContextAccess.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/bind.hpp>

// CUDA
#include <cuda.h>

// std
#include <memory>

namespace scai
{

using tasking::SyncToken;
using tasking::CUDAStreamSyncToken;

namespace hmemo
{

/**  static variables *****************************************************/

SCAI_LOG_DEF_LOGGER( CUDAMemory::logger, "Memory.CUDAMemory" )

/**  constructor  *********************************************************/

CUDAMemory::CUDAMemory( common::shared_ptr<const CUDAContext> cudaContext )

    : Memory( memtype::CUDAMemory ), 
      mCUDAContext( cudaContext )

{
    SCAI_ASSERT( cudaContext, "NULL context for CUDA memory" )

    SCAI_LOG_DEBUG( logger, "construct CUDAMemory for context " << cudaContext )

    mNumberOfAllocatedBytes = 0;
    mNumberOfAllocates = 0;
    mMaxNumberOfAllocatedBytes = 0;
}

/**  destructor   *********************************************************/

CUDAMemory::~CUDAMemory()
{
    SCAI_LOG_INFO( logger, "~CUDAMemory: " << *this )

    if ( mNumberOfAllocates > 0 )
    {
        SCAI_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " allocate without free" )
    }
    else if ( mNumberOfAllocates < 0 )
    {
        SCAI_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " free without allocate" )
    }

    if ( mNumberOfAllocatedBytes != 0 )
    {
        SCAI_LOG_ERROR( logger, *this << ": number of allocated bytes = " << mNumberOfAllocatedBytes 
                                << ", mismatch of free/allocate sizes" )
    }
}

/* ----------------------------------------------------------------------------- */

int CUDAMemory::getDeviceNr() const
{
    return mCUDAContext->getDeviceNr();
}

/* ----------------------------------------------------------------------------- */

ContextPtr CUDAMemory::getContextPtr() const
{
    return mCUDAContext;
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::writeAt( std::ostream& stream ) const
{
    stream << "CUDAMemory( @ " << *mCUDAContext << " )";
}

/* ----------------------------------------------------------------------------- */

void* CUDAMemory::allocate( const size_t size ) const
{
    // SCAI_REGION( "CUDA.allocate" )
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    SCAI_ASSERT( size > 0, "should not call allocate for size = " << size )
    SCAI_LOG_TRACE( logger, *this << ": allocate " << size << " bytes" )
    CUdeviceptr pointer = 0;
    SCAI_CUDA_DRV_CALL(
        cuMemAlloc( &pointer, size ),
        "cuMemAlloc( size = " << size << " ) failed. This allocation would require a total of " << mMaxNumberOfAllocatedBytes + size << " bytes global memory." )
    SCAI_LOG_DEBUG( logger, *this << ": allocated " << size << " bytes, ptr = " << ( ( void* ) pointer ) )
    mNumberOfAllocatedBytes += size;
    mNumberOfAllocates++;
    mMaxNumberOfAllocatedBytes = std::max( mNumberOfAllocatedBytes, mMaxNumberOfAllocatedBytes );
    return ( void* ) pointer;
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::free( void* pointer, const size_t size ) const
{
    // SCAI_REGION( "CUDA.free" )
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    SCAI_LOG_DEBUG( logger, *this << ": free " << size << " bytes, ptr = " << pointer )
    SCAI_CUDA_DRV_CALL( cuMemFree( ( CUdeviceptr ) pointer ), "cuMemFree( " << pointer << " ) failed" )
    mNumberOfAllocatedBytes -= size;
    mNumberOfAllocates--;
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpy( void* dst, const void* src, const size_t size ) const
{
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    SCAI_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (device) " )
    SCAI_CUDA_DRV_CALL( cuMemcpyDtoD( ( CUdeviceptr ) dst, ( CUdeviceptr ) src, size ),
                        "cuMemcpyDtoD( " << dst << ", " << src << ", " << size << " ) failed" )
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpyToCUDA( const CUDAMemory& dstMemory, void* dst, const void* src, const size_t size ) const
{
    SCAI_CONTEXT_ACCESS( mCUDAContext )

    unsigned int flags = 0;  // not any meaning now

    CUcontext dstCUcontext = dstMemory.mCUDAContext->getCUcontext();

    SCAI_CUDA_DRV_CALL( cuCtxEnablePeerAccess( dstCUcontext, flags ), "cuCtxEnablePeerAccess" )

    // unified adressing makes this possible

    SCAI_LOG_INFO( logger, "copy " << size << " bytes to " << dst << " @ " << dstMemory
                                   << " from " << src << " @ " << *this )

    SCAI_CUDA_DRV_CALL( cuMemcpyDtoD( ( CUdeviceptr ) dst, ( CUdeviceptr ) src, size ),
                        "cuMemcpyDtoD( " << dst << ", " << src << ", " << size << " ) failed" )

    SCAI_CUDA_DRV_CALL( cuCtxDisablePeerAccess( dstCUcontext ), "cuCtxDisablePeerAccess" )
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpyFromCUDA( void* dst, const CUDAMemory& srcMemory, const void* src, const size_t size ) const
{
    SCAI_CONTEXT_ACCESS( mCUDAContext )

    unsigned int flags = 0;  // not any meaning now

    CUcontext srcCUcontext = srcMemory.mCUDAContext->getCUcontext();

    SCAI_CUDA_DRV_CALL( cuCtxEnablePeerAccess( srcCUcontext, flags ), "cuCtxEnablePeerAccess" )

    // unified adressing makes this possible

    SCAI_LOG_INFO( logger, "copy " << size << " bytes from " << src << " @ " << srcMemory  
                                   << " to " << dst << " @ " << *this )

    SCAI_CUDA_DRV_CALL( cuMemcpyDtoD( ( CUdeviceptr ) dst, ( CUdeviceptr ) src, size ),
                        "cuMemcpyDtoD( " << dst << ", " << src << ", " << size << " ) failed" )

    SCAI_CUDA_DRV_CALL( cuCtxDisablePeerAccess( srcCUcontext ), "cuCtxDisablePeerAccess" )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAMemory::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    // SCAI_REGION( "CUDA.memcpyDtoDAsync" )
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    // use auto pointer so memory will be freed in case of exceptions
    SCAI_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (device) to " << dst << " (device) " )

    CUDAStreamSyncToken* token = mCUDAContext->getTransferSyncToken();

    SCAI_CUDA_DRV_CALL( cuMemcpyDtoDAsync( ( CUdeviceptr ) dst, ( CUdeviceptr ) src, size, token->getCUDAStream() ),
                        "cuMemcpyDtoDAsync( " << dst << ", " << src << ", " << size << ") failed " )

    // sync token should not synchronize on the full stream but only on the transfer, so add event

    CUevent event;

    SCAI_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ), "Could not create event " )

    SCAI_CUDA_DRV_CALL( cuEventRecord( event, token->getCUDAStream() ), "cuEventRecord failed for CUevent " << event << '.' )

    token->setEvent( event );

    return token;
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpyFromHost( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (host) to " << dst << " (device) " )

    SCAI_CONTEXT_ACCESS( mCUDAContext )

    SCAI_CUDA_DRV_CALL( cuMemcpyHtoD( ( CUdeviceptr ) dst, src, size ),
                        "cuMemcpyHToD( " << dst << ", " << src << ", " << size << ") failed " )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAMemory::memcpyAsyncFromHost( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "async copy " << size << " bytes from " << src << " (host) to " << dst << " (device) " )

    const size_t THRESHOLD_SIZE = 16 * 1024;   // number of bytes where new thread might be useful

    if ( size > THRESHOLD_SIZE )
    {
        return new tasking::TaskSyncToken( common::bind( &CUDAMemory::memcpyFromHost, this, dst, src, size ) );
    }
    else
    {
        memcpyFromHost( dst, src, size );
        return NULL;
    }
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpyToHost( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (host) " )

    SCAI_CONTEXT_ACCESS( mCUDAContext )

    SCAI_CUDA_DRV_CALL( cuMemcpyDtoH( dst, ( CUdeviceptr ) src, size ),
                        "cuMemcpyDToH( " << dst << ", " << src << ", " << size << ") failed " )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAMemory::memcpyAsyncToHost( void* dst, const void* src, const size_t size ) const
{
    SCAI_LOG_INFO( logger, "async copy " << size << " bytes from " << src << " (device) to " << dst << " (host) " )

    const size_t THRESHOLD_SIZE = 16 * 1024;   // number of bytes where new thread might be useful

    if ( size > THRESHOLD_SIZE )
    {
        return new tasking::TaskSyncToken( common::bind( &CUDAMemory::memcpyToHost, this, dst, src, size ) );
    }
    else
    {
        memcpyToHost( dst, src, size );
        return NULL;
    }
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpyFromCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // SCAI_REGION( "CUDA.memcpyCUDAHost->Dev")
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    SCAI_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (host) to " << dst << " (device) " )
    SCAI_CUDA_DRV_CALL( cuMemcpyHtoD( ( CUdeviceptr ) dst, src, size ),
                        "cuMemcpyHToD( " << dst << ", " << src << ", " << size << ") failed " )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAMemory::memcpyAsyncFromCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // SCAI_REGION( "CUDA.memcpyHtoDAsync" )
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    SCAI_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (host) to " << dst << " (device) " )

    CUDAStreamSyncToken* token = mCUDAContext->getTransferSyncToken();

    SCAI_CUDA_DRV_CALL( cuMemcpyHtoDAsync( ( CUdeviceptr ) dst, src, size, token->getCUDAStream() ),
                        "cuMemcpyHtoDAsync( " << dst << ", " << src << ", " << size << ") failed " )

    return token;
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpyToCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // SCAI_REGION( "CUDA.memcpyDev->CUDAHost" )
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    SCAI_LOG_INFO( logger, "copy " << size << " bytes from " << src << " (device) to " << dst << " (cuda host) " )
    SCAI_CUDA_DRV_CALL( cuMemcpyDtoH( dst, ( CUdeviceptr ) src, size ),
                        "cuMemcpyDToH( " << dst << ", " << src << ", " << size << ") failed " )
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAMemory::memcpyAsyncToCUDAHost( void* dst, const void* src, const size_t size ) const
{
    // SCAI_REGION( "CUDA.memcpyDtoHAsync" )
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    SCAI_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (device) to " << dst << " (host) " )

    CUDAStreamSyncToken* token = mCUDAContext->getTransferSyncToken();

    SCAI_CUDA_DRV_CALL(
        cuMemcpyDtoHAsync( dst, ( CUdeviceptr ) src, size, token->getCUDAStream() ),
        "cuMemcpyDtoHAsync( " << dst << ", " << src << ", " << size << ", " << token->getCUDAStream() << ") failed " )
    // sync token should not synchronize on the full stream but only on the transfer, so add event
    CUevent event;
    SCAI_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ), "Could not create event " )
    SCAI_CUDA_DRV_CALL( cuEventRecord( event, token->getCUDAStream() ), "cuEventRecord failed for CUevent " << event << '.' )

    token->setEvent( event );

    return token;
}

/* ----------------------------------------------------------------------------- */

bool CUDAMemory::canCopyFrom( const Memory& other ) const
{
    bool supported = false;

    memtype::MemoryType otherType = other.getType();

    if ( otherType == memtype::HostMemory )
    {
        // CUDADevice -> Host is supported

        supported = true;
    }
    else if ( otherType == memtype::CUDAHostMemory )
    {
        // CUDADevice -> CUDA Host is supported
        // Note: slower but okay if CUDA Host memory does not belong to this device 

        supported = true;
    }
    else if ( otherType == memtype::CUDAMemory )
    {
        const CUDAMemory* otherCUDAMem = dynamic_cast<const CUDAMemory*>( &other );

        SCAI_ASSERT( otherCUDAMem, "dynamic_cast<CUDAMemory*> failed" )

        supported = canCopyCUDA( *otherCUDAMem );
    }

    SCAI_LOG_DEBUG( logger, "canCopyFrom " << other << " to this " << *this << ", supported = " << supported )

    return supported;
}

/* ----------------------------------------------------------------------------- */

bool CUDAMemory::canCopyCUDA( const CUDAMemory& other ) const
{
    bool supported = false;

    if ( other.getDeviceNr() == getDeviceNr() )
    {
        supported = true;
    }
    else
    {
        // Check for the access capability

        SCAI_CONTEXT_ACCESS( mCUDAContext )

	    int accessCapability = 0;

        SCAI_CUDA_DRV_CALL(
            cuDeviceCanAccessPeer( &accessCapability, getDeviceNr(), other.getDeviceNr() ),
            "cuDeviceCanAccessPeer failed" );

        if ( accessCapability >= 1 )
        {
            supported = true;
        }
    }

    return supported;
}

/* ----------------------------------------------------------------------------- */

bool CUDAMemory::canCopyTo( const Memory& other ) const
{
    bool supported = false;

    memtype::MemoryType otherType = other.getType();

    if ( otherType == memtype::HostMemory )
    {
        // CUDAMemory -> HostMemory is supported

        supported = true;
    }
    else if ( otherType == memtype::CUDAHostMemory )
    {
        // CUDAMemory -> CUDA Host is supported

        supported = true;
    }
    else if ( otherType == memtype::CUDAMemory )
    {
        const CUDAMemory* otherCUDA = dynamic_cast<const CUDAMemory*>( &other );

        SCAI_ASSERT( otherCUDA, "dynamic_cast<CUDAMemory*> failed" )

        supported = canCopyCUDA( *otherCUDA );
    }

    SCAI_LOG_DEBUG( logger, "canCopyTo " << other << " from this " << *this << ", supported = " << supported )

    return supported;
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const
{
    if ( srcMemory.getType() == memtype::HostMemory )
    {
        memcpyFromHost( dst, src, size );
    }
    else if ( srcMemory.getType() == memtype::CUDAHostMemory )
    {
        memcpyFromCUDAHost( dst, src, size );
    }
    else if ( srcMemory.getType() == memtype::CUDAMemory )
    {
        const CUDAMemory* srcCUDAMemory = dynamic_cast<const CUDAMemory*>( &srcMemory );

        SCAI_ASSERT( srcCUDAMemory, "dynamic_cast<CUDAMemory*> failed" )

        memcpyFromCUDA( dst, *srcCUDAMemory, src, size );
    }
    else
    {
        SCAI_LOG_ERROR( logger, "copy from " << srcMemory << " to " << *this << " not supported" )
        COMMON_THROWEXCEPTION( "copy from " << srcMemory << " to " << *this << " not supported" )
    }
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAMemory::memcpyFromAsync( void* dst, const Memory& srcMemory, const void* src, size_t size ) const
{
    if ( srcMemory.getType() == memtype::HostMemory )
    {
        // here we have no asynchronous transfer

        return memcpyAsyncFromHost( dst, src, size );
    }
    else if ( srcMemory.getType() == memtype::CUDAHostMemory )
    {
        return memcpyAsyncFromCUDAHost( dst, src, size );
    }
    else if ( srcMemory.getType() == memtype::CUDAMemory )
    {
        const CUDAMemory* srcCUDAMemory = dynamic_cast<const CUDAMemory*>( &srcMemory );

        SCAI_ASSERT( srcCUDAMemory, "dynamic_cast<CUDAMemory*> failed" )

        memcpyFromCUDA( dst, *srcCUDAMemory, src, size );
    }
    else
    {
        SCAI_LOG_ERROR( logger, "copy from " << srcMemory << " to " << *this << " not supported" )
        COMMON_THROWEXCEPTION( "copy from " << srcMemory << " to " << *this << " not supported" )
    }

    return NULL;
}

/* ----------------------------------------------------------------------------- */

void CUDAMemory::memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const
{
    if ( dstMemory.getType() == memtype::HostMemory )
    {
        memcpyToHost( dst, src, size );
    }
    else if ( dstMemory.getType() == memtype::CUDAHostMemory )
    {
        memcpyToCUDAHost( dst, src, size );
    }
    else if ( dstMemory.getType() == memtype::CUDAMemory )
    {
        const CUDAMemory* dstCUDAMemory = dynamic_cast<const CUDAMemory*>( &dstMemory );

        SCAI_ASSERT( dstCUDAMemory, "dynamic_cast<CUDAMemory*> failed" )

        memcpyToCUDA( *dstCUDAMemory, dst, src, size );
    }
    else
    {
        SCAI_LOG_ERROR( logger, "copy to " << dstMemory << " from " << *this << " not supported" )
        COMMON_THROWEXCEPTION( "copy to " << dstMemory << " from " << *this << " not supported" )
    }
}

/* ----------------------------------------------------------------------------- */

SyncToken* CUDAMemory::memcpyToAsync( const Memory& dstMemory, void* dst, const void* src, size_t size ) const
{
    if ( dstMemory.getType() == memtype::HostMemory )
    {
        return memcpyAsyncToHost( dst, src, size );
    }
    else if ( dstMemory.getType() == memtype::CUDAHostMemory )
    {
        return memcpyAsyncToCUDAHost( dst, src, size );
    }
    else if ( dstMemory.getType() == memtype::CUDAMemory )
    {
        const CUDAMemory* dstCUDAMemory = dynamic_cast<const CUDAMemory*>( &dstMemory );

        SCAI_ASSERT( dstCUDAMemory, "dynamic_cast<CUDAMemory*> failed" )

        memcpyToCUDA( *dstCUDAMemory, dst, src, size );
    }
    else
    {
        SCAI_LOG_ERROR( logger, "copy to " << dstMemory << " from " << *this << " not supported" )
        COMMON_THROWEXCEPTION( "copy to " << dstMemory << " from " << *this << " not supported" )
    }

    return NULL;
}

/* ----------------------------------------------------------------------------- */

} /* end namespace hmemo */

} /* end namespace scai */
