/**
 * @file CUDAHostMemory.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Implementation of methods for page-locked memory management to enable
 *        fast DMA transfers to CUDA devices.
 * @author Thomas Brandes
 * @date 16.07.2015
 */

// hpp
#include <scai/hmemo/cuda/CUDAHostMemory.hpp>

// local library

#include <scai/hmemo/ContextAccess.hpp>

// internal scai libraries

#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/macros/assert.hpp>

// std
#include <cstring> // import ::memcpy
#include <memory>
#include <functional>

namespace scai
{

using tasking::SyncToken;
using tasking::CUDAStreamSyncToken;

namespace hmemo
{

SCAI_LOG_DEF_LOGGER( CUDAHostMemory::logger, "Memory.CUDAHostMemory" );

CUDAHostMemory::CUDAHostMemory( std::shared_ptr<const CUDAContext> cudaContext ) :

    Memory( MemoryType::CUDAHostMemory ),
    mCUDAContext( cudaContext ),
    mNumberOfAllocates( 0 ),
    mNumberOfAllocatedBytes( 0 ),
    mMaxAllocatedBytes( 0 )

{
    SCAI_ASSERT( cudaContext, "CUDAHostMemory requires valid CUDAContext, is NULL" )
    SCAI_LOG_INFO( logger, "CUDAHostMemory created, allows faster transfer HOST <-> " << *mCUDAContext )
}

CUDAHostMemory::~CUDAHostMemory()
{
    SCAI_LOG_INFO( logger, "~CUDAHostMemory for " << *mCUDAContext )
}

void CUDAHostMemory::writeAt( std::ostream& stream ) const
{
    stream << "CUDAHostMemory( <-> " << *mCUDAContext << " )";
}

void* CUDAHostMemory::allocate( const size_t size ) const
{
    SCAI_LOG_TRACE( logger, *this << ": allocate " << size << " bytes" )
    void* pointer = 0;
    SCAI_CONTEXT_ACCESS( mCUDAContext );
    SCAI_CUDA_DRV_CALL( cuMemAllocHost( &pointer, size ), "cuMemAllocHost( size = " << size << " ) failed" )
    SCAI_LOG_DEBUG( logger, *this << ": allocated " << size << " bytes, pointer = " << pointer )
    unsigned int flags = 0;
    void* pDevice = NULL;
    // check if we can use HostMemory also for device computations
    SCAI_CUDA_RT_CALL( cudaHostGetDevicePointer( &pDevice, pointer, flags ), "cudaHostGetDevicePointer" )
    SCAI_ASSERT_EQUAL( pDevice, pointer, "Not yet supported: pointer conversion for different context" )

    mNumberOfAllocatedBytes += size;
    mMaxAllocatedBytes = std::max( mMaxAllocatedBytes, mNumberOfAllocatedBytes );
    mNumberOfAllocates++;

    return pointer;
}

void CUDAHostMemory::free( void* pointer, const size_t size ) const
{
    // SCAI_REGION( "CUDAHostMemory::free" )
    // Be careful: do not use
    // ContextAccess useCUDA( ContextPtr( mCUDAContext ) );
    // as this defines a function and not a variable
    // General rule: never use shared_ptr temporaries implicitly
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    SCAI_CUDA_DRV_CALL_NOTHROW( cuMemFreeHost( pointer ), "cuMemFreeHost( " << pointer << ", " << size << " ) failed" )
    SCAI_LOG_DEBUG( logger, *this << ": freed " << size << " bytes, pointer = " << pointer )
    mNumberOfAllocatedBytes -= size;
    mNumberOfAllocates--;
}

void CUDAHostMemory::memcpy( void* dst, const void* src, const size_t size ) const
{
    ::memcpy( dst, src, size );
}

void CUDAHostMemory::memset( void* dst, const int val, const size_t size ) const
{
    ::memset( dst, val, size );
}

SyncToken* CUDAHostMemory::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    SCAI_CONTEXT_ACCESS( mCUDAContext )
    std::unique_ptr<CUDAStreamSyncToken> syncToken( mCUDAContext->getTransferSyncToken() );
    SCAI_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (host) to " << dst << " (host) " )
    SCAI_CUDA_RT_CALL(
        cudaMemcpyAsync( dst, src, size, cudaMemcpyHostToHost, syncToken->getCUDAStream() ),
        "cudaMemcpyAsync( " << dst << ", " << src << ", " << size << ", "
        << cudaMemcpyHostToHost << ", " << syncToken->getCUDAStream() << ") failed " )
    CUevent event;
    SCAI_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ), "Could not create event " )
    SCAI_CUDA_DRV_CALL( cuEventRecord( event, syncToken->getCUDAStream() ),
                        "cuEventRecord failed for CUevent " << event << '.' )
    syncToken->setEvent( event );
    return syncToken.release();
}

ContextPtr CUDAHostMemory::getContextPtr() const
{
    // Currently Host device should do operations on Host memory
    // Possible extension: the corresponding CUDA device can also access the host memory
    //                     with limited PCIe bandwidth (Zero Copy, e.g. on Tegra K1)
    ContextPtr host = Context::getContextPtr( common::ContextType::Host );
    return host;
}

/* ----------------------------------------------------------------------------- */

bool CUDAHostMemory::canCopyFrom( const Memory& other ) const
{
    bool supported = false;
    MemoryType otherType = other.getType();

    if ( otherType == MemoryType::HostMemory )
    {
        // CUDHostMem -> Host is supported
        supported = true;
    }
    else if ( otherType == MemoryType::CUDAHostMemory )
    {
        supported = true;
    }

    return supported;
}

bool CUDAHostMemory::canCopyTo( const Memory& other ) const
{
    bool supported = false;
    MemoryType otherType = other.getType();

    if ( otherType == MemoryType::HostMemory )
    {
        // CUDHostMem -> Host is supported
        supported = true;
    }
    else if ( otherType == MemoryType::CUDAHostMemory )
    {
        supported = true;
    }

    return supported;
}

/* ----------------------------------------------------------------------------- */

void CUDAHostMemory::memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const
{
    // all kind of Host <-> CUDAHost is supported
    if ( srcMemory.getType() == MemoryType::HostMemory )
    {
        ::memcpy( dst, src, size );
    }
    else if ( srcMemory.getType() == MemoryType::CUDAHostMemory )
    {
        ::memcpy( dst, src, size );
    }
    else
    {
        SCAI_LOG_ERROR( logger, "copy from " << srcMemory << " to " << *this << " not supported" )
        COMMON_THROWEXCEPTION( "copy from " << srcMemory << " to " << *this << " not supported" )
    }
}

/* ----------------------------------------------------------------------------- */

void CUDAHostMemory::memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const
{
    // all kind of Host <-> CUDAHost is supported
    if ( dstMemory.getType() == MemoryType::HostMemory )
    {
        ::memcpy( dst, src, size );
    }
    else if ( dstMemory.getType() == MemoryType::CUDAHostMemory )
    {
        ::memcpy( dst, src, size );
    }
    else
    {
        SCAI_LOG_ERROR( logger, "copy to " << dstMemory << " from " << *this << " not supported" )
        COMMON_THROWEXCEPTION( "copy to " << dstMemory << " from " << *this << " not supported" )
    }
}

/* ----------------------------------------------------------------------------- */

size_t CUDAHostMemory::maxAllocatedBytes() const
{
    return mMaxAllocatedBytes;
}

/* ----------------------------------------------------------------------------- */

} /* end namespace hmemo */

} /* end namespace scai */
