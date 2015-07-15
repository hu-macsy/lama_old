/**
 * @file CUDAHostMemory.cpp
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
 * @brief Implementation of methods for page-locked memory management to enable 
 *        fast DMA transfers to CUDA devices.
 *
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#include <cudamem/CUDAHostMemory.hpp>
#include <cudamem/CUDAError.hpp>
#include <cudamem/CUDAStreamSyncToken.hpp>

#include <memory/ContextAccess.hpp>

#include <common/Exception.hpp>

#include <common/bind.hpp>

// import ::memcpy 
#include <cstring>

using tasking::SyncToken;

namespace memory
{

LAMA_LOG_DEF_LOGGER( CUDAHostMemory::logger, "Memory.CUDAHostMemory" );

CUDAHostMemory::CUDAHostMemory( common::shared_ptr<const CUDAContext> cudaContext ) : 

    Memory( memtype::CUDAHostMemory ),
    mCUDAContext( cudaContext )

{
    COMMON_ASSERT( cudaContext, "CUDAHostMemory requires valid CUDAContext, is NULL" )
    LAMA_LOG_INFO( logger, "CUDAHostMemory created, allows faster transfer HOST <-> " << *mCUDAContext )
}

CUDAHostMemory::~CUDAHostMemory()
{
    LAMA_LOG_INFO( logger, "~CUDAHostMemory for " << *mCUDAContext )
}

void CUDAHostMemory::writeAt( std::ostream& stream ) const
{
    stream << "CUDAHostMemory( <-> " << *mCUDAContext << " )";
}

void* CUDAHostMemory::allocate( const size_t size ) const
{
    // LAMA_REGION( "CUDAHostMemory::allocate" )
    LAMA_LOG_TRACE( logger, *this << ": allocate " << size << " bytes" )

    void* pointer = 0;

    LAMA_CONTEXT_ACCESS( mCUDAContext );

    LAMA_CUDA_DRV_CALL( cuMemAllocHost( &pointer, size ), "cuMemAllocHost( size = " << size << " ) failed" )

    LAMA_LOG_DEBUG( logger, *this << ": allocated " << size << " bytes, pointer = " << pointer )

    return pointer;
}

void CUDAHostMemory::free( void* pointer, const size_t size ) const
{
    // LAMA_REGION( "CUDAHostMemory::free" )
    // Be careful: do not use
    // ContextAccess useCUDA( ContextPtr( mCUDAContext ) );
    // as this defines a function and not a variable
    // General rule: never use shared_ptr temporaries implicitly

    LAMA_CONTEXT_ACCESS( mCUDAContext )

    LAMA_CUDA_DRV_CALL( cuMemFreeHost( pointer ), "cuMemFreeHost( " << pointer << ", " << size << " ) failed" )

    LAMA_LOG_DEBUG( logger, *this << ": freed " << size << " bytes, pointer = " << pointer )
}

void CUDAHostMemory::memcpy( void* dst, const void* src, const size_t size ) const
{
    ::memcpy( dst, src, size );
}

SyncToken* CUDAHostMemory::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    LAMA_CONTEXT_ACCESS( mCUDAContext )

    std::auto_ptr<CUDAStreamSyncToken> syncToken( mCUDAContext->getTransferSyncToken() );

    LAMA_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (host) to " << dst << " (host) " )

    LAMA_CUDA_RT_CALL(
        cudaMemcpyAsync( dst, src, size, cudaMemcpyHostToHost, syncToken->getCUDAStream() ),
        "cudaMemcpyAsync( " << dst << ", " << src << ", " << size << ", " 
        << cudaMemcpyHostToHost << ", " << syncToken->getCUDAStream() << ") failed " )

    CUevent event;

    LAMA_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ), "Could not create event " )

    LAMA_CUDA_DRV_CALL( cuEventRecord( event, syncToken->getCUDAStream() ),
                        "cuEventRecord failed for CUevent " << event << '.' )

    syncToken->setEvent( event );

    return syncToken.release();
}

ContextPtr CUDAHostMemory::getContext() const
{
    // Currently Host device should do operations on Host memory
    // Possible extension: the corresponding CUDA device can also access the host memory
    //                     with limited PCIe bandwidth (Zero Copy, e.g. on Tegra K1)

    ContextPtr host = Context::getContext( context::Host );
    return host;
}

} // namespace

