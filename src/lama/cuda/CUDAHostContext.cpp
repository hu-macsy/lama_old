/**
 * @file CUDAHostContext.cpp
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
 * @brief Definition of a host context that allocates pinned memory so that
 *        memory transfer for CUDA context is faster.
 * @author Thomas Brandes
 * @date 16.07.2011
 * @since 1.0.0
 */

// hpp
#include <lama/cuda/CUDAHostContext.hpp>

// others
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>

#include <lama/ContextAccess.hpp>

#include <lama/exception/LAMAAssert.hpp>

// tracing
#include <lama/tracing.hpp>

#include <boost/bind.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDAHostContext::logger, "Context.CUDAHostContext" );

CUDAHostContext::CUDAHostContext( boost::shared_ptr<const CUDAContext> cudaContext )
    : mCUDAContext( cudaContext )
{
    if( !cudaContext )
    {
        LAMA_THROWEXCEPTION( "CUDAHostContext requires valid CUDAContext, is NULL" )
    }

    LAMA_LOG_INFO( logger, "CUDAHostContext created, allows faster transfer HOST <-> " << *mCUDAContext )
}

CUDAHostContext::~CUDAHostContext()
{
}

void CUDAHostContext::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "CUDAHostContext(" << mCUDAContext->getDeviceNr() << ":lama_host_Alloc/Free_cuda)";
}

void* CUDAHostContext::allocate( const size_t size ) const
{
    LAMA_REGION( "CUDAHostContext::allocate" )
    LAMA_LOG_TRACE( logger, *this << ": allocate " << size << " bytes" )

    void* pointer = 0;

    LAMA_CONTEXT_ACCESS( mCUDAContext );

    LAMA_CUDA_DRV_CALL( cuMemAllocHost( &pointer, size ), "cuMemAllocHost( size = " << size << " ) failed" )

    LAMA_LOG_DEBUG( logger, *this << ": allocated " << size << " bytes, pointer = " << pointer )

    return pointer;
}

void CUDAHostContext::allocate( ContextData& contextData, const size_t size ) const
{
    contextData.pointer = allocate( size );
    contextData.setPinned();
}

void CUDAHostContext::free( void* pointer, const size_t size ) const
{
    LAMA_REGION( "CUDAHostContext::free" )
    // Be careful: do not use
    // ContextAccess useCUDA( ContextPtr( mCUDAContext ) );
    // as this defines a function and not a variable
    // General rule: never use shared_ptr temporaries implicitly

    LAMA_CONTEXT_ACCESS( mCUDAContext )

    LAMA_CUDA_DRV_CALL( cuMemFreeHost( pointer ), "cuMemFreeHost( " << pointer << ", " << size << " ) failed" )

    LAMA_LOG_DEBUG( logger, *this << ": freed " << size << " bytes, pointer = " << pointer )
}

void CUDAHostContext::free( ContextData& contextData ) const
{
    LAMA_ASSERT_EQUAL_ERROR( contextData.context->getType(), getType() )
    free( contextData.pointer, contextData.size );
}

void CUDAHostContext::memcpy( void* dst, const void* src, const size_t size ) const
{
    ::memcpy( dst, src, size );
}

SyncToken* CUDAHostContext::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    LAMA_CONTEXT_ACCESS( mCUDAContext )

    std::auto_ptr<CUDAStreamSyncToken> syncToken( mCUDAContext->getTransferSyncToken() );

    LAMA_LOG_INFO( logger, "copy async " << size << " bytes from " << src << " (host) to " << dst << " (host) " )

    LAMA_CUDA_RT_CALL(
        cudaMemcpyAsync( dst, src, size, cudaMemcpyHostToHost, syncToken->getCUDAStream() ),
        "cudaMemcpyAsync( " << dst << ", " << src << ", " << size << ", " << cudaMemcpyHostToHost << ", " << syncToken->getCUDAStream() << ") failed " )

    CUevent event;

    LAMA_CUDA_DRV_CALL( cuEventCreate( &event, CU_EVENT_DEFAULT | CU_EVENT_DISABLE_TIMING ), "Could not create event " )

    LAMA_CUDA_DRV_CALL( cuEventRecord( event, syncToken->getCUDAStream() ),
                        "cuEventRecord failed for CUevent " << event << '.' )

    syncToken->setEvent( event );

    return syncToken.release();
}

bool CUDAHostContext::cancpy( const ContextData& dst, const ContextData& src ) const
{
    return dst.context->getType() == getType() && src.context->getType() == getType();
}

void CUDAHostContext::memcpy( ContextData& dst, const ContextData& src, const size_t size ) const
{
    LAMA_ASSERT_ERROR( dst.context->getType() == getType() && src.context->getType() == getType(),
                       "Can not copy from "<< *(src.context) << " to " << *(dst.context) )
    memcpy( dst.pointer, src.pointer, size );
}

SyncToken* CUDAHostContext::memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const
{
    LAMA_ASSERT_ERROR( dst.context->getType() == getType() && src.context->getType() == getType(),
                       "Can not copy from "<< *(src.context) << " to " << *(dst.context) )
    return memcpyAsync( dst.pointer, src.pointer, size );
}

HostContext::HostContextType CUDAHostContext::getHostType() const
{
    return CUDAHost;
}

} // namespace lama

