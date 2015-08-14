/**
 * @file CUDAHostMemory.hpp
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
 * @brief Host memory management for page-locked memory for fast transfer to CUDA device
 * @author Thomas Brandes, Jiri Kraus
 * @date 04.07.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/memory/Memory.hpp>
#include <scai/tasking/SyncToken.hpp>

#include <scai/memory/cuda/CUDAContext.hpp>

#include <scai/logging.hpp>

#include <scai/common/weak_ptr.hpp>

namespace memory
{

/** Alternative memory to HostMemory so that memory will be allocated
 *  on pinned memory that allows faster transfer to a certain CUDA device.
 *
 *  As the CUDA Host memory will only allow faster transfer to a certain device
 *  and requires CUDA initialization, we will store also a shared pointer to the
 *  corresponding CUDA context.
 *
 *  Note: copy routines between CUDA Host and CUDA device are already provided 
 *        by CUDA memory class and are not required here.
 */

class COMMON_DLL_IMPORTEXPORT CUDAHostMemory: 

    public Memory 
{

public:

    CUDAHostMemory( common::shared_ptr<const CUDAContext> cudaContext );

    virtual ~CUDAHostMemory();

    virtual void* allocate( const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual bool canCopyFrom( const Memory& other ) const;

    virtual bool canCopyTo( const Memory& other ) const;

    virtual void memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    virtual void memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    /** On CUDAHostMemory we work usually with the HostContext. */

    virtual ContextPtr getContextPtr() const;

    const CUDAContext& getCUDAContext() const { return *mCUDAContext; }

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    virtual void writeAt( std::ostream& stream ) const;

    common::shared_ptr<const CUDAContext> mCUDAContext;   // fast DMA transfer to this device
};

}
