/**
 * @file CUDAHostContext.hpp
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
 * @brief Host context where memory is pinned for fast transfer to CUDA device
 * @author Thomas Brandes, Jiri Kraus
 * @date 04.07.2011
 * @since 1.0.0
 */
#pragma once

// for dll_import
#include <common/config.hpp>

// others
#include <memory/HostContext.hpp>
#include <memory/SyncToken.hpp>

#include <cudamem/CUDAContext.hpp>

#include <logging/logging.hpp>

// boost
#include <boost/weak_ptr.hpp>

namespace memory
{

/** Alternative context to DefaultHostContext so that memory will be allocated
 *  on pinned memory that allows faster transfer to a certain CUDA device.
 *
 *  As the CUDA Host memory will only allow faster transfer to a certain device
 *  and requires CUDA initialization, we will store also a shared pointer to the
 *  corresponding CUDA context.
 */

class COMMON_DLL_IMPORTEXPORT CUDAHostContext: public HostContext
{
    friend class CUDAContext; // can only create this host context

public:

    virtual ~CUDAHostContext();

    virtual void* allocate( const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual HostContextType getHostType() const;

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    CUDAHostContext( boost::shared_ptr<const CUDAContext> cudaContext );

    virtual void writeAt( std::ostream& stream ) const;

    boost::shared_ptr<const CUDAContext> mCUDAContext;
};

}
