/**
 * @file CUDAHostMemory.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Host memory management for page-locked memory for fast transfer to CUDA device
 * @author Thomas Brandes, Jiri Kraus
 * @date 04.07.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Memory.hpp>

// local library
#include <scai/hmemo/cuda/CUDAContext.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <memory>

namespace scai
{

namespace hmemo
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

    CUDAHostMemory( std::shared_ptr<const CUDAContext> cudaContext );

    virtual ~CUDAHostMemory();

    virtual void* allocate( const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual void memset( void* dst, const int val, const size_t size ) const;

    virtual tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual bool canCopyFrom( const Memory& other ) const;

    virtual bool canCopyTo( const Memory& other ) const;

    virtual void memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    virtual void memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    /** On CUDAHostMemory we work usually with the HostContext. */

    virtual ContextPtr getContextPtr() const;

    const CUDAContext& getCUDAContext() const
    {
        return *mCUDAContext;
    }

    size_t maxAllocatedBytes() const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    virtual void writeAt( std::ostream& stream ) const;

    std::shared_ptr<const CUDAContext> mCUDAContext;   // fast DMA transfer to this device

private:

    mutable size_t mNumberOfAllocates;
    mutable size_t mNumberOfAllocatedBytes;
    mutable size_t mMaxAllocatedBytes;
};

} /* end namespace hmemo */

} /* end namespace scai */
