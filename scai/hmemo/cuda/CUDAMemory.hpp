/**
 * @file CUDAMemory.hpp
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
 * @brief Definition of memory class for CUDA devices.
 * @author Thomas Brandes, Jiri Kraus
 * @date 08.07.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Memory.hpp>

// CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

// std
#include <string>

namespace scai
{

namespace tasking
{
class COMMON_DLL_IMPORTEXPORT CUDAStreamSyncToken;
}

namespace hmemo
{

/**
 * @brief CUDAMemory initializes the CUDA device with the given number.
 *
 * CUDAMemory keeps always a shared pointer to the CUDA context as most
 * operations require also access to the context.
 */
class COMMON_DLL_IMPORTEXPORT CUDAMemory:

    public Memory
{

public:

    /**
     * @brief Constructor for the CUDA memory management.
     */
    CUDAMemory( std::shared_ptr<const class CUDAContext> cudaContext );

    /**
     * @brief The destructor destroys this CUDA device, and frees the initialized
     *        CUDA device if needed.
     */
    virtual ~CUDAMemory();

    int getDeviceNr() const;

    virtual bool canCopyFrom( const Memory& other ) const;

    virtual bool canCopyTo( const Memory& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    virtual void* allocate( const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual void memset( void* dst, const int val, const size_t size ) const;

    virtual void memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    /** Overrides Memory::memcpyFromAsync */

    virtual tasking::SyncToken* memcpyFromAsync( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    virtual void memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    /** Overrides Memory::memcpyFromAsync */

    virtual tasking::SyncToken* memcpyToAsync( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    virtual tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual ContextPtr getContextPtr() const;

    /** Implementation of Memory::maxAllocatedbytes */

    virtual size_t maxAllocatedBytes() const;

private:

    bool canCopyCUDA( const CUDAMemory& other ) const;

    std::shared_ptr<const CUDAContext> mCUDAContext;

    void memcpyFromHost( void* dst, const void* src, const size_t size ) const;
    void memcpyToHost( void* dst, const void* src, const size_t size ) const;
    void memcpyFromCUDAHost( void* dst, const void* src, const size_t size ) const;
    void memcpyToCUDAHost( void* dst, const void* src, const size_t size ) const;

    void memcpyToCUDA( const CUDAMemory& dstMemory, void* dst, const void* src, const size_t size ) const;
    void memcpyFromCUDA( void* dst, const CUDAMemory& srcMemory, const void* src, const size_t size ) const;

    tasking::SyncToken* memcpyAsyncFromHost( void* dst, const void* src, const size_t size ) const;
    tasking::SyncToken* memcpyAsyncToHost( void* dst, const void* src, const size_t size ) const;
    tasking::SyncToken* memcpyAsyncFromCUDAHost( void* dst, const void* src, const size_t size ) const;
    tasking::SyncToken* memcpyAsyncToCUDAHost( void* dst, const void* src, const size_t size ) const;

    mutable size_t mNumberOfAllocates; //!< variable counts allocates
    mutable size_t mNumberOfAllocatedBytes; //!< variable counts allocated bytes on device
    mutable size_t mMaxAllocatedBytes; //!< variable counts the maximum allocated bytes

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace hmemo */

} /* end namespace scai */
