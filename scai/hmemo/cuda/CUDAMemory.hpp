/**
 * @file CUDAMemory.hpp
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
 * @brief Definition of memory class for CUDA devices.
 * @author Thomas Brandes, Jiri Kraus
 * @date 08.07.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Memory.hpp>

// internal scai libraries
#include <scai/common/Thread.hpp>

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
    CUDAMemory( common::shared_ptr<const class CUDAContext> cudaContext );

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

    virtual void memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    /** Overrides Memory::memcpyFromAsync */

    virtual tasking::SyncToken* memcpyFromAsync( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    virtual void memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    /** Overrides Memory::memcpyFromAsync */

    virtual tasking::SyncToken* memcpyToAsync( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    virtual tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual ContextPtr getContextPtr() const;

private:

    bool canCopyCUDA( const CUDAMemory& other ) const;

    common::shared_ptr<const CUDAContext> mCUDAContext;

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

    mutable int mNumberOfAllocates; //!< variable counts allocates
    mutable long long mNumberOfAllocatedBytes; //!< variable counts allocated bytes on device
    mutable long long mMaxNumberOfAllocatedBytes; //!< variable counts the maximum allocated bytes

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace hmemo */

} /* end namespace scai */
