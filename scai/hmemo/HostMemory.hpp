/**
 * @file HostMemory.hpp
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
 * @brief Definion of memory class for usual Host/CPU memory.
 * @author Thomas Brandes
 * @date 14.07.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Memory.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/Thread.hpp>
#include <scai/common/shared_ptr.hpp>

namespace scai
{

namespace hmemo
{

/** @brief This class implements the default HOST memory.
 *
 *  This class is implemented as a singleton, only one default host
 *  memory is available.
 *
 *  The host memory allocates/frees data in the usual way.
 */

class COMMON_DLL_IMPORTEXPORT HostMemory: public Memory
{

public:

    HostMemory( common::shared_ptr<const class HostContext> hostContext );   

    virtual ~HostMemory();

    virtual void writeAt( std::ostream& stream ) const;

    virtual void* allocate( const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual void memset( void* dst, const int val, const size_t size ) const;

    /** This routine implements Context::memcpyAsync  */

    virtual tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual ContextPtr getContextPtr() const;

    /** This routine returns the singleton instance of the HostMemory. */

    static MemoryPtr getIt();

private:

    common::shared_ptr<const HostContext> mHostContextPtr;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    mutable size_t mNumberOfAllocates; //!< variable counts allocates

    mutable size_t mNumberOfAllocatedBytes;//!< variable counts allocated bytes

    mutable common::Thread::RecursiveMutex allocate_mutex;// needed to make allocate/free thread-safe
};

} /* end namespace hmemo */

} /* end namespace scai */
