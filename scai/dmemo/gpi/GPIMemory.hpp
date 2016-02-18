/**
 * @file GPIMemory.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Host context where memory is pinned for fast transfer via remote read/write
 * @author Thomas Brandes
 * @date 16.05.2014
 * @since 1.0.1
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/hmemo/Memory.hpp>
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

namespace scai
{

namespace dmemo
{

/** 
 *  GPIMemory is just pinned memory that can be used on the host like usual memory.
 */

class COMMON_DLL_IMPORTEXPORT GPIMemory: public hmemo::Memory
{
public:

    GPIMemory();

    virtual ~GPIMemory();

    virtual void* allocate( const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual bool canCopyFrom( const Memory& other ) const;

    virtual bool canCopyTo( const Memory& other ) const;

    virtual void memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    virtual void memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    virtual hmemo::ContextPtr getContextPtr() const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    virtual void writeAt( std::ostream& stream ) const;
};

}

}
