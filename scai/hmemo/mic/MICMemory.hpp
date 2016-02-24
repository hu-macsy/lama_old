/**
 * @file MICMemory.hpp
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
 * @brief Definition of context class for MIC devices and a context manager class.
 * @author Thomas Brandes
 * @date 01.07.2013
 * @since 1.1.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Memory.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/shared_ptr.hpp>

// std
#include <string>

namespace scai
{

namespace hmemo
{

/**
 * @brief MICContext initializes the MIC device with the given number.
 *
 */
class COMMON_DLL_IMPORTEXPORT MICMemory :

    public Memory
{

public:

    /**
     * @brief Constructor for the MIC memory management.
     */
    MICMemory( common::shared_ptr<const class MICContext> micContext );

    /**
     * @brief The destructor tests for all data freed.
     */
    virtual ~MICMemory();

    int getDeviceNr() const;

    virtual void writeAt( std::ostream& stream ) const;

    virtual void* allocate( const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual bool canCopyFrom( const Memory& other ) const;

    virtual bool canCopyTo( const Memory& other ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual void memset( void* dst, const int val, const size_t size ) const;

    virtual void memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    virtual void memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    virtual tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual ContextPtr getContextPtr() const;

private:

    virtual void memcpyToHost( void* dst, const void* src, const size_t size ) const;

    virtual void memcpyFromHost( void* dst, const void* src, const size_t size ) const;

    common::shared_ptr<const MICContext> mMICContext;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace hmemo */

} /* end namespace scai */
