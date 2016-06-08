/**
 * @file MICMemory.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Definition of context class for MIC devices and a context manager class.
 * @author Thomas Brandes
 * @date 01.07.2013
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
