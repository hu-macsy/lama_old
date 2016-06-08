/**
 * @file GPIMemory.hpp
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
 * @brief Host context where memory is pinned for fast transfer via remote read/write
 * @author Thomas Brandes
 * @date 16.05.2014
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
