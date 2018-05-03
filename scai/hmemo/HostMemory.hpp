/**
 * @file HostMemory.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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

#include <memory>
#include <mutex>

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

    HostMemory( std::shared_ptr<const class HostContext> hostContext );

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

    /** Implementation of pure method Memory::maxAllocatedBytes() */

    virtual size_t maxAllocatedBytes() const;

private:

    std::shared_ptr<const HostContext> mHostContextPtr;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    mutable size_t mNumberOfAllocates; //!< variable counts allocates

    mutable size_t mNumberOfAllocatedBytes;//!< variable counts allocated bytes

    mutable size_t mMaxAllocatedBytes;//!< variable counts max allocated bytes

    mutable std::recursive_mutex allocate_mutex;// needed to make allocate/free thread-safe
};

} /* end namespace hmemo */

} /* end namespace scai */
