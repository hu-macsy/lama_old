/**
 * @file HostMemory.cpp
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
 * @brief Implementation of methods for Host memory by using stdlib routines.
 * @author Thomas Brandes
 * @date 17.07.2015
 */

// hpp
#include <scai/hmemo/HostMemory.hpp>
#include <scai/hmemo/exception/MemoryException.hpp>

// local library
#include <scai/hmemo/HostContext.hpp>
#include <scai/hmemo/Context.hpp>

// internal scai libraries

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/safer_memcpy.hpp>

// std
#include <cstring>
#include <functional>

using scai::common::safer_memcpy;

namespace scai
{

namespace hmemo
{

SCAI_LOG_DEF_LOGGER( HostMemory::logger, "Memory.HostMemory" )

HostMemory::HostMemory( std::shared_ptr<const HostContext> hostContextPtr ) :

    Memory( MemoryType::HostMemory ),
    mHostContextPtr( hostContextPtr )
{
    SCAI_LOG_INFO( logger, "HostMemory created" )
}

HostMemory::~HostMemory()
{
    Memory::checkAllFreed();
}

void HostMemory::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "HostMemory( " << *mHostContextPtr << " )";
}

void* HostMemory::allocate( const size_t size )
{
    SCAI_ASSERT( size > 0, "allocate with size = " << size << " should not be done" )
    void* pointer = malloc( size );

    if ( pointer == NULL )
    {
        SCAI_THROWEXCEPTION( MemoryException, "malloc failed for size = " << size )
    }

    // allocate must be thread-safe in case where multiple threads use LAMA arrays
    std::unique_lock<std::recursive_mutex> lock( allocate_mutex );

    Memory::setAllocated( size );

    SCAI_LOG_DEBUG( logger, "allocated " << pointer << ", size = " << size )
    return pointer;
}

void HostMemory::free( void* pointer, const size_t size )
{
    SCAI_LOG_DEBUG( logger, "free " << pointer << ", size = " << size )

    ::free( pointer );
    std::unique_lock<std::recursive_mutex> lock( allocate_mutex );
    Memory::setFreed( size );
}

void HostMemory::memcpy( void* dst, const void* src, const size_t size ) const
{
    SCAI_REGION( "Memory.Host_memcpy" )
    SCAI_LOG_DEBUG( logger, "memcpy: " << dst << " <- " << src << ", size = " << size )
    safer_memcpy( dst, src, size );
}

void HostMemory::memset( void* dst, const int val, const size_t size ) const
{
    SCAI_REGION( "Memory.Host_memset" )
    ::memset( dst, val, size );
}

tasking::SyncToken* HostMemory::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    return new tasking::TaskSyncToken( std::bind( &::memcpy, dst, src, size ) );
}

ContextPtr HostMemory::getContextPtr() const
{
    return mHostContextPtr;
}

MemoryPtr HostMemory::getIt()
{
    static std::shared_ptr<HostMemory> instancePtr;

    if ( !instancePtr.get() )
    {
        SCAI_LOG_DEBUG( logger, "Create instance for HostMemory" )
        ContextPtr contextPtr = Context::getContextPtr( common::ContextType::Host );
        std::shared_ptr<const HostContext> hostContextPtr = std::dynamic_pointer_cast<const HostContext>( contextPtr );
        SCAI_ASSERT( hostContextPtr.get(), "Serious: dynamic cast failed" )
        instancePtr.reset( new HostMemory( hostContextPtr ) );
    }

    return instancePtr;
}

} /* end namespace hmemo */

} /* end namespace scai */
