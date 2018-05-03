/**
 * @file HostMemory.cpp
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
    mHostContextPtr( hostContextPtr ),
    mNumberOfAllocates( 0 ),
    mNumberOfAllocatedBytes( 0 ),
    mMaxAllocatedBytes( 0 )

{
    SCAI_LOG_INFO( logger, "HostMemory created" )
}

HostMemory::~HostMemory()
{
    if ( mNumberOfAllocates > 0 )
    {
        SCAI_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " allocate without free" )
    }

    if ( mNumberOfAllocatedBytes != 0 )
    {
        SCAI_LOG_ERROR( logger,
                        *this << ": number of allocated bytes = " << mNumberOfAllocatedBytes
                        << ", should be 0, so mismatch of free/allocate sizes" )
    }

    SCAI_LOG_INFO( logger, "~HostMemory" )
}

void HostMemory::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "HostMemory( " << *mHostContextPtr << " )";
}

void* HostMemory::allocate( const size_t size ) const
{
    SCAI_ASSERT( size > 0, "allocate with size = " << size << " should not be done" )
    void* pointer = malloc( size );

    if ( pointer == NULL )
    {
        SCAI_THROWEXCEPTION( MemoryException, "malloc failed for size = " << size )
    }

    // allocate must be thread-safe in case where multiple threads use LAMA arrays
    std::unique_lock<std::recursive_mutex> lock( allocate_mutex );

    mNumberOfAllocatedBytes += size;
    mMaxAllocatedBytes = std::max( mNumberOfAllocatedBytes, mMaxAllocatedBytes );
    mNumberOfAllocates++;

    SCAI_LOG_DEBUG( logger, "allocated " << pointer << ", size = " << size )
    return pointer;
}

void HostMemory::free( void* pointer, const size_t size ) const
{
    SCAI_LOG_DEBUG( logger, "free " << pointer << ", size = " << size )
    SCAI_ASSERT( mNumberOfAllocates >= 1, "Invalid free, because there are no open allocates." )
    ::free( pointer );
    std::unique_lock<std::recursive_mutex> lock( allocate_mutex );
    mNumberOfAllocatedBytes -= size;
    mNumberOfAllocates--;
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

size_t HostMemory::maxAllocatedBytes() const
{
    return mMaxAllocatedBytes;
}

} /* end namespace hmemo */

} /* end namespace scai */
