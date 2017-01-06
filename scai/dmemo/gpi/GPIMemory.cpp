/**
 * @file GPIMemory.cpp
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
 * @brief Definition of a GASPI communication context that allocates pinned memory so that
 *        copy in and copy out to GASPI segments is no more needed
 * @author Thomas Brandes
 * @date 16.05.2014
 */

// hpp
#include <scai/dmemo/gpi/GPIMemory.hpp>

#include <scai/dmemo/gpi/GPIMemManager.hpp>
#include <scai/dmemo/gpi/GPIUtils.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/tasking/NoSyncToken.hpp>

// tracing
#include <scai/tracing.hpp>

#include <cstring>

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GPIMemory::logger, "Memory.GPIMemory" );

GPIMemory::GPIMemory() : hmemo::Memory( hmemo::memtype::GPIMemory )
{
    SCAI_LOG_DEBUG( logger, "GPIMemory created, allocates data as GASPI segment" )
}

GPIMemory::~GPIMemory()
{
    SCAI_LOG_DEBUG( logger, "~GPIMemory" )
}

void GPIMemory::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "GPIMemory";
}

void* GPIMemory::allocate( const size_t size ) const
{
    SCAI_REGION( "GPIMemory::allocate" )
    SCAI_LOG_DEBUG( logger, *this << ": allocate " << size << " bytes" )
    void* pointer = 0;
    gaspi_segment_id_t id;
    gaspi_pointer_t ptr;
    gaspi_offset_t offset;
    GPIMemManager::getSegmentData( id, ptr, offset, size );
    pointer = ptr;
    SCAI_LOG_DEBUG( logger, *this << ": allocated " << size << " bytes, pointer = " << pointer )
    // here only for first test
    GPIMemManager::findSegment( id, offset, ptr );
    const int outID = static_cast<int>( id );  // avoid output as char
    SCAI_LOG_DEBUG( logger, *this << ": ptr " << ptr << " belongs to segment " << outID << ", offset = " << offset )
    return pointer;
}

void GPIMemory::free( void* pointer, const size_t size ) const
{
    SCAI_REGION( "GPIMemory::free" )
    gaspi_segment_id_t id;
    gaspi_offset_t offset;
    bool found = GPIMemManager::findSegment( id, offset, pointer );
    SCAI_ASSERT_ERROR( found, "free of non GASPI segment data" )
    GPIMemManager::releaseSegmentData( id, offset );
    SCAI_LOG_DEBUG( logger, *this << ": freed " << size << " bytes, pointer = " << pointer )
}

void GPIMemory::memcpy( void* dst, const void* src, const size_t size ) const
{
    ::memcpy( dst, src, size );
}

tasking::SyncToken* GPIMemory::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    ::memcpy( dst, src, size );
    return new tasking::NoSyncToken();
}

bool GPIMemory::canCopyFrom( const Memory& other ) const
{
    bool supported = false;
    hmemo::memtype::MemoryType otherType = other.getType();

    if ( otherType == hmemo::memtype::HostMemory )
    {
        // GPIMem -> Host is supported
        supported = true;
    }
    else if ( otherType == hmemo::memtype::GPIMemory )
    {
        supported = true;
    }

    return supported;
}

bool GPIMemory::canCopyTo( const Memory& other ) const
{
    // copy from - copy to has symmetric behavior
    return canCopyFrom( other );
}

void GPIMemory::memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const
{
    SCAI_ASSERT_DEBUG( canCopyFrom( srcMemory ), "Cannot copy from " << srcMemory << " to " << *this )
    ::memcpy( dst, src, size );
}

void GPIMemory::memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const
{
    SCAI_ASSERT_DEBUG( canCopyTo( dstMemory ), "Cannot copy to " << dstMemory << " from " << *this )
    ::memcpy( dst, src, size );
}

hmemo::ContextPtr GPIMemory::getContextPtr() const
{
    // Host device does all operations on GPI memory
    return hmemo::Context::getHostPtr();
}

} // namespace dmemo

} // namespace scai
