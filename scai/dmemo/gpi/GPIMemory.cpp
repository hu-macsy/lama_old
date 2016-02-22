/**
 * @file GPIMemory.cpp
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
 * @brief Definition of a GASPI communication context that allocates pinned memory so that
 *        copy in and copy out to GASPI segments is no more needed
 * @author Thomas Brandes
 * @date 16.05.2014
 * @since 1.0.0
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
