/**
 * @file Memory.cpp
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
 * @brief Implementation of methods for base class Memory.
 * @author Thomas Brandes
 * @date 10.07.2015
 */

// hpp
#include <scai/hmemo/Memory.hpp>
#include <scai/hmemo/exception/MemoryException.hpp>

// internal scai libraries
#include <scai/common/macros/throw.hpp>

// std
#include <map>

namespace scai
{

namespace hmemo
{

/* ---------------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( Memory::logger, "Memory" )

/* ---------------------------------------------------------------------------------*/

Memory::Memory( memtype::MemoryType type ) : mMemoryType( type )
{
    SCAI_LOG_DEBUG( logger, "Memory( type = " << mMemoryType << " )" )
}

Memory::~Memory()
{
    SCAI_LOG_DEBUG( logger, "~Memory( type = " << mMemoryType << " )" )
}

/* ---------------------------------------------------------------------------------*/

void Memory::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "Memory";
}

/* ---------------------------------------------------------------------------------*/

bool Memory::canCopyFrom( const Memory& srcMemory ) const
{
    bool supported = false;

    if ( &srcMemory == this )
    {
        supported = true;
    }

    return supported;
}

bool Memory::canCopyTo( const Memory& dstMemory ) const
{
    bool supported = false;

    if ( &dstMemory == this )
    {
        supported = true;
    }

    return supported;
}

void Memory::memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const
{
    if ( &srcMemory == this )
    {
        memcpy( dst, src, size );
        return;
    }

    SCAI_THROWEXCEPTION( MemoryException,
                         "copy from " << srcMemory << " to this " << *this << " is UNSUPPORTED, "
                           << "dst = " << dst << ", src = " << src << ", size = " << size )
}

void Memory::memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const
{
    if ( &dstMemory == this )
    {
        memcpy( dst, src, size );
        return;
    }

    SCAI_THROWEXCEPTION( MemoryException,
                         "copy to " << dstMemory << " from this " << *this << " is UNSUPPORTED, "
                         << "dst = " << dst << ", src = " << src << ", size = " << size )
}

/* ---------------------------------------------------------------------------------*/

tasking::SyncToken* Memory::memcpyAsync( void* dst, const void* src, size_t size ) const
{
    memcpy( dst, src, size );
    return NULL;
}

tasking::SyncToken* Memory::memcpyFromAsync( void* dst, const Memory& srcMemory, const void* src, size_t size ) const
{
    memcpyFrom( dst, srcMemory, src, size );
    return NULL;
}

tasking::SyncToken* Memory::memcpyToAsync( const Memory& dstMemory, void* dst, const void* src, size_t size ) const
{
    memcpyTo( dstMemory, dst, src, size );
    return NULL;
}

/* ---------------------------------------------------------------------------------*/

namespace memtype
{

std::ostream& operator<<( std::ostream& stream, const MemoryType& type )
{
    switch ( type )
    {
        case HostMemory :
            stream << "HostMemory";
            break;

        case CUDAMemory :
            stream << "CUDAMemory";
            break;

        case CUDAHostMemory :
            stream << "CUDAHostMemory";
            break;

        case GPIMemory :
            stream << "GPIMemory";
            break;

        case UserMemory :
            stream << "UserMemory";
            break;

        default:
            stream << "MemoryType_" << (int) type;
    }

    return stream;
}

} /* end namespace memtype */

/* ---------------------------------------------------------------------------------*/

} /* end namespace hmemo */

} /* end namespace scai */
