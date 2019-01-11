/**
 * @file Memory.cpp
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

Memory::Memory( MemoryType type ) : mMemoryType( type )
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

std::ostream& operator<<( std::ostream& stream, const MemoryType& type )
{
    switch ( type )
    {
        case MemoryType::HostMemory :
            stream << "HostMemory";
            break;

        case MemoryType::CUDAMemory :
            stream << "CUDAMemory";
            break;

        case MemoryType::CUDAHostMemory :
            stream << "CUDAHostMemory";
            break;

        case MemoryType::UserMemory :
            stream << "UserMemory";
            break;

        default:
            stream << "MemoryType_" << ( int ) type;
    }

    return stream;
}

/* ---------------------------------------------------------------------------------*/

} /* end namespace hmemo */

} /* end namespace scai */
