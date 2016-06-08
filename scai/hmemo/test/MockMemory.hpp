/**
 * @file MockMemory.hpp
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
 * @endlicense
 *
 * @brief Definition of a Context class used for mock objects during tests.
 * @author Thomas Brandes
 * @date 05.07.2015
 */

#include <scai/hmemo/Memory.hpp>
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/weak_ptr.hpp>

/** Exampes of a new memory class that implements all relevant routines. */

class MockMemory: public scai::hmemo::Memory
{
private: 

    // Member variables

    int mDeviceNr;     // MockMemory with different device numbers are not equal

public:

    // MockMemory uses the type UserContext as its type

    MockMemory( int deviceNr ) : scai::hmemo::Memory( scai::hmemo::memtype::UserMemory )
    {
        mDeviceNr = deviceNr;
    }

    ~MockMemory()
    {
        SCAI_LOG_DEBUG( logger, "~MockMemory: " << *this )
    }

    int getDeviceNr() const
    {
        return mDeviceNr;
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "MockMemory( dev = " << mDeviceNr << " )";
    }

    virtual scai::hmemo::ContextPtr getContextPtr() const
    {
        return scai::hmemo::Context::getContextPtr( scai::common::context::UserContext, mDeviceNr );
    }

    virtual scai::hmemo::memtype::MemoryType getType() const
    {
        return scai::hmemo::memtype::UserMemory;
    }

    virtual void* allocate( const size_t size ) const
    {
        return malloc( size );
    }

    virtual void free( void* pointer, const size_t ) const
    {
        ::free( pointer );
    }

    virtual void memcpy( void* target, const void* source, const size_t size ) const
    {
        ::memcpy( target, source, size );
    }

    virtual void memset( void* target, const int val, const size_t size ) const
    {
        ::memset( target, val, size );
    }

    static scai::tasking::SyncToken* theMemcpyAsync( void* dst, const void* src, const size_t size )
    {
        return new scai::tasking::TaskSyncToken( scai::common::bind( &::memcpy, dst, src, size ) );
    }

    virtual scai::tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const
    {
        return new scai::tasking::TaskSyncToken( scai::common::bind( &::memcpy, dst, src, size ) );
    }

    virtual bool canCopyFrom( const scai::hmemo::Memory& other ) const
    {
        // copy from host to this context should always be supported

        return other.getType() == scai::hmemo::memtype::HostMemory;
    }

    virtual bool canCopyTo( const scai::hmemo::Memory& other ) const
    {
        // copy from this context to host should always be supported

        return other.getType() == scai::hmemo::memtype::HostMemory;
    }

    virtual void memcpyFrom( void* dst, const scai::hmemo::Memory& srcMemory, const void* src, size_t size ) const
    {
        if ( srcMemory.getType() == scai::hmemo::memtype::HostMemory )
        {
            ::memcpy( dst, src, size );
        }
        else
        {
            COMMON_THROWEXCEPTION( "copy from " << srcMemory << " to " << *this << " not supported" )
        }
    }

    virtual void memcpyTo( const scai::hmemo::Memory& dstMemory, void* dst, const void* src, size_t size ) const
    {
        if ( dstMemory.getType() == scai::hmemo::memtype::HostMemory )
        {
            ::memcpy( dst, src, size );
        }
        else
        {
            COMMON_THROWEXCEPTION( "copy to " << dstMemory << " from " << *this << " not supported" )
        }
    }

    virtual scai::tasking::TaskSyncToken* getSyncToken() const
    {
        return new scai::tasking::TaskSyncToken();
    }
};

