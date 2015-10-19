/**
 * @file MockMemory.hpp
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
 * @brief Definition of a Context class used for mock objects during tests.
 * @author: Thomas Brandes
 * @date 05.07.2015
 **/

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
        return scai::hmemo::Context::getContextPtr( scai::hmemo::context::UserContext, mDeviceNr );
    }

    virtual scai::hmemo::MemoryType getType() const
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

