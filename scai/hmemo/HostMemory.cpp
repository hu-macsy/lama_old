/**
 * @file HostMemory.cpp
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

#include <scai/common/macros/assert.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/bind.hpp>

// std
#include <cstring>

namespace scai
{

namespace hmemo
{

SCAI_LOG_DEF_LOGGER( HostMemory::logger, "Memory.HostMemory" )

HostMemory::HostMemory( common::shared_ptr<const HostContext> hostContextPtr ) : 

    Memory( memtype::HostMemory ),
    mHostContextPtr( hostContextPtr )
{
    mNumberOfAllocatedBytes = 0;
    mNumberOfAllocates = 0;

    SCAI_LOG_INFO( logger, "HostMemory created" )
}

HostMemory::~HostMemory()
{
    if( mNumberOfAllocates > 0 )
    {
        SCAI_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " allocate without free" )
    }

    if( mNumberOfAllocatedBytes != 0 )
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

    if( pointer == NULL )
    {
        SCAI_THROWEXCEPTION( MemoryException, "malloc failed for size = " << size )
    }

    // allocate must be thread-safe in case where multiple threads use LAMA arrays

    common::Thread::ScopedLock lock( allocate_mutex );

    mNumberOfAllocatedBytes += size;
    mNumberOfAllocates++;

    SCAI_LOG_DEBUG( logger, "allocated " << pointer << ", size = " << size )

    return pointer;
}

void HostMemory::free( void* pointer, const size_t size ) const
{
    SCAI_LOG_DEBUG( logger, "free " << pointer << ", size = " << size )

    SCAI_ASSERT( mNumberOfAllocates >= 1, "Invalid free, because there are no open allocates." )

    ::free( pointer );

    common::Thread::ScopedLock lock( allocate_mutex );

    mNumberOfAllocatedBytes -= size;
    mNumberOfAllocates--;
}

void HostMemory::memcpy( void* dst, const void* src, const size_t size ) const
{
    SCAI_REGION( "Memory.Host_memcpy" )

    SCAI_LOG_DEBUG( logger, "memcpy: " << dst << " <- " << src << ", size = " << size )

    ::memcpy( dst, src, size );
}

tasking::SyncToken* HostMemory::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    return new tasking::TaskSyncToken( common::bind( &::memcpy, dst, src, size ) );
}

ContextPtr HostMemory::getContextPtr() const
{
    return mHostContextPtr;
}

MemoryPtr HostMemory::getIt()
{
    static common::shared_ptr<HostMemory> instancePtr;

    if ( !instancePtr.get() )
    {
        SCAI_LOG_DEBUG( logger, "Create instance for HostMemory" ) 

        ContextPtr contextPtr = Context::getContextPtr( common::context::Host );
        common::shared_ptr<const HostContext> hostContextPtr = common::dynamic_pointer_cast<const HostContext>( contextPtr );
        SCAI_ASSERT( hostContextPtr.get(), "Serious: dynamic cast failed" )
        instancePtr.reset( new HostMemory( hostContextPtr ) );
    }

    return instancePtr;
}

} /* end namespace hmemo */

} /* end namespace scai */
