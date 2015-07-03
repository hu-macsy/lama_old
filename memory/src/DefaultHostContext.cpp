/**
 * @file DefaultHostContext.cpp
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
 * @brief DefaultHostContext.cpp
 * @author Thomas Brandes
 * @date 11.07.2011
 */

// hpp
#include <memory/DefaultHostContext.hpp>

// others
#include <common/Exception.hpp>
#include <common/OpenMP.hpp>

#include <memory/TaskSyncToken.hpp>

// boost
#include <boost/bind.hpp>

using namespace boost;

namespace memory
{

LAMA_LOG_DEF_LOGGER( DefaultHostContext::logger, "Context.DefaultHostContext" )

DefaultHostContext::DefaultHostContext()
{
    mNumberOfAllocatedBytes = 0;
    mNumberOfAllocates = 0;

    LAMA_LOG_INFO( logger, "DefaultHostContext created" )
}

DefaultHostContext::~DefaultHostContext()
{
    if( mNumberOfAllocates > 0 )
    {
        LAMA_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " allocate without free" )
    }

    if( mNumberOfAllocatedBytes != 0 )
    {
        LAMA_LOG_ERROR( logger,
                        *this << ": number of allocated bytes = " << mNumberOfAllocatedBytes 
                         << ", should be 0, so mismatch of free/allocate sizes" )
    }

    LAMA_LOG_INFO( logger, "~DefaultHostContext" )
}

void DefaultHostContext::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    int nThreads = 1;

#pragma omp parallel
    {
#pragma omp master
        {
            nThreads = omp_get_num_threads();
        }
    }

    stream << "DefaultHostContext( #Threads = " << nThreads << " )";
}

void* DefaultHostContext::allocate( const size_t size ) const
{
    COMMON_ASSERT( size > 0, "allocate with size = " << size << " should not be done" )

    void* pointer = malloc( size );

    if( pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "malloc failed for size = " << size )
    }

    // allocate must be thread-safe in case where multiple threads use LAMA arrays

    common::Thread::ScopedLock lock( allocate_mutex );

    mNumberOfAllocatedBytes += size;
    mNumberOfAllocates++;

    LAMA_LOG_DEBUG( logger, "allocated " << pointer << ", size = " << size )

    return pointer;
}

void DefaultHostContext::free( void* pointer, const size_t size ) const
{
    LAMA_LOG_DEBUG( logger, "free " << pointer << ", size = " << size )

    COMMON_ASSERT( mNumberOfAllocates >= 1, "Invalid free, because there are no open allocates." )

    ::free( pointer );

    common::Thread::ScopedLock lock( allocate_mutex );

    mNumberOfAllocatedBytes -= size;
    mNumberOfAllocates--;
}

void DefaultHostContext::memcpy( void* dst, const void* src, const size_t size ) const
{
    LAMA_LOG_DEBUG( logger, "memcpy: " << dst << " <- " << src << ", size = " << size )

    ::memcpy( dst, src, size );
}

SyncToken* DefaultHostContext::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    return new TaskSyncToken( boost::bind( &::memcpy, dst, src, size ) );
}

static boost::weak_ptr<class DefaultHostContext> contextInstance;

ContextPtr DefaultHostContext::getContext( int deviceNr )
{
    boost::shared_ptr<DefaultHostContext> context;

    if ( deviceNr >= 0 )
    {
        LAMA_LOG_WARN( logger, "Context number ignored for DefaultHostContext, deviceNr = " << deviceNr )
    }

    // use the last contextInstance if it is still valid

    if( contextInstance.expired() )
    {
        // create a new instance of DefaultHostContext and keep it for further uses

        context = boost::shared_ptr<DefaultHostContext>( new DefaultHostContext() );

        contextInstance = context;
    }
    else
    {
        // the last context instance is still valid, so we return new shared pointer to it

        context = contextInstance.lock();
    }

    return context;
}

bool DefaultHostContext::init()
{
    Context::addCreator( Context::Host, &DefaultHostContext::getContext );
}

bool DefaultHostContext::initialized = DefaultHostContext::init();

} // namespace

