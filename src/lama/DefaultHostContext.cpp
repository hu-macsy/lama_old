/**
 * @file DefaultHostContext.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

// hpp
#include <lama/DefaultHostContext.hpp>

// others
#include <lama/exception/LAMAAssert.hpp>

#include <lama/task/TaskSyncToken.hpp>

// tracing
#include <lama/tracing.hpp>

// boost
#include <boost/bind.hpp>

using namespace boost;

namespace lama
{

LAMA_LOG_DEF_LOGGER( DefaultHostContext::logger, "Context.DefaultHostContext" );

DefaultHostContext::DefaultHostContext()
{
    mNumberOfAllocatedBytes = 0;
    mNumberOfAllocates = 0;

    LAMA_LOG_INFO( logger, "DefaultHostContext created" );
}

DefaultHostContext::~DefaultHostContext()
{
    if ( mNumberOfAllocates > 0 )
    {
        LAMA_LOG_ERROR( logger, *this << ": " << mNumberOfAllocates << " allocate without free" );
    }

    if ( mNumberOfAllocatedBytes != 0 )
    {
        LAMA_LOG_ERROR( logger,
                        *this << ": number of allocated bytes = " << mNumberOfAllocatedBytes << ", mismatch of free/allocate sizes" );
    }

    LAMA_LOG_INFO( logger, "~DefaultHostContext" );
}

void DefaultHostContext::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "DefaultHostContext";
}

void* DefaultHostContext::allocate( const size_t size ) const
{
    LAMA_ASSERT_ERROR( size > 0, "allocate with size = " << size << " should not be done" );

    void* pointer = malloc( size );

    if ( pointer == NULL )
    {
        LAMA_THROWEXCEPTION( "malloc failed for size = " << size );
    }

    // allocate must be thread-safe in case where multiple threads use LAMA arrays

    boost::recursive_mutex::scoped_lock scoped_lock( allocate_mutex );

    mNumberOfAllocatedBytes += size;
    mNumberOfAllocates++;

    LAMA_LOG_DEBUG( logger, "allocated " << pointer << ", size = " << size );

    return pointer;
}

void DefaultHostContext::allocate( ContextData& contextData, const size_t size ) const
{
    contextData.pointer = allocate( size );
}

void DefaultHostContext::free( void* pointer, const size_t size ) const
{
    LAMA_LOG_DEBUG( logger, "free " << pointer << ", size = " << size );

    LAMA_ASSERT_ERROR( mNumberOfAllocates >= 1, "Invalid Free, because there are no open allocates." );
    ::free( pointer );

    boost::recursive_mutex::scoped_lock scoped_lock( allocate_mutex );

    mNumberOfAllocatedBytes -= size;
    mNumberOfAllocates--;
}

void DefaultHostContext::free( ContextData& contextData ) const
{
    LAMA_ASSERT_EQUAL_ERROR( contextData.context->getType(), getType() );
    contextData.free();
}

void DefaultHostContext::memcpy( void* dst, const void* src, const size_t size ) const
{
    LAMA_LOG_DEBUG( logger, "memcpy: " << dst << " <- " << src << ", size = " << size );
    ::memcpy( dst, src, size );
}

std::auto_ptr<SyncToken> DefaultHostContext::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    return std::auto_ptr<SyncToken>( new TaskSyncToken( boost::bind( &::memcpy, dst, src, size ) ) );
}

bool DefaultHostContext::cancpy( const ContextData& dst, const ContextData& src ) const
{
    return dst.context->getType() == getType() && src.context->getType() == getType();
}

void DefaultHostContext::memcpy( ContextData& dst, const ContextData& src, const size_t size ) const
{
    LAMA_ASSERT_ERROR( dst.context->getType() == getType() && src.context->getType() == getType(),
                       "Can not copy from "<< *(src.context) << " to " << *(dst.context) );
    memcpy( dst.pointer, src.pointer, size );
}

std::auto_ptr<SyncToken> DefaultHostContext::memcpyAsync(
    ContextData& dst,
    const ContextData& src,
    const size_t size ) const
{
    LAMA_ASSERT_ERROR( dst.context->getType() == getType() && src.context->getType() == getType(),
                       "Can not copy from "<< *(src.context) << " to " << *(dst.context) );
    return memcpyAsync( dst.pointer, src.pointer, size );
}

} // namespace

