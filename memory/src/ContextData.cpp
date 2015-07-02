/**
 * @file ContextData.cpp
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
 * @brief Implementation of methods for Contextdata.
 * @author Thomas Brandes
 * @date 11.07.2011
 * @since 1.0.0
 */

// hpp
#include <memory/ContextData.hpp>
#include <memory/Context.hpp>

#include <common/Exception.hpp>
#include <logging/logging.hpp>

namespace memory
{
/* ---------------------------------------------------------------------------------*/

LAMA_LOG_DEF_LOGGER( ContextData::logger, "ContextData" )

/* ---------------------------------------------------------------------------------*/

ContextData::ContextData( ContextPtr d ) :

    context( d ),
    pointer( 0 ),
    size( 0 ),
    allocated( false ),
    valid( false ),
    pinned( false ),
    mCleanFunktion( 0 )
{
    // there a no read/write locks on the context
    lock[Read] = 0;
    lock[Write] = 0;
}

/* ---------------------------------------------------------------------------------*/

ContextData::~ContextData()
{
    LAMA_LOG_DEBUG( logger, "~ContextData @ " << *context << ", size = " << size )

    // free();  
}

/* ---------------------------------------------------------------------------------*/

void ContextData::allocate( const size_t size )
{
    COMMON_ASSERT( 0 == pointer, "ContextData data already given at " << *context )
    context->allocate( *this, size );

    if ( !pointer )
    {
        COMMON_THROWEXCEPTION( "Could not allocate ContextData of size = " << size << " on " << *context )
    }

    this->size = size;
    allocated = true;
    LAMA_LOG_DEBUG( logger, "allocated " << size << " bytes" )
}

/* ---------------------------------------------------------------------------------*/

void ContextData::setRef( void* reference, const size_t size )
{
    COMMON_ASSERT( 0 == pointer, "ContextData data already given at " << *context )
    pointer = reference;
    this->size = size;
    allocated = false;

    if ( !pointer && size )
    {
        COMMON_THROWEXCEPTION( "NULL pointer cannot set be as reference, size = " << size )
    }

    LAMA_LOG_DEBUG( logger, "set ref for " << size << " bytes" )
}

/* ---------------------------------------------------------------------------------*/

void ContextData::free()
{
    // Free will/should only be called on an unlocked array
    LAMA_LOG_TRACE( logger, "free for " << *context )
    COMMON_ASSERT( 0 == lock[Read], "cannot free read locked data on " << *context )
    COMMON_ASSERT( 0 == lock[Write], "cannot free write locked data on " << *context )

    if ( context && pointer )
    {
        if ( mCleanFunktion )
        {
            mCleanFunktion( pointer );
        }

        pinned = false;

        if ( allocated )
        {
            context->free( pointer, size );
        }
    }

    pointer = 0;
    size = 0;
    valid = false;
    // we do not delete the context pointers as output will cause runtime errors
}

/* ---------------------------------------------------------------------------------*/

bool ContextData::isPinned() const
{
    return pinned;
}

/* ---------------------------------------------------------------------------------*/

void ContextData::setPinned() const
{
    pinned = true;
}

/* ---------------------------------------------------------------------------------*/

void ContextData::setCleanFunction( boost::function<void( void* )> cleanFunktion ) const
{
    mCleanFunktion = cleanFunktion;
}

/* ---------------------------------------------------------------------------------*/

void ContextData::realloc( const size_t newSize, const size_t validSize )
{
    // Note: realloc can also be used to shrink the array size

    COMMON_ASSERT( allocated, "Cannot realloc data set by reference" )
    COMMON_ASSERT( context, "no context available for realloc" )

    COMMON_ASSERT_LE( validSize, size, "size of valid data is more than actual size" )
    COMMON_ASSERT_LE( validSize, newSize, "size of valid data is more than new size" )

    void* oldPointer = pointer;

    size_t oldSize = size;

    if ( validSize <= 0 )
    {
        // existent data is no more needed
        context->free( pointer, size );
    }

    // allocate new memory

    context->allocate( *this, newSize );
    size = newSize;

    if ( validSize > 0 )
    {
        // copy the old entries in the new memory befree free of old memory

        context->memcpy( pointer, oldPointer, validSize );
        context->free( oldPointer, oldSize );
    }
}

void ContextData::reserve( const size_t newSize, const size_t validSize )
{
    if ( newSize <= size )
    {
        COMMON_ASSERT_LE( validSize, newSize, "size of valid data is more than new size" )

        // current capacity is sufficient
        return;
    }

    if ( size == 0 )
    {
        COMMON_ASSERT_LE( validSize, size, "size of valid data is more than actual size" )

        // first allocation of the data, validSize is also 0

        context->allocate( *this, newSize );
        size = newSize;
        allocated = true;
        return;
    }

    realloc( newSize, validSize );
}

void ContextData::releaseLock( AccessKind kind )
{
    COMMON_ASSERT( lock[kind] > 0,
                   "Tried to release a non existing access " << kind << " on " << context )
    --lock[kind];
}

void ContextData::writeAt( std::ostream& stream ) const
{
    stream << "ContextData ( size = " << size 
           << ", valid = " << valid << ", readLocks = " << ( ( int ) lock[Read] )
           << ", writeLocks = " << ( ( int ) lock[Write] ) << " ) @ " << *context;
}

void ContextData::copyFrom( const ContextData& other, size_t size )
{
    LAMA_LOG_INFO( logger, "copyFrom " << *other.context << " to " << *context << ", size = " << size )

    if ( *context == *other.context )
    {
        context->memcpy( pointer, other.pointer, size );
    }
    else if ( context->canCopyFrom( *other.context ) )
    {
        context->memcpyFrom( pointer, *other.context, other.pointer, size );
    }
    else if ( other.context->canCopyTo( *context ) )
    {
        other.context->memcpyTo( *context, pointer, other.pointer, size );
    }
    else
    {
        COMMON_THROWEXCEPTION( "copyFrom  "
                               << *other.context << " to " << *context << ", size = " << size  << " NOT SUPPORTED" )
        // Note: calling routine can deal with it by involving ContextData available on host
    }
}

SyncToken* ContextData::copyFromAsync( const ContextData& other, size_t size )
{
    LAMA_LOG_INFO( logger, "copyFrom " << *other.context << " to " << *context << ", size = " << size )

    if ( *context == *other.context )
    {
        return context->memcpyAsync( pointer, other.pointer, size );
    }
    else if ( context->canCopyFrom( *other.context ) )
    {
        return context->memcpyFromAsync( pointer, *other.context, other.pointer, size );
    }
    else if ( other.context->canCopyTo( *context ) )
    {
        return other.context->memcpyToAsync( *context, pointer, other.pointer, size );
    }
    else
    {
        COMMON_THROWEXCEPTION( "copyFrom  "
                               << *other.context << " to " << *context << ", size = " << size  << " NOT SUPPORTED" )

        // Note: calling routine can deal with it by involving ContextData available on host

        return NULL;  // dead code, but avoids warning
    }
}

} // namespace
