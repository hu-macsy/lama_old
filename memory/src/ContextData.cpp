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
    // free allocated data

    LAMA_LOG_INFO( logger, "~ContextData @ " << context << ", size = " << size )

    free();
}

/* ---------------------------------------------------------------------------------*/

void ContextData::allocate( const size_t size )
{
    COMMON_ASSERT( 0 == pointer, "ContextData data already given at " << *context )

    context->allocate( *this, size );

    if( !pointer )
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

    if( !pointer && size )
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

    if( context && pointer )
    {
        if( mCleanFunktion )
        {
            mCleanFunktion( pointer );
        }

        pinned = false;

        if( allocated )
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

void ContextData::realloc( const size_t newSize, const size_t saveSize )
{
    COMMON_ASSERT( allocated, "Cannot realloc data set by reference" )

    // Note: realloc cn also be used to shrink the array size
    COMMON_ASSERT( context, "no context available for realloc" )
    COMMON_ASSERT( saveSize <= size, "cannot save more than current size" )
    COMMON_ASSERT( saveSize <= newSize, "cannot save more than new size" )

    void* oldPointer = pointer;
    size_t oldSize = size;

    if( saveSize <= 0 )
    {
        context->free( pointer, size );
    }

    // allocate new memory

    context->allocate( *this, newSize );
    size = newSize;

    if( saveSize > 0 )
    {
        // copy the old entries in the new memory befree free of old memory
        context->memcpy( pointer, oldPointer, saveSize );
        context->free( oldPointer, oldSize );
    }
}

void ContextData::reserve( const size_t newSize, const size_t validSize )
{
    COMMON_ASSERT( validSize <= size, "cannot more data be valid than allocated" )

    if ( newSize <= size )
    {
        // current capacity is sufficient

        return;
    }

    if ( size == 0 )
    {
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
    stream << "ContextData( size = " << size << ", readLocks = " << ( (int) lock[Read] ) 
           << ", writeLocks = " << ( (int) lock[Write] ) << " )";
}

} // namespace
