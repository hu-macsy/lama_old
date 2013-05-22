/**
 * @file Context.cpp
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
 * @brief Context.cpp
 * @author Thomas Brandes
 * @date 11.07.2011
 * @since 1.0.0
 */

// hpp
#include <lama/Context.hpp>

// others
#include <lama/LAMAInterfaceRegistry.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// tracing
#include <lama/tracing.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( Context::logger, "Context" )

Context::Context( ContextType type )
    : mContextType( type ), mEnabled( false ), mFile( NULL ), mLine( 0 )
{
    LAMA_LOG_DEBUG( logger, "Context( type = " << mContextType << " )" )
}

Context::~Context()
{
    LAMA_LOG_DEBUG( logger, "~Context( type = " << mContextType << " )" )
}

bool Context::operator==( const Context& other ) const
{
    // pointer equality gives always true.

    bool same = this == &other;

    // otherwise: both contexts must have same type and can use data of each other

    if ( !same )
    {
        same = mContextType == other.mContextType;

        LAMA_ASSERT_EQUAL_DEBUG( canUseData( other), other.canUseData( *this ) )

        same = same && canUseData( other ) && other.canUseData( *this );
    }

    LAMA_LOG_TRACE( logger, *this << " == " << other << " is " << same )
    return same;
}

bool Context::operator!=( const Context& other ) const
{
    return !operator==( other );
}

void Context::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "Context";
}

std::ostream& operator<<( std::ostream& stream, const ContextType type )
{
    if ( type == Context::Host )
    {
        stream << "Host";
    }
    else if ( type == Context::CUDA )
    {
        stream << "CUDA";
    }
    else if ( type == Context::OpenCL )
    {
        stream << "OpenCL";
    }
    else
    {
        stream << (int) type;
    }
    return stream;
}

/* -----------------------------------------------------------------------------*/

const LAMAInterface& Context::getInterface() const
{
    const LAMAInterface* lamaInterface = LAMAInterfaceRegistry::getRegistry().getInterface( getType() );

    // Registry throws an exception if no interface is available

    LAMA_ASSERT_DEBUG( lamaInterface, "No lama interface available on " << *this )

    return *lamaInterface;
}

/* -----------------------------------------------------------------------------*/

void Context::enable( const char* file, int line ) const
{
    LAMA_LOG_INFO( logger, file << "( line = " << line << ") : enable " << *this )

    if ( mEnabled )
    {
        LAMA_LOG_INFO( logger,
                       "Enable " << *this << " at " << file << " ( line = " << line << " )" << ", already enabled at " << mFile << " ( line = " << mLine << " )" )
    }

    mEnabled = true;
    mFile = file;
    mLine = line;
}

/* -----------------------------------------------------------------------------*/

void Context::disable( const char* file, int line ) const
{
    LAMA_LOG_INFO( logger, file << "( line = " << line << ") : disable " << *this )

    if ( !mEnabled )
    {
        LAMA_LOG_INFO( logger,
                       "Disable " << *this << " at " << file << " ( line = " << line << " )" << ", context was not enabled before" )
    }

    mEnabled = false;
    mFile = NULL;
    mLine = 0;
}

/* ---------------------------------------------------------------------------------*/

Context::ContextData::ContextData( ContextPtr d )
    : context( d ), pointer( 0 ), size( 0 ), allocated( false ), valid( false ), pinned( false ), mCleanFunktion(
        0 )
{
    // there a no read/write locks on the context

    lock[Read] = 0;
    lock[Write] = 0;
}

/* ---------------------------------------------------------------------------------*/

Context::ContextData::~ContextData()
{
    // free allocated data

    free();
}

/* ---------------------------------------------------------------------------------*/

void Context::ContextData::allocate( const size_t size )
{
    LAMA_ASSERT_DEBUG( 0 == pointer, "ContextData data already given at " << *context )

    this->size = size;

    context->allocate( *this, size );

    if ( !pointer )
    {
        LAMA_THROWEXCEPTION( "Could not allocate ContextData of size = " << size << " on " << *context )
    }

    allocated = true;

    LAMA_LOG_DEBUG( logger, "allocated " << size << " bytes" )
}

/* ---------------------------------------------------------------------------------*/

void Context::ContextData::setRef( void* reference, const size_t size )
{
    LAMA_ASSERT_DEBUG( 0 == pointer, "ContextData data already given at " << *context )

    pointer = reference;

    this->size = size;

    allocated = false;

    if ( !pointer && size )
    {
        LAMA_THROWEXCEPTION( "NULL pointer cannot set be as reference, size = " << size )
    }

    LAMA_LOG_DEBUG( logger, "set ref for " << size << " bytes" )
}

/* ---------------------------------------------------------------------------------*/

void Context::ContextData::free()
{
    // Free will/should only be called on an unlocked array

    LAMA_LOG_TRACE( logger, "free for " << *context )

    LAMA_ASSERT_DEBUG( 0 == lock[Read], "cannot free read locked data on " << *context )
    LAMA_ASSERT_DEBUG( 0 == lock[Write], "cannot free write locked data on " << *context )

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

bool Context::ContextData::isPinned() const
{
    return pinned;
}

/* ---------------------------------------------------------------------------------*/

void Context::ContextData::setPinned() const
{
    pinned = true;
}

/* ---------------------------------------------------------------------------------*/

void Context::ContextData::setCleanFunction( boost::function<void( void* )> cleanFunktion ) const
{
    mCleanFunktion = cleanFunktion;
}

/* ---------------------------------------------------------------------------------*/

void Context::ContextData::realloc( const size_t newSize, const size_t saveSize )
{
    LAMA_ASSERT_ERROR( allocated, "Cannot realloc data set by reference" )

    // Note: realloc cn also be used to shrink the array size
    LAMA_ASSERT_DEBUG( context, "no context available for realloc" )
    LAMA_ASSERT_DEBUG( saveSize <= size, "cannot save more than current size" )
    LAMA_ASSERT_DEBUG( saveSize <= newSize, "cannot save more than new size" )
    void* oldPointer = pointer;
    size_t oldSize = size;

    if ( saveSize <= 0 )
    {
        context->free( pointer, size );
    }

    // allocate new memory
    context->allocate( *this, newSize );
    size = newSize;

    if ( saveSize > 0 )
    {
        // copy the old entries in the new memory befree free of old memory
        context->memcpy( pointer, oldPointer, saveSize );
        context->free( oldPointer, oldSize );
    }
}

}
