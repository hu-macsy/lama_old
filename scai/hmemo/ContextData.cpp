/**
 * @file ContextData.cpp
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
 * @brief Implementation of methods for Contextdata.
 * @author Thomas Brandes
 * @date 11.07.2011
 */

// hpp
#include <scai/hmemo/ContextData.hpp>
#include <scai/hmemo/exception/MemoryException.hpp>

// local library
#include <scai/hmemo/Context.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/macros/assert.hpp>

namespace scai
{

using tasking::SyncToken;

namespace hmemo
{

/* ---------------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( ContextData::logger, "ContextData" )

/* ---------------------------------------------------------------------------------*/

ContextData::ContextData( MemoryPtr memory ) :

    mSize( 0 ),
    mMemory( memory ),
    mPointer( 0 ),
    mValid( false ),
    mAllocated( false )
{
}

/* ---------------------------------------------------------------------------------*/

ContextData::~ContextData()
{
    SCAI_LOG_DEBUG( logger, "~ContextData @ " << *mMemory << ", size = " << mSize )
    free();
}

/* ---------------------------------------------------------------------------------*/

ContextData::ContextData( ContextData&& other ) noexcept
{
    mSize = other.mSize;
    mMemory = other.mMemory;
    mPointer = other.mPointer;
    mValid = other.mValid;
    mAllocated = other.mAllocated;

    // set the data of other in such a way that destructor will not damage anything

    other.mSize = 0;
    other.mPointer = nullptr;
    other.mValid = false;
}

/* ---------------------------------------------------------------------------------*/

void ContextData::allocate( const size_t size )
{
    SCAI_ASSERT_ERROR( mPointer == nullptr, "already allocated at " << *mMemory )

    mPointer = mMemory->allocate( size );

    if ( !mPointer )
    {
        SCAI_THROWEXCEPTION( MemoryException,
                             "Could not allocate ContextData of size = " << size << " on " << *mMemory )
    }

    mSize = size;
    mAllocated = true;
    SCAI_LOG_DEBUG( logger, "allocated " << mSize << " bytes" )
}

/* ---------------------------------------------------------------------------------*/

void ContextData::setRef( void* reference, const size_t size )
{
    // we assume that reference pointer really belongs to same memory/context

    SCAI_ASSERT_ERROR( mPointer == nullptr, "setRef, but already set context data " << *mMemory )

    mPointer = reference;
    mSize = size;
    mAllocated = false;
    mValid     = true;   // we assume it as valid data

    if ( !mPointer && mSize )
    {
        COMMON_THROWEXCEPTION( "null pointer cannot set be as reference, size = " << size )
    }

    SCAI_LOG_DEBUG( logger, "set ref for " << size << " bytes" )
}

/* ---------------------------------------------------------------------------------*/

void ContextData::free()
{
    if ( mMemory && mPointer )
    {
        if ( mAllocated )
        {
            SCAI_LOG_DEBUG( logger, "ContextData: free " << mSize )
            mMemory->free( mPointer, mSize );
        }
        else
        {
            // data was only referenced
            SCAI_LOG_DEBUG( logger, "ContextData: no more ref to " << mSize )
        }
    }

    mPointer = nullptr;
    mSize = 0;
    mValid = false;

    // we do not delete the mMemory pointers as output might cause runtime errors
}

/* ---------------------------------------------------------------------------------*/

void ContextData::realloc( const size_t newSize, const size_t validSize )
{
    // Note: realloc can also be used to shrink the array size

    SCAI_ASSERT_ERROR( mAllocated, "Cannot realloc data set by reference" )
    SCAI_ASSERT_ERROR( mMemory, "no mMemory available for realloc" )
    SCAI_ASSERT_LE_ERROR( validSize, mSize, "size of valid data is more than actual size" )
    SCAI_ASSERT_LE_ERROR( validSize, newSize, "size of valid data is more than new size" )

    void* oldPointer = mPointer;
    size_t oldSize = mSize;

    if ( validSize == 0 )
    {
        // existent data is no more needed
        mMemory->free( mPointer, mSize );
    }

    // allocate new memory

    mPointer = mMemory->allocate( newSize );

    mSize = newSize;

    if ( validSize > 0 )
    {
        // copy the old entries in the new memory befree free of old memory
        mMemory->memcpy( mPointer, oldPointer, validSize );
        mMemory->free( oldPointer, oldSize );
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextData::reserve( const size_t newSize, const size_t validSize, bool inUse )
{
    if ( newSize <= mSize )
    {
        SCAI_ASSERT_LE( validSize, newSize, "size of valid data is more than new size" )
        // current capacity is sufficient
        return;
    }

    SCAI_ASSERT( !inUse, "reserve/reallocate required on array that is already in use" )

    if ( mSize == 0 )
    {
        SCAI_ASSERT_LE( validSize, mSize, "size of valid data is more than actual size" )
        // first allocation of the data, validSize is also 0
        mPointer = mMemory->allocate( newSize );
        mSize = newSize;
        mAllocated = true;
        return;
    }

    realloc( newSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

void ContextData::writeAt( std::ostream& stream ) const
{
    stream << "ContextData ( ";

    if ( mSize == 0 )
    {
        stream << "null";
    }
    else if ( mAllocated )
    {
        stream << "allocated, size = " << mSize;
    }
    else
    {
        stream << "ref, size = " << mSize;
    }

    stream << ", valid = " << mValid << " ) @ " << *mMemory;
}

/* ---------------------------------------------------------------------------------*/

void ContextData::copyFrom( const ContextData& other, size_t size )
{
    SCAI_LOG_INFO( logger, "copyFrom " << *other.mMemory << " to " << *mMemory << ", size = " << size )

    if ( mMemory.get() == other.mMemory.get() )
    {
        SCAI_LOG_INFO( logger, "copy on same context" )
        mMemory->memcpy( mPointer, other.mPointer, size );
    }
    else if ( mMemory->canCopyFrom( *other.mMemory ) )
    {
        SCAI_LOG_INFO( logger, "copy from" )
        mMemory->memcpyFrom( mPointer, *other.mMemory, other.mPointer, size );
    }
    else if ( other.mMemory->canCopyTo( *mMemory ) )
    {
        SCAI_LOG_INFO( logger, "copy to" )
        other.mMemory->memcpyTo( *mMemory, mPointer, other.mPointer, size );
    }
    else
    {
        SCAI_LOG_WARN( logger, "copyFrom " << *other.mMemory << " to " << *mMemory << " UNSUPPORTED" )
        COMMON_THROWEXCEPTION( "copyFrom  "
                               << *other.mMemory << " to " << *mMemory << ", size = " << size  << " UNSUPPORTED" )
        // Note: calling routine can deal with it by involving ContextData available on host
    }
}

/* ---------------------------------------------------------------------------------*/

SyncToken* ContextData::copyFromAsync( const ContextData& other, size_t size )
{
    SCAI_ASSERT_GE_DEBUG( mSize, size, "insufficient allocated memory to copy data from other context to here" )

    SyncToken* token = nullptr;    // default value avoids compiler warning due to exception

    SCAI_LOG_INFO( logger, "copyFrom " << *other.mMemory << " to " << *mMemory << ", size = " << size )

    if ( mMemory.get() == other.mMemory.get() )
    {
        // pointer equality implies it is the same context
        token = mMemory->memcpyAsync( mPointer, other.mPointer, size );
    }
    else if ( mMemory->canCopyFrom( *other.mMemory ) )
    {
        token = mMemory->memcpyFromAsync( mPointer, *other.mMemory, other.mPointer, size );
    }
    else if ( other.mMemory->canCopyTo( *mMemory ) )
    {
        token = other.mMemory->memcpyToAsync( *mMemory, mPointer, other.mPointer, size );
    }
    else
    {
        COMMON_THROWEXCEPTION( "copyFrom  "
                               << *other.mMemory << " to " << *mMemory << ", size = " << size  << " NOT SUPPORTED" )
        // Note: calling routine can deal with it by involving ContextData available on host
    }

    return token;
}

} /* end namespace hmemo */

} /* end namespace scai */
