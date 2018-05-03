/**
 * @file _HArray.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Implementations of class _HArray.
 * @author Thomas Brandes
 * @date 03.07.2015
 */

#include <scai/hmemo/_HArray.hpp>

// local library

#include <scai/hmemo/HostMemory.hpp>

namespace scai
{
namespace hmemo
{

/* ---------------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( _HArray::logger, "HArray" )

/* ---------------------------------------------------------------------------------*/

_HArray::_HArray( const _HArray& other ) :

    mSize( other.mSize ),
    mValueSize( other.mValueSize ),
    constFlag( false )                    // copy can be modified

{
    mContextDataManager.copyAllValidEntries( other.mContextDataManager, mSize * mValueSize );

    SCAI_LOG_DEBUG( logger, "_HArray copy constructor: this = " << *this );
}

/* ---------------------------------------------------------------------------------*/

_HArray::_HArray( _HArray&& other ) noexcept :

    mSize( other.mSize ),
    mValueSize( other.mValueSize ),
    constFlag( other.constFlag  ), 
    mContextDataManager( std::move( other.mContextDataManager ) )

{
    other.mSize = 0;
    other.constFlag = false;

    SCAI_LOG_DEBUG( logger, "_HArray move constructor: now other = " << other << ", this = " << *this );
}

/* ---------------------------------------------------------------------------------*/

_HArray& _HArray::operator=( const _HArray& other )
{
    SCAI_LOG_DEBUG( logger, "copy assignment, other = " << other )

    if ( &other == this )
    {
        return *this;
    }

    SCAI_ASSERT_EQ_DEBUG( mValueSize, other.mValueSize, "assignment must only be called for same value type" )

    mSize      = other.mSize;
    constFlag  = false;

    SCAI_ASSERT( !mContextDataManager.locked(), "assign to a locked array (read/write access)" )

    // ToDo: we might add an exception on same thread: only valid write location is copied
    SCAI_ASSERT( !other.mContextDataManager.locked( common::AccessKind::Write ), "assign of a write locked array" )
    mContextDataManager.invalidateAll();
    // Now the same stuff as in copy constructor
    mContextDataManager.copyAllValidEntries( other.mContextDataManager, mSize * mValueSize );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

void _HArray::assign( const _HArray& other, ContextPtr context )
{
    if ( &other == this )
    {
        mContextDataManager.setValidData( context, mContextDataManager, mSize * mValueSize );
        return;
    }

    SCAI_ASSERT_EQ_DEBUG( mValueSize, other.mValueSize, "assign must only be called for same value type" )

    SCAI_ASSERT( !mContextDataManager.locked(), "assign to a locked array (read/write access)" )

    mSize     = other.mSize;
    constFlag = false;

    SCAI_ASSERT( !other.mContextDataManager.locked( common::AccessKind::Write ), "assign of a write locked array" )
    mContextDataManager.invalidateAll();
    mContextDataManager.setValidData( context, other.mContextDataManager, mSize * mValueSize );
    SCAI_LOG_DEBUG( logger, *this << " has now been assigned at " << *context )
}

/* ---------------------------------------------------------------------------------*/

_HArray& _HArray::operator=( _HArray&& other ) noexcept
{
    // Because of the noexcept keyword, we can not throw an exception here without possibly
    // causing compiler warnings, so standard SCAI_ assertions are off the table. Since this
    // method should not be possible to call for a user, we print an error directly and
    // immediately terminate instead.
#ifdef SCAI_ASSERT_LEVEL_DEBUG
    if (mValueSize != other.mValueSize)
    {
        std::cerr << "INTERNAL ERROR: _HArray move constructor called for incompatible types." << std::endl;
        std::terminate();
    }
#endif

    SCAI_LOG_DEBUG( logger, "_HArray move assign, other = " << other << ", this = " << *this );

    mSize = other.mSize;
    constFlag = other.constFlag;

    other.mSize = 0;
    other.constFlag = false;

    mContextDataManager = std::move( other.mContextDataManager );

    SCAI_LOG_DEBUG( logger, "_HArray move assign done, other = " << other << ", this = " << *this );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

void _HArray::swap( _HArray& other )
{
    SCAI_ASSERT_EQUAL( mValueSize, other.mValueSize, "value size mismatch, different types for swap" )

    // we cannot swap if there is any access for any array

    SCAI_ASSERT_EQUAL( 0, other.mContextDataManager.locked(), "swap: other array locked: " << other )
    SCAI_ASSERT_EQUAL( 0, mContextDataManager.locked(), "this array locked: " << *this )

    mContextDataManager.swap( other.mContextDataManager );

    std::swap( mSize, other.mSize );

    SCAI_LOG_DEBUG( logger, *this << ": has been swapped with other = " << other )
}

/* ---------------------------------------------------------------------------------*/

void _HArray::setHostRef( const IndexType size, void* pointer )
{
    ContextData& host = mContextDataManager[ HostMemory::getIt() ];

    // dealing with const references in ContextData is not supported

    host.setRef( pointer, size * mValueSize );

    // Take care of const awareness by setting a flag
 
    constFlag = false;

    mSize = size;
}

void _HArray::setHostRef( const IndexType size, const void* pointer )
{
    ContextData& host = mContextDataManager[ HostMemory::getIt() ];

    // dealing with const references in ContextData is not supported

    host.setRef( const_cast<void*>( pointer ), size * mValueSize );

    // Take care of const awareness by setting a flag
 
    constFlag = true;

    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

void _HArray::reserveWithIndex( ContextDataIndex index, const IndexType size ) const
{
    if ( size <= mSize )
    {
        return;   // nothing to do
    }

    bool inUse =  mContextDataManager.locked() > 1;   // further accesses on this array
    ContextData& entry = mContextDataManager[index];

    size_t allocSize = static_cast<size_t>( size ) * static_cast<size_t>( mValueSize );
    size_t validSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );

    entry.reserve( allocSize, validSize, inUse );

    // Note: mSize does not change by the reserve
}

/* ---------------------------------------------------------------------------------*/

void _HArray::resizeWithIndex( ContextDataIndex index, const IndexType size )
{
    ContextData& entry = mContextDataManager[index];

    bool inUse =  mContextDataManager.locked() > 1;   // further accesses on this array

    // static cast to have multiplication with 64 bit values for memory addresses

    size_t allocSize = static_cast<size_t>( size )  * static_cast<size_t>( mValueSize );
    size_t validSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );

    if ( validSize > allocSize )
    {
        validSize = allocSize;   // some entries are no more needed
    }

    SCAI_LOG_INFO( logger, *this << ": resize, needed = " << allocSize << " bytes, used = "
                   << validSize << " bytes, capacity = " << entry.capacity() << " bytes" )
    entry.reserve( allocSize, validSize, inUse );
    // capacity is now sufficient for size elements
    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

void _HArray::clearWithIndex( const ContextDataIndex index )
{   
    // make sure that we have exactly one write access at this context

    ContextData& data = mContextDataManager[index];
    SCAI_ASSERT_EQUAL( 1, mContextDataManager.locked( common::AccessKind::Write ), "multiple write access for clear" << data )
    SCAI_ASSERT_EQUAL( 0, mContextDataManager.locked( common::AccessKind::Read ), "further read access, cannot clear " << data )
    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

void _HArray::writeAt( std::ostream& stream ) const
{
    // instead of value type we can only print number of bytes for one element

    stream << "_HArray<";
    stream << mValueSize;
    stream << ">(" << mSize;
    stream << ") ";
    stream << mContextDataManager;
}

void _HArray::writeAtTyped( std::ostream& stream, const char* type ) const
{
    // instead of value type we can only print number of bytes for one element

    stream << "HArray<";
    stream << type;
    stream << ">(" << mSize;
    stream << ") ";
    stream << mContextDataManager;
}

/* ---------------------------------------------------------------------------------*/

} /* end namespace hmemo */

} /* end namespace scai */
