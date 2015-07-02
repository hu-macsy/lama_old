/**
 * @file LAMAArray.cpp
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
 * @brief Implementation of methods for class LAMAArray.
 * @author Thomas Brandes
 * @date 14.03.2011
 * @since 1.0.0
 */

// hpp
#include <memory/LAMAArray.hpp>

// others
#include <common/Exception.hpp>

#include <memory/SyncToken.hpp>

#include <memory>
#include <limits>

namespace memory
{

typedef ContextData::AccessKind AccessKind;

/* ---------------------------------------------------------------------------------*/

LAMA_LOG_DEF_LOGGER( _LAMAArray::logger, "LAMAArray" )

/* -------------------------------------------------------------------------- */

_LAMAArray* _LAMAArray::create( const Scalar::ScalarType valueType )
{
    switch( valueType )
    {
        case Scalar::FLOAT:
        {
            return new LAMAArray<float>();
        }

        case Scalar::DOUBLE:
        {
            return new LAMAArray<double>();
        }

        default:
        {
            COMMON_THROWEXCEPTION( "Unsupported ValueType " << valueType )
        }
    }
}

/* ---------------------------------------------------------------------------------*/

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, LAMAArray<ValueType>::logger, "LAMAArray" )

/* ================================================================================= */

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray()
                : _LAMAArray( 0, sizeof( ValueType ) )
{
    // just make an entry for host context

    ContextPtr host = Context::getContext( Context::Host );

    ContextDataIndex data = mContextManager.getContextData( host );

    LAMA_LOG_DEBUG( logger, "created new LAMA array, mSize = " << mSize )
    LAMA_LOG_DEBUG( logger, "created new LAMA array, mValueSize = " << mValueSize )
    
    // LAMA_LOG_DEBUG( logger, "created new context array: " << *this )
}

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n ) : _LAMAArray( n, sizeof( ValueType) )
{
    ContextPtr host = Context::getContext( Context::Host );

    // just allocate the data at host context with a corresponding write access

    size_t validSize = 0;   // no valid data availalbe, so even don't search for it

    ContextDataIndex index = mContextManager.acquireAccess( host, ContextData::Write, mSize * mValueSize, validSize );

    ContextData& data = mContextManager[index];
   
    data.releaseLock( ContextData::Write );
}

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n, const ValueType& value ) : _LAMAArray( n, sizeof( ValueType ) )

{
    // In constructor of the LAMA array lock of accesses is not required 

    ContextPtr host = Context::getContext( Context::Host );

    ContextDataIndex index = mContextManager.acquireAccess( host, ContextData::Write, mSize * mValueSize, 0 );

    ContextData& data = mContextManager[index];

    if ( n > 0 )
    {
        ValueType* host_pointer = static_cast<ValueType*>( data.pointer );

#pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

        for ( int i = 0; i < mSize; ++i )
        {
            host_pointer[i] = value;
        }
    }

    data.releaseLock( ContextData::Write );

    LAMA_LOG_DEBUG( logger, "constructed: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const LAMAArray<ValueType>& other )

    : _LAMAArray( other.mSize, sizeof( ValueType ) )
{
    mContextManager.copyAllValidEntries( other.mContextManager, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::~LAMAArray()
{
    // destructor of ContextManager does all the release/check stuff

    LAMA_LOG_DEBUG( logger, "~LAMAArray = " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Scalar::ScalarType LAMAArray<ValueType>::getValueType() const
{
    // Note: this is implementation of the pure method of base class _LAMAArray.

    return Scalar::getType<ValueType>();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>& LAMAArray<ValueType>::operator=( const LAMAArray<ValueType>& other )
{
    LAMA_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    if ( other == *this )
    {
        return *this;
    }

    mSize      = other.mSize;
    mValueSize = other.mValueSize;

    COMMON_ASSERT( !mContextManager.locked(), "assign to a locked array (read/write access)" )

    // ToDo: we might add an exception on same thread: only valid write location is copied

    COMMON_ASSERT( !other.mContextManager.writeLocked(), "assign of a write locked array" )

    mContextManager.invalidateAll();

    // Now the same stuff as in copy constructor

    mContextManager.copyAllValidEntries( other.mContextManager, mSize * mValueSize );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::assign( const LAMAArray<ValueType>& other, ContextPtr context )
{
    LAMA_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    if ( other == *this )
    {
         mContextManager.setValidData( context, mContextManager, mSize * mValueSize );
         return;
    }

    mSize      = other.mSize;
    mValueSize = other.mValueSize;

    COMMON_ASSERT( !mContextManager.locked(), "assign to a locked array (read/write access)" )

    // ToDo: we might add an exception on same thread: only valid write location is copied

    COMMON_ASSERT( !other.mContextManager.writeLocked(), "assign of a write locked array" )

    mContextManager.invalidateAll();

    mContextManager.setValidData( context, other.mContextManager, mSize * mValueSize );
 
    LAMA_LOG_DEBUG( logger, *this << " has now been assigned at " << *context )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::swap( LAMAArray<ValueType>& other )
{
    LAMA_LOG_DEBUG( logger, *this << ": swap with other = " << other )

    // we cannot swap if there is any access for any array

    COMMON_ASSERT( !other.mContextManager.locked(), "swap: other array locked" )
    COMMON_ASSERT( mContextManager.locked(), "this array locked" )

    std::swap( mSize, other.mSize );

    mContextManager.swap( other.mContextManager );

    LAMA_LOG_DEBUG( logger, *this << ": has been swapped with other = " << other )
}

/* ---------------------------------------------------------------------------------*/

void _LAMAArray::prefetch( ContextPtr context ) const
{
    mContextManager.prefetch( context, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

bool _LAMAArray::isValid( ContextPtr context ) const
{
    return mContextManager.isValid( context );
}

/* ---------------------------------------------------------------------------------*/

IndexType _LAMAArray::capacity( ContextPtr context ) const
{
    return mContextManager.capacity( context ) / mValueSize;
}

/* ---------------------------------------------------------------------------------*/

ContextPtr _LAMAArray::getValidContext( const Context::ContextType preferredType ) const
{
    return mContextManager.getValidContext( preferredType );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( ContextPtr context, const IndexType capacity )
{
    mContextManager.reserve( context, capacity * mValueSize, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear()
{
    COMMON_ASSERT( !mContextManager.locked(), "Tried to clear a locked LAMAArray " << *this )

    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::purge()
{
    mContextManager.purge();

    mSize = 0;
 
    LAMA_LOG_DEBUG( logger, *this << " purged" )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType* LAMAArray<ValueType>::get( ContextDataIndex index )
{
    return static_cast<ValueType*>( mContextManager[index].pointer );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const ValueType* LAMAArray<ValueType>::get( ContextDataIndex index ) const
{
    return static_cast<const ValueType*>( mContextManager[index].pointer );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear( const ContextDataIndex index )
{
    ContextData& data = mContextManager[index];

    COMMON_ASSERT( data.locked( ContextData::Write ), "clear illegal here " << data )

    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::resize( ContextDataIndex index, const IndexType size )
{
    ContextData& entry = mContextManager[index];

    COMMON_ASSERT( entry.locked( ContextData::Write ), "resize illegal here " << entry )

    size_t allocSize = size * mValueSize;

    size_t validSize = mSize * mValueSize;

    if ( validSize > allocSize )
    {
        validSize = allocSize;   // some entries are no more needed
    }

    LAMA_LOG_INFO( logger, *this << ": resize, needed = " << allocSize << " bytes, used = " 
                         << validSize << " bytes, capacity = " << entry.capacity() << " bytes" )

    entry.reserve( allocSize, validSize );

    // capacity is now sufficient for size elements

    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( ContextDataIndex index, const IndexType size ) const
{
    COMMON_THROWEXCEPTION( "not available yet" )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType LAMAArray<ValueType>::capacity( ContextDataIndex index ) const
{
    ContextData& entry = mContextManager[index];

    return static_cast<IndexType>( entry.capacity() / sizeof( ValueType ) );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ContextDataIndex LAMAArray<ValueType>::acquireReadAccess( ContextPtr context ) const
{
    // common::Thread::ScopedLock lock( mAccessMutex );

    LAMA_LOG_DEBUG( logger, "acquireReadAccess for " << *this );

    size_t allocSize = mSize * mValueSize;
    size_t validSize = allocSize;                   // read access needs valid data in any case

    return mContextManager.acquireAccess( context, ContextData::Read, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseReadAccess( ContextDataIndex index ) const
{
    // common::Thread::ScopedLock lock( mAccessMutex );

    mContextManager[index].releaseLock( ContextData::Read );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ContextDataIndex LAMAArray<ValueType>::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    // common::Thread::ScopedLock lock( mAccessMutex );

    size_t allocSize = mSize * sizeof( ValueType );
    size_t validSize = keepFlag ? allocSize : 0 ;

    return mContextManager.acquireAccess( context, ContextData::Write, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseWriteAccess( ContextDataIndex index )
{
    // common::Thread::ScopedLock lock( mAccessMutex );

    ContextData& data = mContextManager[index];

    data.releaseLock( ContextData::Write );

    LAMA_LOG_INFO( logger, "releaseWriteAccess: " << data );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "LAMAArray<";
    stream << Scalar::getType<ValueType>();
    stream << ">(" << mSize; stream << ")";
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( ValueType* pointer, IndexType size )
                : LAMAArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

/*
    ContextData& host = *mContextData[0];
    host.setRef( pointer, size * sizeof(ValueType) );
*/

    mSize = size;
}

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( const ValueType* pointer, IndexType size )
                : LAMAArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

/*
    ContextData& host = *mContextData[0];
    host.setRef( const_cast<ValueType*>( pointer ), size * sizeof(ValueType) );

    constFlag = true; // makes sure that we cannot have a WriteAccess

    mSize = size;
*/

}

/* ---------------------------------------------------------------------------------*/

template class COMMON_DLL_IMPORTEXPORT LAMAArray<float> ;   
template class COMMON_DLL_IMPORTEXPORT LAMAArrayRef<float> ;

template class COMMON_DLL_IMPORTEXPORT LAMAArray<double> ;   
template class COMMON_DLL_IMPORTEXPORT LAMAArrayRef<double> ;

template class COMMON_DLL_IMPORTEXPORT LAMAArray<IndexType> ;   
template class COMMON_DLL_IMPORTEXPORT LAMAArrayRef<IndexType> ;

}// namespace LAMA
