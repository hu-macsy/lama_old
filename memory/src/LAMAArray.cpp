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

// Help function needed due to bug in Intel compiler

static size_t maxSize()
{
    return std::numeric_limits<size_t>::max();
}

size_t _LAMAArray::nContextIndex = maxSize();

// size_t _LAMAArray::nContextIndex = std::numeric_limits<size_t>::max();

/* ================================================================================= */

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray()
                : _LAMAArray( 0, sizeof( ValueType ) )
{
    // just make an entry for host context

    ContextPtr host = Context::getContext( Context::Host );

    ContextDataRef data = mContextManager.getContextData( host, ContextData::Write );

    LAMA_LOG_DEBUG( logger, "created new context array: " << *this )
}

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n ) : _LAMAArray( n, sizeof( ValueType) )
{
    ContextPtr host = Context::getContext( Context::Host );

    // just allocate the data at host context with a corresponding write access

    size_t validSize = 0;   // no valid data availalbe, so even don't search for it

    ContextDataRef index = mContextManager.acquireAccess( host, ContextData::Write, n * sizeof( ValueType ), validSize );

    ContextData& data = mContextManager[index];
   
    data.releaseLock( ContextData::Write );
}

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n, const ValueType& value ) : _LAMAArray( n, sizeof( ValueType ) )

{
    // In constructor of the LAMA array lock of accesses is not required 

    ContextPtr host = Context::getContext( Context::Host );

    ContextDataRef index = mContextManager.acquireAccess( host, ContextData::Write, n * sizeof( ValueType ), 0 );

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
    ContextPtr host = Context::getContext( Context::Host );

    ContextData& data = mContextManager[ mContextManager.getContextData( host, ContextData::Write ) ];

    // find a valid location of the other array in any context

    LAMA_LOG_DEBUG( logger, other << " will be copied to " << *this )

    COMMON_THROWEXCEPTION( "not available yet" )

    operator=( other );
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

    COMMON_ASSERT( this != &other, "self assign for LAMAArray forbidden" )

/*
    wait();
    mSize = other.mSize;
    other.wait();

    // this array should have no read/write access

    for( size_t i = 0; i < mContextData.size(); i++ )
    {
        ContextData& entry = *mContextData[i];

        if( entry.lock[ContextData::Read] || entry.lock[ContextData::Write] )
        {
            COMMON_THROWEXCEPTION( *this << ": cannot copy/assign to locked array" )
        }

        // make all entries invalid because we will overwrite them

        if( entry.valid )
        {
            entry.valid = false;
            LAMA_LOG_DEBUG( logger, *this << ": invalidated at index = " << i << " for " << entry.context )
        }
    }

    // each valid data of the other array will be copied into the same context for this array
    size_t nOtherContexts = other.mContextData.size();

    for( size_t i = 0; i < nOtherContexts; i++ )
    {
        LAMA_LOG_TRACE( logger, "Other Context index = " << i )
        const ContextData& otherEntry = *other.mContextData[i];

        if( !otherEntry.valid )
        {
            LAMA_LOG_TRACE( logger,
                            other << ": context " << i << " of " << nOtherContexts << " for " << *otherEntry.context << " is invalid, will not be copied" )
            continue; // do not copy any invalid data
        }

        size_t contextIndex;
        size_t validIndex;
        getAccess( contextIndex, validIndex, otherEntry.context, ContextData::Write );
        LAMA_LOG_TRACE( logger,
                        *this << ": write access for " << *otherEntry.context << ", contextIndex = " << contextIndex << ", validIndex = " << validIndex )
        ContextData& contextEntry = *mContextData[contextIndex];
        // now make reservation for enough memory, copy not needed
        reserve( contextIndex, mSize, false );

        // and then copy the data within the same context

        if( mSize > 0 )
        {
            fetch( contextEntry, otherEntry );
        }

        contextEntry.valid = true;
    }

*/
    COMMON_THROWEXCEPTION( "not available yet" )

    LAMA_LOG_DEBUG( logger, *this << " has now been assigned all valid data" )
    // copy of data in other contexts will be done only on demand later
    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::assign( const LAMAArray<ValueType>& other, ContextPtr context )
{
    LAMA_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    COMMON_THROWEXCEPTION( "not available yet" )

/*
    COMMON_ASSERT( this != &other, "self assign for LAMAArray forbidden" )

    wait();
    mSize = other.mSize;
    other.wait();

    // this array should have no read/write access

    for( size_t i = 0; i < mContextData.size(); i++ )
    {
        ContextData& entry = *mContextData[i];

        if( entry.lock[ContextData::Read] || entry.lock[ContextData::Write] )
        {
            COMMON_THROWEXCEPTION( *this << ": cannot copy/assign to locked array" )
        }

        // make all entries invalid because we will overwrite them

        if( entry.valid )
        {
            entry.valid = false;
            LAMA_LOG_DEBUG( logger, *this << ": invalidated at index = " << i << " for " << entry.context )
        }
    }

    //make other available at context
    other.prefetch( context );
    other.wait();

    size_t nOtherContexts = other.mContextData.size();

    bool copyDone = false;

    for( size_t i = 0; i < nOtherContexts; i++ )
    {
        LAMA_LOG_TRACE( logger, "Other Context index = " << i )
        const ContextData& otherEntry = *other.mContextData[i];

        if( otherEntry.context != context )
        {
            continue;
        }

        size_t contextIndex;
        size_t validIndex;
        getAccess( contextIndex, validIndex, otherEntry.context, ContextData::Write );
        LAMA_LOG_TRACE( logger,
                        *this << ": write access for " << *otherEntry.context << ", contextIndex = " << contextIndex << ", validIndex = " << validIndex )
        ContextData& contextEntry = *mContextData[contextIndex];
        // now make reservation for enough memory, copy not needed
        reserve( contextIndex, mSize, false );

        // and then copy the data within the same context

        if( mSize > 0 )
        {
            fetch( contextEntry, otherEntry );
        }

        contextEntry.valid = true;
        copyDone = true;
        break;
    }

    COMMON_ASSERT( copyDone, "assignment failed" )

*/

    LAMA_LOG_DEBUG( logger, *this << " has now been assigned at " << *context )
    // copy of data in other contexts will be done only on demand later
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

template<typename ValueType>
void LAMAArray<ValueType>::prefetch( ContextPtr context ) const
{
    COMMON_THROWEXCEPTION( "not available yet" )
/*
    // public routine -> always check
    COMMON_ASSERT( context, "NULL pointer for context" )
    size_t contextIndex;
    size_t validIndex;
    LAMA_LOG_DEBUG( logger, *this << ": prefetch on " << *context )
    getAccess( contextIndex, validIndex, context, ContextData::Read );
    LAMA_LOG_TRACE( logger,
                    "prefetch on " << *context << ": contextIndex = " << contextIndex << ", validIndex = " << validIndex )
    ContextData& contextEntry = *mContextData[contextIndex];

    if( contextEntry.valid || mSize == 0 )
    {
        return;
    }

    wait();
    COMMON_ASSERT( validIndex != nContextIndex, "no valid context for " << *this )
    const ContextData& validEntry = *mContextData[validIndex];
    reserve( contextIndex, mSize, false ); //  take care for sufficient memory
    mSyncToken.reset( fetchAsync( contextEntry, validEntry ) );
    // mSyncToken->wait();  // To be deleted
    contextEntry.valid = true;
*/

}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
bool LAMAArray<ValueType>::isAvailableAt( ContextPtr context ) const
{
    COMMON_ASSERT( context, "NULL pointer for context" )

    LAMA_LOG_DEBUG( logger, *this << ": check availability on " << *context )

    ContextDataRef index = mContextManager.getContextData( context, ContextData::Read );

    ContextData& data = mContextManager[ index ];

    return data.valid;
}

template<typename ValueType>
ContextPtr LAMAArray<ValueType>::getValidContext( const Context::ContextType preferredType /*= Context::Host*/) const
{
    ContextPtr validContext;

/*
    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = *mContextData[i];

        if( entry.valid )
        {
            validContext = entry.context;
        }
        else
        {
            COMMON_ASSERT( 0 == entry.lock[ContextData::Read], "read access on non valid location" )
            COMMON_ASSERT( 0 == entry.lock[ContextData::Write], "write access on non valid location" )
            continue;
        }

        if( preferredType == validContext->getType() )
        {
            break;
        }
    }
*/

    return validContext;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( ContextPtr context, const IndexType capacity )
{
    COMMON_THROWEXCEPTION( "not supported yet" )

    // ToDo: lock

/*
    size_t contextIndex = nContextIndex;

    // search for an available entry

    for( size_t i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = *mContextData[i];

        if( *entry.context == *context )
        {
            contextIndex = i;
        }
        else if( entry.lock[ContextData::Write] )
        {
            COMMON_THROWEXCEPTION( "no further access on write locked array " << *this )
        }
    }

    // if context is used first time, make new entry for context data

    if( contextIndex == nContextIndex )
    {
        contextIndex = mContextData.size();
        mContextData.push_back( new ContextData( context ) );
        LAMA_LOG_DEBUG( logger, "new context data entry for " << *context << ", index = " << contextIndex )

        // pointer = ..., size = , ... allocated, valid, lock ...
    }

    // Error if entry is locked

    // now resize

    bool copyFlag = mContextData[contextIndex]->valid;

    reserve( contextIndex, capacity, copyFlag );
*/

}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::fetch( ContextData& target, const ContextData& source ) const
{
/*
    LAMA_LOG_INFO( logger,
                   *this << ": fetch (size = " << mSize << ") from " << *source.context << " to " << *target.context )
    COMMON_ASSERT( source.valid, "fetch from invalid source" )
    COMMON_ASSERT( !target.valid, "fetch to valid target" )
    COMMON_ASSERT( mSize, "size = 0, no fetch needed" )
    COMMON_ASSERT( target.size >= mSize * sizeof(ValueType),
                       *this << ": fetch has insufficient capacity on target context " << target.context )
    size_t transferSize = mSize * sizeof(ValueType);

    if( target.context->getType() == Context::Host && source.context->getType() != Context::Host )
    {
        ValueType* target_pointer = static_cast<ValueType*>( target.pointer );

        // to ensure first touch
#pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

        for( IndexType i = 0; i < mSize; ++i )
        {
            target_pointer[i] = 0;
        }
    }

    if( source.context->getType() == Context::Host && target.context->getType() == Context::Host )
    {
        ValueType* target_pointer = static_cast<ValueType*>( target.pointer );

        ValueType* source_pointer = static_cast<ValueType*>( source.pointer );

        // to preserve first touch
#pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

        for( IndexType i = 0; i < mSize; ++i )
        {
            target_pointer[i] = source_pointer[i];
        }
    }
    else if( source.context->canUseData( *target.context ) )
    {
        LAMA_LOG_DEBUG( logger,
                        "same use context transfer to " << *target.context << " from " << *source.context << ", size = " << transferSize )
        source.context->memcpy( target.pointer, source.pointer, transferSize );
    }
    else if( target.context->cancpy( target, source ) )
    {
        LAMA_LOG_DEBUG( logger,
                        "transfer to " << *target.context << " from " << *source.context << ", size = " << transferSize )
        target.context->memcpy( target, source, transferSize );
    }
    else if( source.context->cancpy( target, source ) )
    {
        LAMA_LOG_DEBUG( logger,
                        "transfer from " << *source.context << " to " << *target.context << ", size = " << transferSize )
        source.context->memcpy( target, source, transferSize );
    }
    else
    {
        COMMON_ASSERT( target.context->getType() != Context::Host, "memcpyToThis for host unsupported" )
        COMMON_ASSERT( source.context->getType() != Context::Host, "memcpyFromThis for host unsupported" )
        LAMA_LOG_INFO( logger,
                       "transfer " << *target.context << " <- " << *source.context << " not supported, try via host" )
        ContextData& host = *mContextData[0];
        COMMON_ASSERT( host.context->getType() == Context::Host, "context 0 is not host" )
        fetch( host, source );
        host.valid = true; // must be valid now, otherwise next fetch fails
        fetch( target, host );
    }
*/
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
SyncToken* LAMAArray<ValueType>::fetchAsync( ContextData& target, const ContextData& source ) const
{
/*
    LAMA_LOG_INFO( logger,
                   *this << ": async fetch (size = " << mSize << ") from " << *source.context << " to " << *target.context )
    COMMON_ASSERT( source.valid, "async fetch from invalid source" )
    COMMON_ASSERT( !target.valid, "async fetch to valid target" )
    COMMON_ASSERT( mSize, "size = 0, no fetch needed" )
    COMMON_ASSERT( target.size >= mSize * sizeof(ValueType),
                   *this << ": fetch has insufficient capacity on target context " << target.context )

    size_t transferSize = mSize * sizeof(ValueType);

    if( source.context->canUseData( *target.context ) )
    {
        LAMA_LOG_DEBUG( logger,
                        "same use context async transfer to " << *target.context << " from " << *source.context << ", size = " << transferSize )
        return source.context->memcpyAsync( target.pointer, source.pointer, transferSize );
    }
    else if( target.context->cancpy( target, source ) )
    {
        LAMA_LOG_INFO( logger,
                       "async transfer from " << *source.context << " to " << *target.context << ", size = " << transferSize )
        return target.context->memcpyAsync( target, source, transferSize );
    }
    else if( source.context->cancpy( target, source ) )
    {
        LAMA_LOG_INFO( logger,
                       "async transfer from " << *source.context << " to " << *target.context << ", size = " << transferSize )
        return source.context->memcpyAsync( target, source, transferSize );
    }
    else
    {
        COMMON_THROWEXCEPTION(
                        "no async memory transfer from " << *source.context << " to " << *target.context << " supported" );
    }
*/

    return NULL;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::wait() const
{
/*
    LAMA_LOG_TRACE( logger, "wait" )

    if( 0 != mSyncToken.get() )
    {
        LAMA_LOG_DEBUG( logger, "Waiting for SyncToken: " << *mSyncToken )

        mSyncToken.reset(); // waits for transfer and frees resources
    }
*/

}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear()
{
    wait();

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
ValueType* LAMAArray<ValueType>::get( ContextDataRef index )
{
    return static_cast<ValueType*>( mContextManager[index].pointer );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const ValueType* LAMAArray<ValueType>::get( ContextDataRef index ) const
{
    return static_cast<const ValueType*>( mContextManager[index].pointer );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear( const ContextDataRef index )
{
    ContextData& data = mContextManager[index];

    COMMON_ASSERT( data.locked( ContextData::Write ), "clear illegal here " << data )

    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::resize( ContextDataRef index, const IndexType size )
{
    ContextData& entry = mContextManager[index];

    COMMON_ASSERT( entry.locked( ContextData::Write ), "resize illegal here " << entry )

    size_t allocSize = size * sizeof( ValueType );
    size_t validSize = mSize * sizeof( ValueType );

    LAMA_LOG_INFO( logger, *this << ": resize, needed = " << allocSize << " bytes, used = " 
                         << validSize << " bytes, capacity = " << entry.capacity() << " bytes" )

    entry.reserve( allocSize, validSize );

    // capacity is now sufficient for size elements

    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( ContextDataRef index, const IndexType size ) const
{
    COMMON_THROWEXCEPTION( "not available yet" )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType LAMAArray<ValueType>::capacity( ContextDataRef index ) const
{
    ContextData& entry = mContextManager[index];

    return static_cast<IndexType>( entry.size / sizeof( ValueType ) );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ContextDataRef LAMAArray<ValueType>::acquireReadAccess( ContextPtr context ) const
{
    common::Thread::ScopedLock lock( mAccessMutex );

    size_t allocSize = mSize * sizeof( ValueType );
    size_t validSize = allocSize;                   // read access needs valid data in any case

    return mContextManager.acquireAccess( context, ContextData::Read, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseReadAccess( ContextDataRef index ) const
{
    common::Thread::ScopedLock lock( mAccessMutex );

    mContextManager[index].releaseLock( ContextData::Read );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ContextDataRef LAMAArray<ValueType>::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    common::Thread::ScopedLock lock( mAccessMutex );

    size_t allocSize = mSize * sizeof( ValueType );
    size_t validSize = keepFlag ? allocSize : 0 ;

    return mContextManager.acquireAccess( context, ContextData::Write, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseWriteAccess( ContextDataRef index )
{
    common::Thread::ScopedLock lock( mAccessMutex );

    ContextData& data = mContextManager[index];

    data.releaseLock( ContextData::Write );

    LAMA_LOG_INFO( logger, "releaseWriteAccess: " << data );
}

/* ---------------------------------------------------------------------------------*/

/*
template<typename ValueType>
void LAMAArray<ValueType>::copy( const int toIndex, const int fromIndex ) const
{
    COMMON_THROWEXCEPTION( "copy from " << fromIndex << " to " << toIndex << " not available yet " )
}
*/

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "LAMAArray<";
    stream << Scalar::getType<ValueType>();
    stream << ">(" << mSize;

/*
    for( size_t i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = *mContextData[i];
        // unsigned char for locks must be casted for output stream
        stream << ", " << *entry.context << ": c=" << entry.size / sizeof(ValueType) << ", v=" << entry.valid << ", #r="
                        << static_cast<unsigned short>( entry.lock[ContextData::Read] ) << ", #w="
                        << static_cast<unsigned short>( entry.lock[ContextData::Write] );
    }

*/

    stream << ")";
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
