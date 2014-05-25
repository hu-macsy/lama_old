/**
 * @file LAMAArray.cpp
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
 * @brief Implementation of methods for class LAMAArray.
 * @author Thomas Brandes
 * @date 14.03.2011
 * @since 1.0.0
 */

// hpp
#include <lama/LAMAArray.hpp>

// others
#include <lama/exception/Exception.hpp>
#include <lama/exception/LAMAAssert.hpp>

#include <lama/ContextFactory.hpp>
#include <lama/SyncToken.hpp>

// boost
#include <boost/thread/recursive_mutex.hpp>
#include <boost/preprocessor.hpp>

#include <memory>
#include <limits>

namespace lama
{

typedef Context::ContextData ContextData;
typedef Context::ContextData::AccessKind AccessKind;

/* -------------------------------------------------------------------------- */

_LAMAArray* _LAMAArray::create( const Scalar::ScalarType valueType )
{
    switch ( valueType )
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
        LAMA_THROWEXCEPTION( "Unsupported ValueType " << valueType )
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
    : _LAMAArray( 0 ), mSyncToken( 0 )
{
    setHostContext();
    LAMA_LOG_DEBUG( logger, "created new context array: " << *this )
}

template<typename ValueType>
void LAMAArray<ValueType>::setHostContext()
{
    LAMA_ASSERT_DEBUG( mContextData.size() == 0, "setHostContext only on new LAMA array" )
    // good alternative: ContextFactory::maxContexts()
    mContextData.reserve( LAMA_MAX_CONTEXTS );
    // first location is always on Host Context.
    ContextData* hostContext = new ContextData( ContextFactory::getContext( Context::Host ) );
    hostContext->valid = true; // an empty array is also valid
    mContextData.push_back( hostContext );
    LAMA_LOG_TRACE( logger, "Host context set for LAMA array: " << *this )
}

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n )
    : _LAMAArray( n ), mSyncToken( 0 )
{
    setHostContext();

    if ( n > 0 )
    {
        ContextData& host = *mContextData[0];
        host.allocate( mSize * sizeof(ValueType) );
    }

    LAMA_LOG_DEBUG( logger, "constructed: " << *this )
}

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n, const ValueType& value )
    : _LAMAArray( n ), mSyncToken( 0 )
{
    setHostContext();

    if ( n <= 0 )
    {
        LAMA_LOG_INFO( logger, "Construtor LAMAArray, size = 0, value = " << value )
        return;
    }

    ContextData& host = *mContextData[0];
    host.allocate( mSize * sizeof(ValueType) );
    LAMA_LOG_DEBUG( logger, "constructed: " << *this )

    ValueType* host_pointer = static_cast<ValueType*>( host.pointer );

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( int i = 0; i < mSize; ++i )
    {
        host_pointer[i] = value;
    }

    host.valid = true;
    LAMA_LOG_DEBUG( logger, "constructed: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const LAMAArray<ValueType>& other )
    : _LAMAArray( other.mSize ), mSyncToken( 0 )
{
    // find a valid location of the other array in any context
    setHostContext(); // make default initialization
    LAMA_LOG_DEBUG( logger, other << " will be copied to " << *this )
    operator=( other );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::~LAMAArray()
{
    wait();

    LAMA_LOG_DEBUG( logger, "~LAMAArray = " << *this )

    // explicit free for all context data needed

    bool locked = false;

    // check for existing read / write lock

    if ( LAMA_LOG_WARN_ON(logger) )
    {
        for ( size_t i = 0; i < mContextData.size(); i++ )
        {
            if ( mContextData[i]->lock[ContextData::Read] || mContextData[i]->lock[ContextData::Write] )
            {
                locked = true;
                break;
            }
        }

    }

    // give a warning, never throw an exception in destructor ( might crash )

    if ( locked )
    {
        LAMA_LOG_WARN( logger, "Destructor on read/write locked array: " << *this )
    }

    {
        for ( size_t i = 0; i < mContextData.size(); i++ )
        {
            delete mContextData[i];
        }
    }
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

    LAMA_ASSERT_ERROR( this != &other, "self assign for LAMAArray forbidden" )

    wait();
    mSize = other.mSize;
    other.wait();

    // this array should have no read/write access

    for ( size_t i = 0; i < mContextData.size(); i++ )
    {
        ContextData& entry = *mContextData[i];

        if ( entry.lock[ContextData::Read] || entry.lock[ContextData::Write] )
        {
            LAMA_THROWEXCEPTION( *this << ": cannot copy/assign to locked array" )
        }

        // make all entries invalid because we will overwrite them

        if ( entry.valid )
        {
            entry.valid = false;
            LAMA_LOG_DEBUG( logger, *this << ": invalidated at index = " << i << " for " << entry.context )
        }
    }

    // each valid data of the other array will be copied into the same context for this array
    size_t nOtherContexts = other.mContextData.size();

    for ( size_t i = 0; i < nOtherContexts; i++ )
    {
        LAMA_LOG_TRACE( logger, "Other Context index = " << i )
        const ContextData& otherEntry = *other.mContextData[i];

        if ( !otherEntry.valid )
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

        if ( mSize > 0 )
        {
            fetch( contextEntry, otherEntry );
        }

        contextEntry.valid = true;
    }

    LAMA_LOG_DEBUG( logger, *this << " has now been assigned all valid data" )
    // copy of data in other contexts will be done only on demand later
    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::assign( const LAMAArray<ValueType>& other, ContextPtr context )
{
    LAMA_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    LAMA_ASSERT_ERROR( this != &other, "self assign for LAMAArray forbidden" )

    wait();
    mSize = other.mSize;
    other.wait();

    // this array should have no read/write access

    for ( size_t i = 0; i < mContextData.size(); i++ )
    {
        ContextData& entry = *mContextData[i];

        if ( entry.lock[ContextData::Read] || entry.lock[ContextData::Write] )
        {
            LAMA_THROWEXCEPTION( *this << ": cannot copy/assign to locked array" )
        }

        // make all entries invalid because we will overwrite them

        if ( entry.valid )
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

    for ( size_t i = 0; i < nOtherContexts; i++ )
    {
        LAMA_LOG_TRACE( logger, "Other Context index = " << i )
        const ContextData& otherEntry = *other.mContextData[i];

        if ( otherEntry.context != context )
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

        if ( mSize > 0 )
        {
            fetch( contextEntry, otherEntry );
        }

        contextEntry.valid = true;
        copyDone = true;
        break;
    }

    LAMA_ASSERT_ERROR( copyDone, "assignment failed" )

    LAMA_LOG_DEBUG( logger, *this << " has now been assigned at " << *context )
    // copy of data in other contexts will be done only on demand later
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::swap( LAMAArray<ValueType>& other )
{
    LAMA_LOG_DEBUG( logger, *this << ": swap with other = " << other )

    // we cannot swap if there is any access for any array

    for ( size_t i = 0; i < mContextData.size(); i++ )
    {
        LAMA_ASSERT_ERROR( 0 == mContextData[i]->lock[ContextData::Read],
                           *this << ": cannot be swapped, " << " read lock on " << *mContextData[i]->context )
        LAMA_ASSERT_ERROR( 0 == mContextData[i]->lock[ContextData::Write],
                           *this << ": cannot be swapped, " << " write lock on " << *mContextData[i]->context )
    }

    for ( size_t i = 0; i < other.mContextData.size(); i++ )
    {
        LAMA_ASSERT_ERROR( 0 == other.mContextData[i]->lock[ContextData::Read],
                           other << ": cannot be swapped, " << " read lock on " << *other.mContextData[i]->context )
        LAMA_ASSERT_ERROR( 0 == other.mContextData[i]->lock[ContextData::Write],
                           other << ": cannot be swapped, " << " write lock on " << *other.mContextData[i]->context )
    }

    std::swap( mSize, other.mSize );
    std::swap( mContextData, other.mContextData );
    std::swap( mSyncToken, other.mSyncToken );
    LAMA_LOG_DEBUG( logger, *this << ": has been swapped with other = " << other )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::prefetch( ContextPtr context ) const
{
    // public routine -> always check
    LAMA_ASSERT_ERROR( context, "NULL pointer for context" )
    size_t contextIndex;
    size_t validIndex;
    LAMA_LOG_DEBUG( logger, *this << ": prefetch on " << *context )
    getAccess( contextIndex, validIndex, context, ContextData::Read );
    LAMA_LOG_TRACE( logger,
                    "prefetch on " << *context << ": contextIndex = " << contextIndex << ", validIndex = " << validIndex )
    ContextData& contextEntry = *mContextData[contextIndex];

    if ( contextEntry.valid || mSize == 0 )
    {
        return;
    }

    wait();
    LAMA_ASSERT_DEBUG( validIndex != nContextIndex, "no valid context for " << *this )
    const ContextData& validEntry = *mContextData[validIndex];
    reserve( contextIndex, mSize, false ); //  take care for sufficient memory
    mSyncToken.reset( fetchAsync( contextEntry, validEntry ) );
    // mSyncToken->wait();  // To be deleted
    contextEntry.valid = true;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
bool LAMAArray<ValueType>::isAvailableAt( ContextPtr context ) const
{
    LAMA_ASSERT_ERROR( context, "NULL pointer for context" )
    size_t contextIndex;
    size_t validIndex;
    LAMA_LOG_DEBUG( logger, *this << ": check availability on " << *context )
    getAccess( contextIndex, validIndex, context, ContextData::Read );
    return mContextData[contextIndex]->valid;
}

template<typename ValueType>
ContextPtr LAMAArray<ValueType>::getValidContext( const Context::ContextType preferredType /*= Context::Host*/) const
{
    ContextPtr validContext;
    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = *mContextData[i];

        if ( entry.valid )
        {
            validContext = entry.context;
        }
        else
        {
            LAMA_ASSERT_DEBUG( 0 == entry.lock[ContextData::Read], "read access on non valid location" )
            LAMA_ASSERT_DEBUG( 0 == entry.lock[ContextData::Write], "write access on non valid location" )
            continue;
        }

        if ( preferredType == validContext->getType() )
        {
            break;
        }
    }
    return validContext;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::fetch( ContextData& target, const ContextData& source ) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": fetch (size = " << mSize << ") from " << *source.context << " to " << *target.context )
    LAMA_ASSERT_DEBUG( source.valid, "fetch from invalid source" )
    LAMA_ASSERT_DEBUG( !target.valid, "fetch to valid target" )
    LAMA_ASSERT_DEBUG( mSize, "size = 0, no fetch needed" )
    LAMA_ASSERT_DEBUG( target.size >= mSize*sizeof(ValueType),
                       *this << ": fetch has insufficient capacity on target context " << target.context )
    size_t transferSize = mSize * sizeof(ValueType);

    if ( target.context->getType() == Context::Host && source.context->getType() != Context::Host )
    {
        ValueType* target_pointer = static_cast<ValueType*>( target.pointer );

        // to ensure first touch
        #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
        for ( IndexType i = 0; i < mSize; ++i )
        {
            target_pointer[i] = 0;
        }
    }

    if ( source.context->getType() == Context::Host && target.context->getType() == Context::Host )
    {
        ValueType* target_pointer = static_cast<ValueType*>( target.pointer );

        ValueType* source_pointer = static_cast<ValueType*>( source.pointer );

        // to preserve first touch
        #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
        for ( IndexType i = 0; i < mSize; ++i )
        {
            target_pointer[i] = source_pointer[i];
        }
    }
    else if ( source.context->canUseData( *target.context ) )
    {
        LAMA_LOG_DEBUG( logger,
                        "same use context transfer to " << *target.context << " from " << *source.context << ", size = " << transferSize )
        source.context->memcpy( target.pointer, source.pointer, transferSize );
    }
    else if ( target.context->cancpy( target, source ) )
    {
        LAMA_LOG_DEBUG( logger,
                        "transfer to " << *target.context << " from " << *source.context << ", size = " << transferSize )
        target.context->memcpy( target, source, transferSize );
    }
    else if ( source.context->cancpy( target, source ) )
    {
        LAMA_LOG_DEBUG( logger,
                        "transfer from " << *source.context << " to " << *target.context << ", size = " << transferSize )
        source.context->memcpy( target, source, transferSize );
    }
    else
    {
        LAMA_ASSERT( target.context->getType() != Context::Host, "memcpyToThis for host unsupported" )
        LAMA_ASSERT( source.context->getType() != Context::Host, "memcpyFromThis for host unsupported" )
        LAMA_LOG_INFO( logger,
                       "transfer " << *target.context << " <- " << *source.context << " not supported, try via host" )
        ContextData& host = *mContextData[0];
        LAMA_ASSERT( host.context->getType() == Context::Host, "context 0 is not host" )
        fetch( host, source );
        host.valid = true; // must be valid now, otherwise next fetch fails
        fetch( target, host );
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
SyncToken* LAMAArray<ValueType>::fetchAsync( ContextData& target, const ContextData& source ) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": async fetch (size = " << mSize << ") from " << *source.context << " to " << *target.context )
    LAMA_ASSERT_DEBUG( source.valid, "async fetch from invalid source" )
    LAMA_ASSERT_DEBUG( !target.valid, "async fetch to valid target" )
    LAMA_ASSERT_DEBUG( mSize, "size = 0, no fetch needed" )
    LAMA_ASSERT_DEBUG( target.size >= mSize*sizeof(ValueType),
                       *this << ": fetch has insufficient capacity on target context " << target.context )

    size_t transferSize = mSize * sizeof(ValueType);

    if ( source.context->canUseData( *target.context ) )
    {
        LAMA_LOG_DEBUG( logger,
                        "same use context async transfer to " << *target.context << " from " << *source.context << ", size = " << transferSize )
        return source.context->memcpyAsync( target.pointer, source.pointer, transferSize );
    }
    else if ( target.context->cancpy( target, source ) )
    {
        LAMA_LOG_INFO( logger,
                       "async transfer from " << *source.context << " to " << *target.context << ", size = " << transferSize )
        return target.context->memcpyAsync( target, source, transferSize );
    }
    else if ( source.context->cancpy( target, source ) )
    {
        LAMA_LOG_INFO( logger,
                       "async transfer from " << *source.context << " to " << *target.context << ", size = " << transferSize )
        return source.context->memcpyAsync( target, source, transferSize );
    }
    else
    {
        LAMA_THROWEXCEPTION(
            "no async memory transfer from " << *source.context << " to " << *target.context << " supported" );
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::wait() const
{
    LAMA_LOG_TRACE( logger, "wait" )

    if ( 0 != mSyncToken.get() )
    {
        LAMA_LOG_DEBUG( logger, "Waiting for SyncToken: " << *mSyncToken )

        mSyncToken.reset();  // waits for transfer and frees resources
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear()
{
    LAMA_LOG_DEBUG( logger, *this << ": clear" )
    wait();

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        LAMA_ASSERT( 0 == mContextData[i]->lock[ContextData::Read], "Tried to clear a locked LAMAArray " << *this )
        LAMA_ASSERT( 0 == mContextData[i]->lock[ContextData::Write], "Tried to clear a locked LAMAArray " << *this )
    }

    mSize = 0;

    LAMA_LOG_TRACE( logger, *this << " cleared" )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::purge()
{
    LAMA_LOG_DEBUG( logger, *this << " will be purged" )
    wait();

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        ContextData& entry = *mContextData[i];
        LAMA_ASSERT( entry.lock[ContextData::Read] == 0, "Tried to purge a locked LAMAArray " << *this )
        LAMA_ASSERT( entry.lock[ContextData::Write] == 0, "Tried to purge a locked LAMAArray " << *this )
        entry.free();
    }

    mSize = 0;
    LAMA_LOG_TRACE( logger, *this << " purged" )
    // we still keep the context location entries
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType* LAMAArray<ValueType>::get( size_t index )
{
    LAMA_ASSERT_DEBUG( index < mContextData.size(), "Invalid context index = " << index )
    ContextData& entry = *mContextData[index];
    LAMA_ASSERT_DEBUG( entry.valid, "Tried to get pointer for invalid data on context " << *entry.context )
    //if we have a write access it is not possible to have a transfer running
    LAMA_ASSERT_DEBUG( 0 == mSyncToken.get(), "Although a write access exists we have a running transfer." )
    LAMA_LOG_TRACE( logger, *this << ": get at " << *entry.context << ", ptr = " << entry.pointer )
    return static_cast<ValueType*>( entry.pointer );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const ValueType* LAMAArray<ValueType>::get( size_t index ) const
{
    LAMA_ASSERT_DEBUG( index < mContextData.size(), "Invalid context index = " << index )
    const ContextData& entry = *mContextData[index];
    LAMA_ASSERT_DEBUG( entry.valid, "Tried to get pointer for invalid data on context " << *entry.context )
    LAMA_LOG_TRACE( logger, *this << ": const get at " << *entry.context << ", ptr = " << entry.pointer )
    return static_cast<const ValueType*>( entry.pointer );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear( const size_t index )
{
    LAMA_ASSERT_DEBUG( index < mContextData.size(), "Invalid context index = " << index )
    LAMA_ASSERT_ERROR( 1 == mContextData[index]->lock[ContextData::Write],
                       *this << ": tried to clear on " << *mContextData[index]->context << " without WriteAccess" )
    mSize = 0; // just reset the size, no free of any memory
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::resize( const size_t index, const IndexType size )
{
    LAMA_ASSERT_DEBUG( index < mContextData.size(), "Invalid context index = " << index )
    ContextData& entry = *mContextData[index];
    LAMA_ASSERT_ERROR(
        1 == entry.lock[ContextData::Write],
        *this << ": tried to resize on " << *entry.context << ", new size = " << size << " without WriteAccess" )

    if ( size * sizeof(ValueType) > entry.size )
    {
        reserve( index, size, true ); // copies old data
    }

    // changing the size does not imply any initialization
    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( const size_t index, const IndexType capacity, const bool copy ) const
{
    LAMA_ASSERT_DEBUG( index < mContextData.size(), "Invalid context index = " << index )
    // Note: reserve is also possible on non-valid contexts and without write locks
    ContextData& entry = *mContextData[index];
    LAMA_ASSERT_ERROR( 0 == mSyncToken.get(), "Although a write access exists we have a running transfer." )

    if ( capacity * sizeof(ValueType) > entry.size )
    {
        LAMA_LOG_TRACE( logger,
                        "growing capacity of " << *this << " from " << entry.size/sizeof(ValueType) << " to " << capacity << " at " << *entry.context )
        LAMA_ASSERT_DEBUG( !copy || entry.valid, "reserve with copy on non-valid context" )

        if ( entry.pointer )
        {
            IndexType copySize = 0;

            if ( copy )
            {
                copySize = mSize;
            }

            entry.realloc( capacity * sizeof(ValueType), copySize * sizeof(ValueType) );
        }
        else
        {
            entry.allocate( capacity * sizeof(ValueType) );
        }
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType LAMAArray<ValueType>::capacity( const size_t index ) const
{
    LAMA_ASSERT_DEBUG( index < mContextData.size(), "Invalid context index = " << index )

    // Note: reserve is also possible on non-valid contexts and without write locks

    ContextData& entry = *mContextData[index];

    return static_cast<IndexType>( entry.size / sizeof(ValueType) );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
int LAMAArray<ValueType>::acquireReadAccess( ContextPtr context ) const
{
    boost::recursive_mutex::scoped_lock scoped_lock( access_mutex );

    LAMA_ASSERT_DEBUG( context, "NULL pointer for context" )
    size_t contextIndex; // for entry that belongs to context
    size_t validIndex;
    LAMA_LOG_TRACE( logger, "acquire read access on " << *context )
    getAccess( contextIndex, validIndex, context, ContextData::Read );
    LAMA_LOG_TRACE( logger,
                    "read access on " << *context << ": contextIndex = " << contextIndex << ", validIndex = " << validIndex )
    ContextData& contextEntry = *mContextData[contextIndex];
    wait();

    // fetch only if size > 0, there might be no valid location for mSize == 0

    if ( !contextEntry.valid && ( mSize > 0 ) )
    {
        LAMA_ASSERT_DEBUG( validIndex != nContextIndex, "no valid context for " << *this )
        LAMA_LOG_TRACE( logger, "data not valid at " << *contextEntry.context )
        const ContextData& validEntry = *mContextData[validIndex];
        // make sure that we have enough memory on the target context
        // old data is invalid so it must not be saved.
        reserve( contextIndex, mSize, false );
        fetch( contextEntry, validEntry );
    }
    else if ( mSize > 0 )
    {
        LAMA_LOG_TRACE( logger, "data already valid for " << *contextEntry.context )
    }

    contextEntry.valid = true;
    ++contextEntry.lock[ContextData::Read];
    LAMA_LOG_TRACE( logger,
                    "ready acquire read access on " << *contextEntry.context << ", #read locks = " << static_cast<unsigned short>( contextEntry.lock[ContextData::Read] ) )
    return static_cast<int>( contextIndex );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseReadAccess( const size_t index ) const
{
    boost::recursive_mutex::scoped_lock scoped_lock( access_mutex );

    LAMA_ASSERT_DEBUG( index < mContextData.size(), "index out of range" )
    ContextData& entry = *mContextData[index];
    LAMA_LOG_TRACE( logger,
                    "release read access on " << *entry.context << ", #read locks = " << static_cast<unsigned short>( entry.lock[ContextData::Read] ) )
    LAMA_ASSERT( entry.lock[ContextData::Read] > 0,
                 "Tried to release a non existing ReadAccess on " << *entry.context )
    --entry.lock[ContextData::Read];
    LAMA_LOG_TRACE( logger,
                    "ready release read access on " << *entry.context << ", #read locks = " << static_cast<unsigned short>( entry.lock[ContextData::Read] ) )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::getAccess(
    size_t& contextIndex,
    size_t& validIndex,
    ContextPtr context,
    AccessKind kind ) const
{
    contextIndex = nContextIndex;
    validIndex = nContextIndex;
    LAMA_LOG_TRACE( logger, "check access for " << *context )

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = *mContextData[i];

        if ( context->canUseData( *entry.context ) )
        {
            LAMA_LOG_TRACE( logger, "can use context at index " << i )
            contextIndex = i;
        }

        if ( !entry.valid )
        {
            LAMA_ASSERT_DEBUG( 0 == entry.lock[ContextData::Read], "read access on non valid location" )
            LAMA_ASSERT_DEBUG( 0 == entry.lock[ContextData::Write], "write access on non valid location" )
            continue;
        }

        //  so we have found a context with valid data

        if ( validIndex == nContextIndex )
        {
            validIndex = i;
        }

        if ( entry.lock[ContextData::Read] )
        {
            LAMA_ASSERT_DEBUG( 0 == entry.lock[ContextData::Write], "write and read access" )

            if ( kind == ContextData::Write )
            {
                LAMA_THROWEXCEPTION( "try to get write access on read locked array " << *this )
            }
        }
        else if ( entry.lock[ContextData::Write] )
        {
            LAMA_THROWEXCEPTION( "no further access on write locked array " << *this )
        }
    }

    // if context is used first time, make new entry for context data

    if ( contextIndex == nContextIndex )
    {
        contextIndex = mContextData.size();
        mContextData.push_back( new ContextData( context ) );
        LAMA_LOG_DEBUG( logger, "new context data entry for " << *context << ", index = " << contextIndex )
    }

    LAMA_LOG_TRACE( logger,
                    "get access on " << *context << ": contextIndex = " << contextIndex << ", validIndex = " << validIndex )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
int LAMAArray<ValueType>::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    boost::recursive_mutex::scoped_lock scoped_lock( access_mutex );

    LAMA_ASSERT_DEBUG( context, "NULL pointer for context" )
    LAMA_LOG_TRACE( logger, "acquire write access on " << *context )
    size_t contextIndex; // index that belongs to context
    size_t validIndex; // index for a valid entry
    getAccess( contextIndex, validIndex, context, ContextData::Write );
    LAMA_LOG_TRACE( logger,
                    "write access on " << *context << ": contextIndex = " << contextIndex << ", validIndex = " << validIndex )
    ContextData& contextEntry = *mContextData[contextIndex];
    wait();

    // fetch only if size > 0, there might be no valid location for mSize == 0

    if ( !contextEntry.valid && mSize > 0 )
    {
        LAMA_LOG_TRACE( logger, "data not valid at " << *contextEntry.context )
        // make sure that we have enough memory on the target context
        // old data is invalid so it must not be saved.
        reserve( contextIndex, mSize, 0 );

        if ( keepFlag )
        {
            LAMA_ASSERT_DEBUG( validIndex != nContextIndex, "no valid context for " << *this )
            const ContextData& validEntry = *mContextData[validIndex];
            fetch( contextEntry, validEntry );
        }
    }

    size_t noContexts = mContextData.size();
    LAMA_LOG_TRACE( logger, "invalidate for " << noContexts << " context locations" )

    for ( size_t i = 0; i < noContexts; ++i )
    {
        mContextData[i]->valid = false;
    }

    contextEntry.valid = true;
    contextEntry.lock[ContextData::Write] = true;
    LAMA_LOG_TRACE( logger, *this << ": ready acquire write access on " << *contextEntry.context )
    return static_cast<int>( contextIndex );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
int LAMAArray<ValueType>::acquireWriteAccess()
{
    boost::recursive_mutex::scoped_lock scoped_lock( access_mutex );
    size_t contextIndex = nContextIndex;

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = *mContextData[i];

        if ( entry.valid )
        {
            contextIndex = acquireWriteAccess( entry.context, true );
            break;
        }
    }

    return static_cast<int>( contextIndex );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseWriteAccess( const size_t index )
{
    boost::recursive_mutex::scoped_lock scoped_lock( access_mutex );

    LAMA_ASSERT_DEBUG( index < mContextData.size(), "index out of range" )
    ContextData& entry = *mContextData[index];
    LAMA_ASSERT( entry.lock[ContextData::Write] > 0,
                 "Tried to release a non existing WriteAccess on " << *entry.context )
    --entry.lock[ContextData::Write];
    LAMA_LOG_TRACE( logger, *this << ": released WriteAccess for " << *entry.context )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::copy( const int toIndex, const int fromIndex ) const
{
    LAMA_THROWEXCEPTION( "copy from " << fromIndex << " to " << toIndex << " not available yet " )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "LAMAArray<";
    stream << Scalar::getType<ValueType>();
    stream << ">(" << mSize;

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = *mContextData[i];
        // unsigned char for locks must be casted for output stream
        stream << ", " << *entry.context << ": c=" << entry.size / sizeof(ValueType) << ", v=" << entry.valid << ", #r="
               << static_cast<unsigned short>( entry.lock[ContextData::Read] ) << ", #w="
               << static_cast<unsigned short>( entry.lock[ContextData::Write] );
    }

    stream << ")";
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( ValueType* pointer, IndexType size )
    : LAMAArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if ( size != 0 && pointer == NULL )
    {
        LAMA_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

    ContextData& host = *mContextData[0];
    host.setRef( pointer, size * sizeof(ValueType) );

    mSize = size;
}

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( const ValueType* pointer, IndexType size )
    : LAMAArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if ( size != 0 && pointer == NULL )
    {
        LAMA_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

    ContextData& host = *mContextData[0];
    host.setRef( const_cast<ValueType*>( pointer ), size * sizeof(ValueType) );

    constFlag = true; // makes sure that we cannot have a WriteAccess

    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

#define LAMA_ARRAY_INSTANTIATE(z, I, _)                               \
template class LAMA_DLL_IMPORTEXPORT LAMAArray<ARRAY_TYPE##I> ;       \
template class LAMA_DLL_IMPORTEXPORT LAMAArrayRef<ARRAY_TYPE##I> ;
 
// template instantiation for the supported data types

BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_ARRAY_INSTANTIATE, _ )

#undef LAMA_ARRAY_INSTANTIATE


} // namespace LAMA
