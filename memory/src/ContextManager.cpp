/**
 * @file ContextManger.cpp
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
 * @brief Implementation of methods for class ContextManager.
 * @author Thomas Brandes
 * @date 14.03.2011
 */

// hpp

#include <memory/ContextManager.hpp>
#include <memory/Context.hpp>

#include <common/Exception.hpp>

namespace memory
{

typedef ContextData::AccessKind AccessKind;

/* ---------------------------------------------------------------------------------*/

LAMA_LOG_DEF_LOGGER( ContextManager::logger, "ContextManager" )

/* ---------------------------------------------------------------------------------*/

ContextManager::ContextManager() :

    mContextData(),
    mSyncToken()

{
    mLock[ContextData::Read] = 0;
    mLock[ContextData::Write] = 0;
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::wait()
{
    if ( 0 != mSyncToken.get() )
    {
        LAMA_LOG_DEBUG( logger, "Waiting for SyncToken: " << *mSyncToken )
        mSyncToken.reset(); // waits for transfer and frees resources
    }
}

/* ---------------------------------------------------------------------------------*/

bool ContextManager::locked() const
{
    return ( mLock[ContextData::Write] > 0 ) || ( mLock[ContextData::Read] > 0 );
}

/* ---------------------------------------------------------------------------------*/

bool ContextManager::locked( ContextData::AccessKind kind ) const
{
    return ( mLock[kind] > 0 );
}

/* ---------------------------------------------------------------------------------*/

bool ContextManager::isAccessContext( ContextPtr context )
{
    if ( !accessContext )
    {
        return false;
    }

    return context.get() == accessContext.get();
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::lockAccess( ContextData::AccessKind kind, ContextPtr context )
{
   //  common::Thread::ScopedLock lock( mAccessMutex );

    if ( ( mLock[ContextData::Read] == 0 ) && ( mLock[ContextData::Write] == 0 ) )
    {
        // first access, so we can set access context

        accessContext = context;
    } 
    else if ( !isAccessContext( context ) )
    {
        // some checks are required

        if ( kind == ContextData::Read )
        {
            COMMON_ASSERT_EQUAL( 0, mLock[ContextData::Write], "write lock, no read access possible" )
        }
        else if ( kind == ContextData::Write )
        {
            COMMON_ASSERT_EQUAL( 0, mLock[ContextData::Read], "no write access possible, is locked" )
            COMMON_ASSERT_EQUAL( 0, mLock[ContextData::Write], "no write access possible, is locked" )
        }

        // multiple reads will set access context to NULL

        accessContext.reset();
    }

    mLock[kind]++;
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::unlockAccess( ContextData::AccessKind kind )
{
   //  common::Thread::ScopedLock lock( mAccessMutex );

    COMMON_ASSERT_LE( 1, mLock[kind], "release access " << kind << ", never acquired" )

    mLock[kind]--;

    if ( ( mLock[ContextData::Write] == 0 ) && ( mLock[ContextData::Write] == 0 ) )
    {
        accessContext.reset();
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::releaseAccess( ContextDataIndex index, ContextData::AccessKind kind )
{
    // we should check that this is really the context data for which access was reserved
    unlockAccess( kind );
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::purge()
{
    // purge frees all data but keeps the ContextData entries
    wait();
    COMMON_ASSERT( !locked(), "purge on array with access" )

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        ContextData& entry = ( *this )[i];
        entry.free();
    }
}

/* ---------------------------------------------------------------------------------*/

ContextManager::~ContextManager()
{
    wait();
    // explicit free for all context data needed
    bool isLocked = false;

    // check for existing read / write lock

    if ( LAMA_LOG_WARN_ON( logger ) )
    {
        isLocked = locked();
    }

    // give a warning, never throw an exception in destructor ( might crash )

    if ( isLocked )
    {
        LAMA_LOG_WARN( logger, "Destructor on read/write locked array: " )
    }

    for ( size_t i = 0; i < mContextData.size(); i++ )
    {
        LAMA_LOG_INFO( logger, "~ContextManager, free " << mContextData[i] )
        mContextData[i].free();
    }

    mContextData.clear();
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextManager::findContextData( ContextPtr context ) const
{
    ContextDataIndex i;

    for ( i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = mContextData[i];

        if ( context->canUseData( *entry.context() ) )
        {
            break;
        }
    }

    return i;
}

/* ---------------------------------------------------------------------------------*/

bool ContextManager::isValid( ContextPtr context ) const
{
    ContextDataIndex index = findContextData( context );

    if ( index < mContextData.size() )
    {
        const ContextData& entry = mContextData[index];
        return entry.isValid();
    }
    else
    {
        return false;
    }
}

/* ---------------------------------------------------------------------------------*/

size_t ContextManager::capacity( ContextPtr context ) const
{
    ContextDataIndex index = findContextData( context );

    if ( index < mContextData.size() )
    {
        const ContextData& entry = mContextData[index];
        return entry.capacity();
    }
    else
    {
        return 0;
    }
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextManager::getContextData( ContextPtr context )
{
    bool found = false;
    size_t contextIndex = findContextData( context );

    // if context is used first time, make new entry for context data

    if ( contextIndex == mContextData.size() )
    {
        if ( contextIndex == 0 )
        {
            mContextData.reserve( LAMA_MAX_CONTEXTS );
        }

        mContextData.push_back( ContextData( context ) );
        LAMA_LOG_DEBUG( logger, "new context data entry for " << *context << ", index = " << contextIndex )
    }

    return contextIndex;
}

/* ---------------------------------------------------------------------------------*/

ContextData& ContextManager::operator[] ( ContextDataIndex index )
{
    COMMON_ASSERT( index < mContextData.size(), "index = " << index << " is illegal index, size = " << mContextData.size() )
    return mContextData[index];
}

/* ---------------------------------------------------------------------------------*/

ContextData& ContextManager::operator[] ( ContextPtr context )
{
    return operator[]( getContextData( context ) );
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextManager::findValidData() const
{
    ContextDataIndex index;

    for ( index = 0; index < mContextData.size(); ++index )
    {
        const ContextData& entry = mContextData[index];
        LAMA_LOG_INFO( logger, "check valid: " << entry )

        if ( entry.isValid() )
        {
            LAMA_LOG_INFO( logger, "found valid entry at index = " << index << ": " << entry )
            break;
        }
    }

    if ( index == mContextData.size() )
    {
        LAMA_LOG_INFO( logger, "no valid data found, index = " << index )
    }

    return index;
}

/* ---------------------------------------------------------------------------------*/

const ContextData& ContextManager::getValidData() const
{
    ContextDataIndex index = findValidData();

    if ( index < mContextData.size() )
    {
        return mContextData[index];
    }
    else
    {
        COMMON_THROWEXCEPTION( "no valid data found (zero size, uninitialized array)" )
    }
}

/* ---------------------------------------------------------------------------------*/

ContextPtr ContextManager::getValidContext( const Context::ContextType preferredType )
{
    ContextPtr result;

    for ( size_t index = 0; index < mContextData.size(); ++index )
    {
        const ContextData& entry = mContextData[index];

        if ( entry.isValid() )
        {
            if ( entry.context()->getType() == preferredType )
            {
                return entry.context();
            }
            else if ( result )
            {
                // do not overwrite first context found
            }
            else
            {
                result = entry.context();
            }
        }
    }

    return result;  // might be NULL
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextManager::acquireAccess( ContextPtr context, AccessKind kind,
        size_t allocSize, size_t validSize )
{
    COMMON_ASSERT( context, "NULL pointer for context" )
    lockAccess( kind, context );
    LAMA_LOG_DEBUG( logger, "acquire access on " << *context << ", kind = " << kind
                    << ", allocSize = " << allocSize << ", validSize = " << validSize )
    ContextDataIndex index = getContextData( context );
    ContextData& data = ( *this )[index];
    wait();

    // fetch only if size > 0, there might be no valid location for mSize == 0

    if ( !data.isValid() && allocSize > 0 )
    {
        LAMA_LOG_DEBUG( logger, "data not valid at " << *data.context() )
        // make sure that we have enough memory on the target context
        // old data is invalid so it must not be saved.
        data.reserve( allocSize, 0 );  // do not save any old values

        if ( validSize )
        {
            const ContextData& validEntry = getValidData();
            LAMA_LOG_INFO( logger, "valid data here: " << validEntry )
            fetch( data, validEntry, validSize );
        }
    }

    if ( kind == ContextData::Write )
    {
        invalidateAll();        // invalidate all entries
        data.setValid( true );  // for next access the data @ context is valid.
    }

    LAMA_LOG_DEBUG( logger, "acquired access :" << data );
    return index;
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::fetch( ContextData& target, const ContextData& source, size_t size )
{
    try
    {
        target.copyFrom( source, size );
    }
    catch ( common::Exception& ex )
    {
        LAMA_LOG_INFO( logger, target << " copy from " << source << " not supported" )
        // try it via host
        ContextPtr hostContext = Context::getContext( Context::Host );

        if ( target.context()->getType() == Context::Host )
        {
            COMMON_THROWEXCEPTION( "unsupported" )
        }

        if ( source.context()->getType() == Context::Host )
        {
            COMMON_THROWEXCEPTION( "unsupported" )
        }

        ContextData& hostEntry = ( *this )[hostContext];

        if ( ! hostEntry.isValid() )
        {
            hostEntry.reserve( size, 0 );  // reserve it
            hostEntry.copyFrom( source, size );
            hostEntry.setValid( true );
        }

        target.copyFrom( hostEntry, size );
    }
}

/* ---------------------------------------------------------------------------------*/

SyncToken* ContextManager::fetchAsync( ContextData& target, const ContextData& source, size_t size )
{
    try
    {
        SyncToken* token = target.copyFromAsync( source, size );
        return token;
    }
    catch ( common::Exception& ex )
    {
        LAMA_LOG_INFO( logger, target << " async copy from " << source << " not supported" )

        ContextPtr hostContext = Context::getContext( Context::Host );

        if ( target.context()->getType() == Context::Host )
        {
            COMMON_THROWEXCEPTION( "unsupported" )
        }

        if ( source.context()->getType() == Context::Host )
        {
            COMMON_THROWEXCEPTION( "unsupported" )
        }

        ContextData& hostEntry = ( *this )[hostContext];

        if ( ! hostEntry.isValid() )
        {
            hostEntry.reserve( size, 0 );  // reserve it
            hostEntry.copyFrom( source, size );
            hostEntry.setValid( true );
        }

        return target.copyFromAsync( hostEntry, size );
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::copyAllValidEntries( const ContextManager& other, const size_t size )
{
    // each valid data of the other array will be copied into the same context for this array
    size_t nOtherContexts = other.mContextData.size();

    for ( size_t i = 0; i < nOtherContexts; i++ )
    {
        const ContextData& otherData = other.mContextData[i];

        if ( !otherData.isValid() )
        {
            continue; // do not copy any invalid data
        }

        ContextData& data = operator[]( getContextData( otherData.context() ) );

        if ( size > 0 )
        {
            data.reserve( size, 0 );
            fetch( data, otherData, size );
        }

        data.setValid( true );
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::setValidData( ContextPtr context, const ContextManager& other, const size_t size )
{
    ContextData& data = operator[]( context );

    if ( size == 0 )
    {
        return;
    }

    // there must be at least one valid entry
    const ContextData& validData = other.getValidData();
    data.reserve( size, 0 );
    fetch( data, validData, size );
    data.setValid( true );
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::invalidateAll()
{
    size_t noContexts = mContextData.size();
    LAMA_LOG_DEBUG( logger, "invalidate for " << noContexts << " context locations" )

    for ( size_t i = 0; i < noContexts; ++i )
    {
        mContextData[i].setValid( false );
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::swap( ContextManager& other )
{
    // there must be no accesses to the swapped arrays as references would be invalid. */
    COMMON_ASSERT( !locked(), "" )
    COMMON_ASSERT( !other.locked(), "" )
    // due to the pointers swap on vectors is okay
    std::swap( mContextData, other.mContextData );
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::prefetch( ContextPtr context, size_t size )
{
    ContextData& data = operator[]( context );

    if ( data.isValid() || size == 0 )
    {
        return;
    }

    wait();
    data.reserve( size, 0 );
    const ContextData& validEntry = getValidData();
    SyncToken* token = fetchAsync( data, validEntry, size );

    if ( token != NULL )
    {
        // save it, so we can wait for it
        mSyncToken.reset( token );
    }

    // we set data already as valid even if transfer is not finished yet.
    // so any query for valid data must wait for token.
    data.setValid( true );
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::reserve( ContextPtr context, const size_t size, const size_t validSize )
{
    ContextData& data = ( *this )[context];

    // ToDo: must have a write access here COMMON_ASSERT( !data.locked( ContextData::Write ), "no reserve on write locked data." )

    if ( data.isValid() )
    {
        data.reserve( size, validSize );
    }
    else
    {
        data.reserve( size, 0 );  // no valid data
    }
}

} // namespace

