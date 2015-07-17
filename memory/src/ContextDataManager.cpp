/**
 * @file ContextDataManager.cpp
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
 * @brief Implementation of methods for class ContextDataManager.
 * @author Thomas Brandes
 * @date 14.03.2011
 */

// hpp

#include <memory/ContextDataManager.hpp>
#include <memory/Context.hpp>

#include <common/Exception.hpp>

namespace memory
{

using tasking::SyncToken;

/* ---------------------------------------------------------------------------------*/

LAMA_LOG_DEF_LOGGER( ContextDataManager::logger, "ContextDataManager" )

/* ---------------------------------------------------------------------------------*/

ContextDataManager::ContextDataManager() :

    mContextData(),
    mSyncToken()

{
    mLock[context::Read] = 0;
    mLock[context::Write] = 0;
    multiContext = false;
    multiThreaded = false;
  
    LAMA_LOG_DEBUG( logger, "ContextDataManager()" )
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::wait()
{
    if ( 0 != mSyncToken.get() )
    {
        LAMA_LOG_DEBUG( logger, "Waiting for SyncToken: " << *mSyncToken )
        mSyncToken.reset(); // waits for transfer and frees resources
    }
}

/* ---------------------------------------------------------------------------------*/

bool ContextDataManager::locked() const
{
    return ( mLock[context::Write] > 0 ) || ( mLock[context::Read] > 0 );
}

/* ---------------------------------------------------------------------------------*/

bool ContextDataManager::locked( AccessKind kind ) const
{
    return ( mLock[kind] > 0 );
}

/* ---------------------------------------------------------------------------------*/

bool ContextDataManager::hasAccessConflict( AccessKind kind ) const
{
    bool conflict = false;

    if ( kind == context::Read )
    {
        conflict = locked( context::Write );
    }
    else if ( kind == context::Write )
    {
        conflict = locked();
    }
    return conflict;
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::lockAccess( AccessKind kind, ContextPtr context )
{
    common::Thread::ScopedLock lock( mAccessMutex );

    common::Thread::Id id = common::Thread::getSelf();

    LAMA_LOG_DEBUG( logger, "lockAccess, kind = " << kind << ", #reads = "
                     << mLock[context::Read] << ", #writes = " << mLock[context::Write] )

    if ( !locked() )
    {
        // first access, so we can set access context

        accessContext  = context;
        accessThread   = id;
        multiContext   = false;
        multiThreaded  = false;

        LAMA_LOG_DEBUG( logger, "first access, set context = " << *context << ", set thread = " << id )

    } 
    else if ( !hasAccessConflict( kind ) )
    {
        // multiple reads

        if ( accessContext.get() != context.get() )
        {
            multiContext = true;
            LAMA_LOG_DEBUG( logger, "multiple Context for read" )
        }

        if ( id != accessThread )
        {
            multiThreaded = true;
            LAMA_LOG_DEBUG( logger, "multiple Thread for read" )
        }
    }
    else if ( multiContext || ( accessContext.get() != context.get() ) )
    {
        // access at different context

        COMMON_THROWEXCEPTION( "Access conflict, kind = " << kind << ", #reads = "
                 << mLock[context::Read] << ", #writes = " << mLock[context::Write] 
                 << ", more than one context" )
    }
    else if ( multiThreaded || id != accessThread )
    {
        // same context but different thread
 
        while ( locked() )
        {
            LAMA_LOG_DEBUG( logger, id << ": wait for free access, blocked by " << accessThread << ", multiple = " << multiThreaded )
            mAccessCondition.wait( lock );
            LAMA_LOG_DEBUG( logger, id << ": have now free access" )
        }

        accessContext  = context;
        accessThread   = id;
        multiContext   = false;
        multiThreaded  = false;
    }
    else
    {
        LAMA_LOG_DEBUG( logger, "same thread, same context, multiThreaded = " << multiThreaded << ", multiContext = " << multiContext )
        // same thread, same context, that is okay for now
    }

    mLock[kind]++;

    LAMA_LOG_DEBUG( logger, "lockAccess done, kind = " << kind << ", #reads = "
                     << mLock[context::Read] << ", #writes = " << mLock[context::Write] )
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::unlockAccess( context::AccessKind kind )
{
    common::Thread::ScopedLock lock( mAccessMutex );

    COMMON_ASSERT_LE( 1, mLock[kind], "release access " << kind << ", never acquired" )

    mLock[kind]--;

    LAMA_LOG_DEBUG( logger, kind << "Access released, #reads = "
                 << mLock[context::Read] << ", #writes = " << mLock[context::Write] )
 
    if ( ( mLock[context::Write] == 0 ) && ( mLock[context::Read] == 0 ) )
    {
        multiContext = false;
        multiThreaded = false;
        LAMA_LOG_DEBUG( logger, common::Thread::getSelf() << ": notify threads waiting for access" )
        mAccessCondition.notifyOne();  // Notify waiting thread
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::releaseAccess( ContextDataIndex index, context::AccessKind kind )
{
    // we should check that this is really the context data for which access was reserved
 
    COMMON_ASSERT_LT( index, mContextData.size(), "Illegal context data index = " << index )

    unlockAccess( kind );
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::purge()
{
    // purge frees all data but keeps the ContextData entries
    wait();
    COMMON_ASSERT( !locked(), "purge on array with access" )

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        ContextData& entry = mContextData[i];
        entry.free();
    }
}

/* ---------------------------------------------------------------------------------*/

ContextDataManager::~ContextDataManager()
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
        LAMA_LOG_INFO( logger, "~ContextDataManager, free " << mContextData[i] )
        mContextData[i].free();
    }

    mContextData.clear();
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextDataManager::findContextData( ContextPtr context ) const
{
    ContextDataIndex i;

    for ( i = 0; i < mContextData.size(); ++i )
    {
        const ContextData& entry = mContextData[i];

        if ( context->canUseMemory( *entry.memory() ) )
        {
            break;
        }
    }

    return i;
}

/* ---------------------------------------------------------------------------------*/

bool ContextDataManager::isValid( ContextPtr context ) const
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

size_t ContextDataManager::capacity( ContextPtr context ) const
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

ContextDataIndex ContextDataManager::getContextData( ContextPtr context )
{
    size_t contextIndex = findContextData( context );

    LAMA_LOG_DEBUG( logger, "contextIndex = " << contextIndex << ", size = " << mContextData.size() )

    // if context is used first time, make new entry for context data

    if ( contextIndex == mContextData.size() )
    {
        if ( contextIndex == 0 )
        {
            mContextData.reserve( MEMORY_MAX_CONTEXTS );
        }

        MemoryPtr memoryPtr = context->getMemory();

        if ( context->getType() == context::Host )
        {
            // For host context we might find more convenient memory

            if ( contextIndex > 0 )
            {
                memoryPtr = mContextData[0].memory()->getContext()->getHostMemory();
            }
        }

        COMMON_ASSERT( memoryPtr, "getMemory failed for context = " << *context )

        mContextData.push_back( ContextData( memoryPtr ) );

        LAMA_LOG_DEBUG( logger, "new context data entry for " << *context << ", index = " << contextIndex )
    }

    return contextIndex;
}

/* ---------------------------------------------------------------------------------*/

ContextData& ContextDataManager::operator[] ( ContextDataIndex index )
{
    COMMON_ASSERT( index < mContextData.size(), "index = " << index << " is illegal index, size = " << mContextData.size() )
    return mContextData[index];
}

/* ---------------------------------------------------------------------------------*/

ContextData& ContextDataManager::operator[] ( ContextPtr context )
{
    return operator[]( getContextData( context ) );
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextDataManager::findValidData() const
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

const ContextData& ContextDataManager::getValidData() const
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

ContextPtr ContextDataManager::getValidContext( const ContextType preferredType )
{
    ContextPtr result;

    for ( size_t index = 0; index < mContextData.size(); ++index )
    {
        const ContextData& entry = mContextData[index];

        if ( entry.isValid() )
        {
            ContextPtr context = entry.memory()->getContext();

            if ( context->getType() == preferredType )
            {
                return context;
            }
            else if ( result )
            {
                // do not overwrite first context found
            }
            else
            {
                result = context;
            }
        }
    }

    return result;  // might be NULL
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextDataManager::acquireAccess( ContextPtr context, AccessKind kind,
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
        LAMA_LOG_DEBUG( logger, "data not valid at " << *data.memory() )
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

    if ( kind == context::Write )
    {
        invalidateAll();        // invalidate all entries
    }

    data.setValid( true );  // for next access the data @ context is valid.

    LAMA_LOG_DEBUG( logger, "acquired access :" << data );
    return index;
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::fetch( ContextData& target, const ContextData& source, size_t size )
{
    try
    {
        target.copyFrom( source, size );
    }
    catch ( common::Exception& ex )
    {
        LAMA_LOG_INFO( logger, target << " copy from " << source << " not supported" )

        // try it via host

        if ( target.memory()->getType() == memtype::HostMemory )
        {
            COMMON_THROWEXCEPTION( "Unsupported: copy to host from: " << *source.memory() )
        }

        if ( source.memory()->getType() == memtype::HostMemory )
        {
            COMMON_THROWEXCEPTION( "Unsupported: copy from host to: " << *target.memory() )
        }

        ContextPtr hostContext = Context::getContext( context::Host );

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

SyncToken* ContextDataManager::fetchAsync( ContextData& target, const ContextData& source, size_t size )
{
    try
    {
        SyncToken* token = target.copyFromAsync( source, size );
        return token;
    }
    catch ( common::Exception& ex )
    {
        LAMA_LOG_INFO( logger, target << " async copy from " << source << " not supported" )

        ContextPtr hostContext = Context::getContext( context::Host );

        if ( target.memory()->getType() == memtype::HostMemory )
        {
            COMMON_THROWEXCEPTION( "unsupported" )
        }

        if ( source.memory()->getType() == memtype::HostMemory )
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

void ContextDataManager::copyAllValidEntries( const ContextDataManager& other, const size_t size )
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

        ContextData& data = operator[]( getContextData( otherData.memory()->getContext() ) );

        if ( size > 0 )
        {
            data.reserve( size, 0 );
            fetch( data, otherData, size );
        }

        data.setValid( true );
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::setValidData( ContextPtr context, const ContextDataManager& other, const size_t size )
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

void ContextDataManager::invalidateAll()
{
    size_t noContexts = mContextData.size();
    LAMA_LOG_DEBUG( logger, "invalidate for " << noContexts << " context locations" )

    for ( size_t i = 0; i < noContexts; ++i )
    {
        mContextData[i].setValid( false );
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::swap( ContextDataManager& other )
{
    // there must be no accesses to the swapped arrays as references would be invalid. */
    COMMON_ASSERT( !locked(), "" )
    COMMON_ASSERT( !other.locked(), "" )
    // due to the pointers swap on vectors is okay
    std::swap( mContextData, other.mContextData );
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::prefetch( ContextPtr context, size_t size )
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

void ContextDataManager::reserve( ContextPtr context, const size_t size, const size_t validSize )
{
    wait();  // valid is checked so outstanding transfers must be finished

    ContextData& data = ( *this )[context];

    // ToDo: must have a write access here COMMON_ASSERT( !data.locked( context::Write ), "no reserve on write locked data." )

    if ( data.isValid() )
    {
        data.reserve( size, validSize );
    }
    else
    {
        data.reserve( size, 0 );  // no valid data
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::resize( const size_t size, const size_t validSize )
{
    COMMON_ASSERT( !locked(), "Array is locked, no resize possible" )

    wait();  // valid is checked so outstanding transfers must be finished

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        ContextData& data = mContextData[i];

        if ( data.isValid() )
        {
            data.reserve( size, validSize );
        }
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::writeAt( std::ostream& stream ) const
{
    if ( mContextData.size() == 0 )
    {
        stream << "no data";
        return;
    }

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        if ( i > 0 )
        {
            stream << ", ";
        }

        stream << mContextData[i];
    }
}

/* ---------------------------------------------------------------------------------*/


} // namespace

