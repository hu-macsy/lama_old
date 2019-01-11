/**
 * @file ContextDataManager.cpp
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
 * @brief Implementation of methods for class ContextDataManager.
 * @author Thomas Brandes
 * @date 14.03.2011
 */

// hpp
#include <scai/hmemo/ContextDataManager.hpp>

// local library
#include <scai/hmemo/Context.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/terminate.hpp>

namespace scai
{

using common::AccessKind;

namespace hmemo
{

using tasking::SyncToken;

/* ---------------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( ContextDataManager::logger, "ContextDataManager" )

/* ---------------------------------------------------------------------------------*/

ContextDataManager::ContextDataManager() :

    mContextData(),
    mSyncToken()

{
    mLock[static_cast<int>( AccessKind::Read  )] = 0;
    mLock[static_cast<int>( AccessKind::Write )] = 0;
    multiContext = false;
    multiThreaded = false;
    SCAI_LOG_DEBUG( logger, "ContextDataManager()" )
}

/* ---------------------------------------------------------------------------------*/

ContextDataManager::ContextDataManager( ContextDataManager&& other ) noexcept :

    mContextData( std::move( other.mContextData ) )

{
    SCAI_LOG_INFO( logger, "move constructor: " 
                          << ", this context data has now " << mContextData.size() << " entries"
                          << ", other context data has now " << other.mContextData.size() << " entries" )

    const int pWrite = static_cast<int>( AccessKind::Write );
    const int pRead = static_cast<int>( AccessKind::Read );

    // check for current accesses, if yes terminate as there is no way to recover

    if ( ( other.mLock[pWrite] != 0 ) )
    {
        SCAI_TERMINATE( "SERIOUS: move operation on write locked array, causes termination" )
    }

    if ( ( other.mLock[pRead] != 0 ) )
    {
        SCAI_TERMINATE( "SERIOUS: move operation on read locked array, causes termination" )
    }

    mLock[pRead] = 0;
    mLock[pWrite] = 0;
    multiContext = false;
    multiThreaded = false;
}

/* ---------------------------------------------------------------------------------*/

ContextDataManager& ContextDataManager::operator=( ContextDataManager&& other ) 
{
    SCAI_LOG_INFO( logger, "move asignment called for context data manager: " 
                          << ", this context data has " << mContextData.size() << " entries"
                          << ", other context data has " << other.mContextData.size() << " entries" )

    const int pWrite = static_cast<int>( AccessKind::Write );
    const int pRead = static_cast<int>( AccessKind::Read );

    // check for current accesses, if yes terminate as there is no way to recover

    if ( ( other.mLock[pWrite] != 0 ) )
    {
        SCAI_TERMINATE( "SERIOUS: move operation on write locked array, causes termination" )
    }

    if ( ( other.mLock[pRead] != 0 ) )
    {
        SCAI_TERMINATE( "SERIOUS: move operation on read locked array, causes termination" )
    }

    mLock[pRead] = 0;
    mLock[pWrite] = 0;
    multiContext = false;
    multiThreaded = false;

    // Note: existing entries in mContextData will be freed

    mContextData = std::move( other.mContextData );

    SCAI_LOG_INFO( logger, "move asignment done for context data manager: " 
                          << ", this context data has " << mContextData.size() << " entries"
                          << ", other context data has " << other.mContextData.size() << " entries" )

    return *this;
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::wait()
{
    if ( 0 != mSyncToken.get() )
    {
        SCAI_LOG_DEBUG( logger, "Waiting for SyncToken: " << *mSyncToken )
        mSyncToken.reset(); // waits for transfer and frees resources
    }
}

/* ---------------------------------------------------------------------------------*/

int ContextDataManager::locked() const
{
    return mLock[static_cast<int>( AccessKind::Write )] + mLock[static_cast<int>( AccessKind::Read )];
}

/* ---------------------------------------------------------------------------------*/

int ContextDataManager::locked( AccessKind kind ) const
{
    return mLock[static_cast<int>( kind )];
}

/* ---------------------------------------------------------------------------------*/

bool ContextDataManager::hasAccessConflict( AccessKind kind ) const
{
    bool conflict = false;

    if ( kind == AccessKind::Read )
    {
        conflict = locked( AccessKind::Write );
    }
    else if ( kind == AccessKind::Write )
    {
        conflict = locked();
    }

    return conflict;
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::lockAccess( AccessKind kind, ContextPtr context )
{
    const int pKind  = static_cast<int>( kind );
    const int pRead  = static_cast<int>( AccessKind::Read );
    const int pWrite = static_cast<int>( AccessKind::Write );

    std::unique_lock<std::mutex> lock( mAccessMutex );
    std::thread::id id = std::this_thread::get_id();
    SCAI_LOG_DEBUG( logger, "lockAccess, kind = " << kind << ", #reads = "
                    << mLock[pRead] << ", #writes = " << mLock[pWrite] )

    if ( !locked() )
    {
        // first access, so we can set access context
        accessContext  = context;
        accessThread   = id;
        multiContext   = false;
        multiThreaded  = false;
        SCAI_LOG_DEBUG( logger, "first access, set context = " << *context << ", set thread = " << id )
    }
    else if ( locked( AccessKind::Write ) > 0 )
    {
        COMMON_THROWEXCEPTION( "Access conflict, no further access after a write access, data might have been reallocated" )
    }
    else if ( !hasAccessConflict( kind ) )
    {
        // multiple reads
        if ( accessContext.get() != context.get() )
        {
            multiContext = true;
            SCAI_LOG_DEBUG( logger, "multiple Context for read" )
        }

        if ( id != accessThread )
        {
            multiThreaded = true;
            SCAI_LOG_DEBUG( logger, "multiple Thread for read" )
        }
    }
    else if ( multiContext || ( accessContext.get() != context.get() ) )
    {
        // access at different context
        COMMON_THROWEXCEPTION( "Access conflict, kind = " << kind << ", #reads = "
                               << mLock[pRead] << ", #writes = " << mLock[pWrite]
                               << ", more than one context" )
    }
    else if ( multiThreaded || id != accessThread )
    {
        // same context but different thread
        while ( locked() )
        {
            SCAI_LOG_DEBUG( logger, id << ": wait for free access, blocked by " << accessThread << ", multiple = " << multiThreaded )
            mAccessCondition.wait( lock );
            SCAI_LOG_DEBUG( logger, id << ": have now free access" )
        }

        accessContext  = context;
        accessThread   = id;
        multiContext   = false;
        multiThreaded  = false;
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "same thread, same context, multiThreaded = " << multiThreaded << ", multiContext = " << multiContext )
        // same thread, same context, that is okay for now
    }

    mLock[pKind]++;
    SCAI_LOG_DEBUG( logger, "lockAccess done, kind = " << kind << ", #reads = "
                    << mLock[pRead] << ", #writes = " << mLock[pWrite] )
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::unlockAccess( AccessKind kind )
{
    const int pKind  = static_cast<int>( kind );
    const int pRead  = static_cast<int>( AccessKind::Read );
    const int pWrite = static_cast<int>( AccessKind::Write );

    std::unique_lock<std::mutex> lock( mAccessMutex );
    SCAI_ASSERT_LE( 1, mLock[pKind], "release access " << kind << ", never acquired" )
    mLock[pKind]--;
    SCAI_LOG_DEBUG( logger, kind << "Access released, #reads = "
                    << mLock[pRead] << ", #writes = " << mLock[pWrite] )

    if ( ( mLock[pWrite] == 0 ) && ( mLock[pRead] == 0 ) )
    {
        multiContext = false;
        multiThreaded = false;
        SCAI_LOG_DEBUG( logger, std::this_thread::get_id() << ": notify threads waiting for access" )
        mAccessCondition.notify_one();  // Notify waiting thread
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::releaseAccess( ContextDataIndex index, AccessKind kind )
{
    // we should check that this is really the context data for which access was reserved
    SCAI_ASSERT_LT( index, mContextData.size(), "Illegal context data index = " << index )
    unlockAccess( kind );
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::purge()
{
    // purge frees all data but keeps the ContextData entries
    wait();
    SCAI_ASSERT( !locked(), "purge on array with access" )

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

    if ( SCAI_LOG_WARN_ON( logger ) )
    {
        isLocked = locked();
    }

    // give a warning, never throw an exception in destructor ( might crash )

    if ( isLocked )
    {
        SCAI_LOG_WARN( logger, "Destructor on read/write locked array: " )
    }

    for ( size_t i = 0; i < mContextData.size(); i++ )
    {
        SCAI_LOG_INFO( logger, "~ContextDataManager, free " << mContextData[i] )
        mContextData[i].free();
    }

    mContextData.clear();
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextDataManager::findMemoryData( MemoryPtr memory ) const
{
    const ContextDataIndex nindex = mContextData.size();
    ContextDataIndex i;

    // search for entry

    for ( i = 0; i < nindex; ++i )
    {
        const ContextData& entry = mContextData[i];

        // comparison for memory is just pointer equality

        if ( entry.getMemoryPtr().get() == memory.get() )
        {
            break;
        }
    }

    return i;
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextDataManager::findContextData( ContextPtr context ) const
{
    const ContextDataIndex nindex = mContextData.size();
    ContextDataIndex i;

    // search for entry

    for ( i = 0; i < nindex; ++i )
    {
        const ContextData& entry = mContextData[i];

        // can this entry be used

        if ( context->canUseMemory( entry.getMemory() ) )
        {
            // First touch policy, so no further search, not even for valid entries
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
    SCAI_LOG_DEBUG( logger, "contextIndex = " << contextIndex << ", size = " << mContextData.size() )

    // if context is used first time, make new entry for context data

    if ( contextIndex == mContextData.size() )
    {
        if ( contextIndex == 0 )
        {
            mContextData.reserve( MEMORY_MAX_CONTEXTS );
        }

        MemoryPtr memoryPtr = context->getMemoryPtr();

        if ( context->getType() == common::ContextType::Host )
        {
            // For host context we might find more convenient memory by first touch
            if ( contextIndex > 0 )
            {
                memoryPtr = mContextData[0].getMemory().getContext().getHostMemoryPtr();
            }
        }

        SCAI_ASSERT( memoryPtr, "getMemory failed for context = " << *context )
        mContextData.push_back( ContextData( memoryPtr ) );
        SCAI_LOG_DEBUG( logger, "new context data entry for " << *context << ", index = " << contextIndex )
    }

    return contextIndex;
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextDataManager::getMemoryData( MemoryPtr memoryPtr )
{
    size_t contextIndex = findMemoryData( memoryPtr );
    SCAI_LOG_DEBUG( logger, "contextIndex = " << contextIndex << ", size = " << mContextData.size() )

    // if this memory is used first time, make new entry for context data

    if ( contextIndex == mContextData.size() )
    {
        if ( contextIndex == 0 )
        {
            mContextData.reserve( MEMORY_MAX_CONTEXTS );
        }

        mContextData.push_back( ContextData( memoryPtr ) );
        SCAI_LOG_DEBUG( logger, "new data entry for " << *memoryPtr << ", index = " << contextIndex )
    }

    return contextIndex;
}

/* ---------------------------------------------------------------------------------*/

ContextData& ContextDataManager::operator[] ( ContextDataIndex index )
{
    SCAI_ASSERT( index < mContextData.size(), "index = " << index << " is illegal index, size = " << mContextData.size() )
    return mContextData[index];
}

/* ---------------------------------------------------------------------------------*/

ContextData& ContextDataManager::operator[] ( ContextPtr context )
{
    return operator[]( getContextData( context ) );
}

/* ---------------------------------------------------------------------------------*/

ContextData& ContextDataManager::operator[] ( MemoryPtr memory )
{
    return operator[]( getMemoryData( memory ) );
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextDataManager::findValidData() const
{
    ContextDataIndex index;

    for ( index = 0; index < mContextData.size(); ++index )
    {
        const ContextData& entry = mContextData[index];
        SCAI_LOG_INFO( logger, "check valid: " << entry )

        if ( entry.isValid() )
        {
            SCAI_LOG_INFO( logger, "found valid entry at index = " << index << ": " << entry )
            break;
        }
    }

    if ( index == mContextData.size() )
    {
        SCAI_LOG_INFO( logger, "no valid data found, index = " << index )
    }

    return index;
}

/* ---------------------------------------------------------------------------------*/

ContextPtr ContextDataManager::getValidContext( const ContextPtr prefContext ) const
{
    ContextPtr result;

    for ( size_t index = 0; index < mContextData.size(); ++index )
    {
        const ContextData& entry = mContextData[index];

        if ( entry.isValid() )
        {
            ContextPtr context = entry.getMemory().getContextPtr();

            if ( context.get() == prefContext.get() )
            {
                return prefContext;
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

    if ( !result.get() )
    {
        // might happen for uninitialized arrays
        SCAI_LOG_INFO( logger, "no valid context found for HArray, take first touch context" )
        return getFirstTouchContextPtr();
    }

    return result;  // must not be NULL
}

/* ---------------------------------------------------------------------------------*/

ContextPtr ContextDataManager::getFirstTouchContextPtr() const
{
    if ( mContextData.size() == 0 )
    {
        return Context::getHostPtr();
    }

    const ContextData& entry = mContextData[0];

    return entry.getMemory().getContextPtr();
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextDataManager::acquireAccess( ContextPtr context, AccessKind kind,
        size_t allocSize, size_t validSize )
{
    SCAI_ASSERT( context, "NULL pointer for context" )
    lockAccess( kind, context );
    SCAI_LOG_DEBUG( logger, "acquire access on " << *context << ", kind = " << kind
                    << ", allocSize = " << allocSize << ", validSize = " << validSize )
    ContextDataIndex index = getContextData( context );
    ContextData& data = ( *this )[index];
    wait();

    SCAI_LOG_DEBUG( logger, "index = " << index << " for data, valid = " << data.isValid() )

    // fetch only if size > 0, there might be no valid location for mSize == 0

    if ( !data.isValid() && allocSize > 0 )
    {
        SCAI_LOG_DEBUG( logger, "data not valid at " << data.getMemory() )
        // make sure that we have enough memory on the target context
        // old data is invalid so it must not be saved.
        size_t validSizeHere = 0;   // no values to save here in this context
        bool inUse = false;         // no other access here, so realloc is okay
        data.reserve( allocSize, validSizeHere, inUse );  // do not save any old values, no other use

        if ( validSize )
        {
            ContextDataIndex validIndex = findValidData();

            if ( validIndex < mContextData.size() )
            {
                const ContextData& validEntry = mContextData[ validIndex ];
                SCAI_LOG_INFO( logger, "valid data here: " << validEntry )
                fetch( data, validEntry, validSize );
            }
            else if ( kind == AccessKind::Read )
            {
                SCAI_LOG_WARN( logger, "acquired read access for uninitialized array" )
            }
        }
    }

    if ( kind == AccessKind::Write )
    {
        invalidateAll();        // invalidate all entries
    }

    data.setValid( true );  // for next access the data @ context is valid.
    SCAI_LOG_DEBUG( logger, "acquired access :" << data );
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
        SCAI_LOG_INFO( logger, target << " copy from " << source << " not supported" )

        // try it via host

        if ( target.getMemory().getType() == MemoryType::HostMemory )
        {
            COMMON_THROWEXCEPTION( "Unsupported: copy to host from: " << source.getMemory() )
        }

        if ( source.getMemory().getType() == MemoryType::HostMemory )
        {
            COMMON_THROWEXCEPTION( "Unsupported: copy from host to: " << target.getMemory() )
        }

        ContextPtr hostContextPtr = Context::getHostPtr();
        ContextData& hostEntry = ( *this )[hostContextPtr];

        if ( ! hostEntry.isValid() )
        {
            hostEntry.reserve( size, 0, false );  // reserve it
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
        SCAI_LOG_ERROR( logger, target << " async copy from " << source << " to " << target
                        << " not supported, or has thrown exception" )

        ContextPtr hostContextPtr = Context::getHostPtr();

        if (     target.getMemory().getType() == MemoryType::HostMemory
                 ||  source.getMemory().getType() == MemoryType::HostMemory )
        {
            COMMON_THROWEXCEPTION( "copyAsync from " << source << " to " << target
                                   << " must be supported with HostMemory, exception = " << ex.what() )
        }

        // as neither source nor target are HostMemory, try it via host

        ContextData& hostEntry = ( *this )[hostContextPtr];

        if ( ! hostEntry.isValid() )
        {
            hostEntry.reserve( size, 0, false );  // reserve it
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

        ContextData& data = operator[]( getContextData( otherData.getMemory().getContextPtr() ) );

        if ( size > 0 )
        {
            data.reserve( size, 0, false );
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

    bool inUse = false;   // no other access
    data.reserve( size, 0, inUse );   // reserve even if no valid data is available
    ContextDataIndex validIndex = other.findValidData();

    if ( validIndex < other.mContextData.size() )
    {
        // there must be at least one valid entry
        const ContextData& validData = other.mContextData[ validIndex ];
        fetch( data, validData, size );
        data.setValid( true );
    }
    else
    {
        SCAI_LOG_WARN( logger, "cannot set valid data as no valid data is available (uninitialized)" )
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::invalidateAll()
{
    size_t noContexts = mContextData.size();
    SCAI_LOG_DEBUG( logger, "invalidate for " << noContexts << " context locations" )

    for ( size_t i = 0; i < noContexts; ++i )
    {
        mContextData[i].setValid( false );
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::swap( ContextDataManager& other )
{
    // there must be no accesses to the swapped arrays as references would be invalid. */
    SCAI_ASSERT( !locked(), "" )
    SCAI_ASSERT( !other.locked(), "" )
    // due to the pointers swap on vectors is okay
    std::swap( mContextData, other.mContextData );
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::prefetch( ContextPtr context, size_t size )
{
    ContextData& data = operator[]( context );
    SCAI_LOG_DEBUG( logger, "prefetch to " << *context << ", valid = " << data.isValid() << ", size = " << size )

    if ( data.isValid() || size == 0 )
    {
        return;
    }

    wait();   // wait on previous transfers
    data.reserve( size, 0, false );
    ContextDataIndex validIndex = findValidData();

    if ( validIndex < mContextData.size() )
    {
        const ContextData& validEntry = mContextData[ validIndex ];
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
    else
    {
        // no valid data is not serious
        // might be worth a warning but happens very often in LAMA
        SCAI_LOG_INFO( logger, "prefetch on array with no valid data" )
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::reserve( ContextPtr context, const size_t size, const size_t validSize )
{
    wait();  // valid is checked so outstanding transfers must be finished
    ContextData& data = ( *this )[context];
    bool inUse = false;

    // ToDo: must have a write access here SCAI_ASSERT( !data.locked( AccessKind::Write ), "no reserve on write locked data." )

    if ( data.isValid() )
    {
        data.reserve( size, validSize, inUse );
    }
    else
    {
        data.reserve( size, 0, inUse );  // no valid data
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextDataManager::resize( const size_t size, const size_t validSize )
{
    bool inUse = locked();   // do not 'really' resize for if there are any write/read accesses
    wait();  // valid is checked so outstanding transfers must be finished

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        ContextData& data = mContextData[i];

        if ( data.isValid() )
        {
            data.reserve( size, validSize, inUse );
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

void ContextDataManager::init( const void* data, const size_t size )
{
    // choose as context for valid entries the first touch context
    ContextPtr ctx = getFirstTouchContextPtr();
    // acquire write only access on this context, validSize = 0;
    ContextDataIndex index = acquireAccess( ctx, AccessKind::Write, size, 0 );
    // make a fictive ContextData entry for host data
    ContextData hostEntry( Context::getHostPtr()->getMemoryPtr() );
    hostEntry.setRef( const_cast<void*>( data ), size );
    mContextData[index].reserve( size, 0, false );
    // copy the host data to the destination
    mContextData[index].copyFrom( hostEntry, size );
    releaseAccess( index, AccessKind::Write );
}

void ContextDataManager::getData( void* data, const size_t offset, const size_t size )
{
    // get a valid context, preferred host

    ContextPtr ctx     = Context::getHostPtr();
    const Memory& hostMem = *ctx->getMemoryPtr();

    ctx = getValidContext( ctx );

    ContextDataIndex index = acquireAccess( ctx, AccessKind::Read, offset + size, offset + size );

    const Memory& mem = mContextData[index].getMemory();

    const char* ptr = static_cast<const char*>( mContextData[index].get() );

    SCAI_LOG_INFO( logger, "getData( offset = " << offset << ", size = " << size << " ), index = " << index << ", mem = " << mem << " copies to " << hostMem )

    mem.memcpyTo( hostMem, data, ptr + offset, size );

    releaseAccess( index, common::AccessKind::Read );
}

void ContextDataManager::setData( const void* data, const size_t offset, const size_t dataSize, const size_t allocSize )
{
    // get a valid context, preferred host

    ContextPtr ctx     = Context::getHostPtr();
    const Memory& hostMem = *ctx->getMemoryPtr();

    ctx = getValidContext( ctx );

    ContextDataIndex index = acquireAccess( ctx, AccessKind::Write, allocSize, allocSize );

    SCAI_LOG_INFO( logger, "setData( offset = " << offset << ", size = " << dataSize << " ), index = " << index << ", data: " << mContextData[index] )

    const Memory& mem = mContextData[index].getMemory();

    char* ptr = static_cast<char*>( mContextData[index].get() );

    SCAI_LOG_DEBUG( logger, "setData( offset = " << offset << ", size = " << dataSize << " ), index = " << index << ", mem = " << mem << " copies from " << hostMem )

    mem.memcpyFrom( ptr + offset, hostMem, data, dataSize );

    releaseAccess( index, common::AccessKind::Write );
}

/* ---------------------------------------------------------------------------------*/

} /* end namespace hmemo */

} /* end namespace scai */
