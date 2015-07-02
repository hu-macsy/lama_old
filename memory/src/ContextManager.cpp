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
 * @since 1.0.0
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

void ContextManager::wait() 
{
    if( 0 != mSyncToken.get() )
    {
        LAMA_LOG_DEBUG( logger, "Waiting for SyncToken: " << *mSyncToken )

        mSyncToken.reset(); // waits for transfer and frees resources
    }
}

/* ---------------------------------------------------------------------------------*/

bool ContextManager::locked() const
{
    bool locked = false;

    // check for existing read / write lock

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        if ( mContextData[i]->locked() )
        {
            return true;
        }
    }

    return false;
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::purge()
{
    // purge frees all data but keeps the ContextData entries

    wait();

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        ContextData& entry = (*this)[i];

        COMMON_ASSERT( !entry.locked(), "purge, but locked " << entry );
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

    if( isLocked )
    {
        LAMA_LOG_WARN( logger, "Destructor on read/write locked array: " )
    }

    for ( size_t i = 0; i < mContextData.size(); i++ )
    {
        LAMA_LOG_INFO( logger, "~ContextManager, delete " << *mContextData[i] )
        delete mContextData[i];
    }

    mContextData.clear();
}

/* ---------------------------------------------------------------------------------*/

ContextDataRef ContextManager::getContextData( ContextPtr context, AccessKind kind )
{
    bool found = false;

    size_t contextIndex;

    LAMA_LOG_DEBUG( logger, "check access for " << *context )

    for( size_t i = 0; i < mContextData.size(); ++i )
    {
        ContextData& entry = *mContextData[i];

        if ( context->canUseData( *entry.context ) )
        {
            contextIndex = i;
            found = true;
        }

        if( !entry.valid )
        {
            // Assertions might be removed later
            // Also be careful if getContextData

            COMMON_ASSERT( !entry.locked(), "read or write access on non valid location" )
        }

        if ( entry.locked( ContextData::Read ) )
        {
            COMMON_ASSERT( !entry.locked( ContextData::Write ), "write and read access" )

            if( kind == ContextData::Write )
            {
                COMMON_THROWEXCEPTION( "try to get write access on read locked array " )
            }
        }
        else if( entry.locked( ContextData::Write ) )
        {
            COMMON_THROWEXCEPTION( "no further access on write locked array " )
        }
    }

    // if context is used first time, make new entry for context data

    if ( !found )
    {
        contextIndex = mContextData.size();

        if ( contextIndex == 0 )
        {
            mContextData.reserve( LAMA_MAX_CONTEXTS );
        }
            
        mContextData.push_back( new ContextData( context ) );

        LAMA_LOG_DEBUG( logger, "new context data entry for " << *context << ", index = " << contextIndex )
    }

    // increase the lock counter already here, so it is clear that there is a reference

    // better not  mContextData[contextIndex].addAccess( kind );

    return contextIndex;
}

/* ---------------------------------------------------------------------------------*/

ContextData& ContextManager::operator[] (ContextDataRef ref )
{
    COMMON_ASSERT( ref < mContextData.size(), "ref = " << ref << " is illegal index, size = " << mContextData.size() )

    return *mContextData[ref];
}

/* ---------------------------------------------------------------------------------*/

ContextDataRef ContextManager::getValidData()
{
    bool found = false;

    for ( size_t i = 0; i < mContextData.size(); ++i )
    {
        ContextData& entry = *mContextData[i];

        if ( !entry.valid )
        {
            continue;
        }

        //  so we have found a context with valid data

        return i;
    }

    COMMON_THROWEXCEPTION( "no valid data found" )
}

/* ---------------------------------------------------------------------------------*/

ContextDataRef ContextManager::acquireAccess( ContextPtr context, AccessKind kind, 
                                              size_t allocSize, size_t validSize )
{
    COMMON_ASSERT( context, "NULL pointer for context" )

    LAMA_LOG_DEBUG( logger, "acquire access on " << *context << ", kind = " << kind 
                            << ", allocSize = " << allocSize << ", validSize = " << validSize )

    ContextDataRef ref = getContextData( context, kind );

    ContextData& data = (*this)[ref];

    wait();

    // fetch only if size > 0, there might be no valid location for mSize == 0

    if( !data.valid && allocSize > 0 )
    {
        LAMA_LOG_DEBUG( logger, "data not valid at " << *data.context )

        // make sure that we have enough memory on the target context
        // old data is invalid so it must not be saved.

        data.reserve( allocSize, 0 );  // do not save any old values

        if ( validSize )
        {
            ContextDataRef validIndex = getValidData();

            ContextData& validEntry = (*this)[validIndex];

            LAMA_LOG_INFO( logger, "valid data here: " << validEntry )

            try
            {
                data.copyFrom( validEntry, validSize );
            }
            catch ( common::Exception& ex )
            {
                LAMA_LOG_INFO( logger, data << " copy from " << validEntry << " not supported" )

                // try it via host
           
                ContextPtr hostContext = Context::getContext( Context::Host );
                
                if ( *data.context == *hostContext )
                {
                     COMMON_THROWEXCEPTION( "unsupported" )
                }
                if ( *validEntry.context == *hostContext )
                {
                     COMMON_THROWEXCEPTION( "unsupported" )
                }

                ContextDataRef hostIndex = getContextData( hostContext, ContextData::Read );
                ContextData& hostEntry = (*this)[hostIndex];

                if ( ! hostEntry.valid )
                {
                    hostEntry.reserve( allocSize, 0 );  // reserve it
                    hostEntry.copyFrom( validEntry, validSize );
                    hostEntry.setValid( true );
                }

                data.copyFrom( hostEntry, validSize );
            }
        }
    }

    invalidateAll(); // invalidate all entries

    data.setValid( true );  // for next access the data @ context is valid.
    data.addLock( kind );

    LAMA_LOG_DEBUG( logger, "acquired access :" << data );

    return ref;
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::copyAllValidEntries( const ContextManager& other, const size_t size )
{
    // each valid data of the other array will be copied into the same context for this array

    size_t nOtherContexts = other.mContextData.size();

    for ( size_t i = 0; i < nOtherContexts; i++ )
    {
        const ContextData& otherEntry = *other.mContextData[i];

        if( !otherEntry.valid )
        {
            continue; // do not copy any invalid data
        }

        ContextData& myEntry = operator[]( getContextData( otherEntry.context, ContextData::Write ) );

        bool keepFlag = false;   // do not keep values of the current data

        // ToDo: myEntry.reserve( size, false );

        // and then copy the data within the same context

        if( size > 0 )
        {
            // ToDo: myEntry.copyFrom( otherEntry, size );   // copy on same device
        }

        myEntry.valid = true;

        // ToDo: myEntry.release( ContextData::Write );
    }
}

/* ---------------------------------------------------------------------------------*/

void ContextManager::invalidateAll()
{
    size_t noContexts = mContextData.size();

    LAMA_LOG_DEBUG( logger, "invalidate for " << noContexts << " context locations" )

    for( size_t i = 0; i < noContexts; ++i )
    {
        mContextData[i]->valid = false;
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

} // namespace

