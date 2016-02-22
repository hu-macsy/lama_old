/**
 * @file SyncToken.cpp
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
 * @brief Implementation of non-pure methods of abstract class SyncToken.
 * @author Thomas Brandes, Jiri Kraus
 * @date 22.03.2011
 */

// hpp
#include <scai/tasking/SyncToken.hpp>

// internal scai libraries

#include <scai/common/Thread.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>

namespace scai
{

using common::shared_ptr;

namespace tasking
{

SCAI_LOG_DEF_LOGGER( SyncToken::logger, "SyncToken" )

/* ------------------------------------------------------------------------ */

SyncToken::CGuard::CGuard()
{
}

SyncToken::CGuard::~CGuard()
{
    // this destructor is called at the end of the program
    // Give an error message if not all SyncTokens have been deleted

    if ( countSyncToken )
    {
        SCAI_LOG_ERROR( logger, "Remaining SyncToken (undeleted) = " << countSyncToken );
    }
}

/* ------------------------------------------------------------------------ */

int SyncToken::countSyncToken = 0;

SyncToken::CGuard SyncToken::cguard;

/* ------------------------------------------------------------------------ */

SyncToken::SyncToken()
    : mSynchronized( false )
{
    SCAI_LOG_DEBUG( logger, "SyncToken constructed" )

    countSyncToken++;
}

SyncToken::~SyncToken()
{
    SCAI_LOG_DEBUG( logger, "~SyncToken" )

    if ( !mSynchronized )
    {
        SCAI_LOG_WARN( logger, "no synchronization called on SyncToken" )
    }

    countSyncToken--;
}

/* ------------------------------------------------------------------------ */

void SyncToken::writeAt( std::ostream& stream ) const
{
    stream << "SyncToken( synchronized = " << mSynchronized << " )";
}

/* ------------------------------------------------------------------------ */

void SyncToken::pushToken( shared_ptr<SyncTokenMember> member )
{
    SCAI_ASSERT( member.get(), "NULL token cannot be pushed for synchronization." )

    if ( mSynchronized )
    {
        SCAI_LOG_DEBUG( logger, *this << ": push token not done, already synchronized" )
    }
    else
    {
        SCAI_LOG_INFO( logger, *this << ": push token, will be freed at synchronization" )

        // take ownership of the token so it is not deleted before synchronization

        mTokens.push_back( member );
    }
}

/* ----------------------------------------------------------------------- */

void SyncToken::pushRoutine( common::function<void()> routine )
{
    mSynchronizedFunctions.push_back( routine );
}

/* ----------------------------------------------------------------------- */

bool SyncToken::isSynchronized() const
{
    return mSynchronized;
}

/* ----------------------------------------------------------------------- */

void SyncToken::setSynchronized()
{
    if ( mSynchronized )
    {
        COMMON_THROWEXCEPTION( *this << " is already synchronized" )
    }

    SCAI_LOG_INFO( logger, "setSynchronized, free " << mTokens.size() << " SyncTokenMember "
                   << " and call " << mSynchronizedFunctions.size() << " clean functions" )

    mSynchronized = true;

    // after synchronization we can give up ownership of tokens

    mTokens.clear();

    for ( size_t i = 0; i < mSynchronizedFunctions.size(); ++i )
    {
        mSynchronizedFunctions[i]();
    }
}

/* ----------------------------------------------------------------------- */

void SyncToken::setCurrent()
{
    SyncToken* current = currentSyncToken.get();

    if ( current != NULL )
    {
        SCAI_LOG_ERROR( logger, "setCurrent: " << *this << ", but current is: " << *current )
    }

    currentSyncToken.set( this );
}

void SyncToken::unsetCurrent()
{
    SyncToken* current = currentSyncToken.get();

    if ( current == NULL )
    {
        SCAI_LOG_WARN( logger, "unset current sync token, not available" )
    } 
    else 
    {
        SCAI_LOG_INFO( logger, "no more current sync token " << *current )
        currentSyncToken.set( NULL );
    }
}

SyncToken* SyncToken::getCurrentSyncToken()
{
    return currentSyncToken.get();
}

// we can rely on the fact that thread-private variable is initialized with NULL 

common::ThreadPrivatePtr<SyncToken> SyncToken::currentSyncToken;

/* ----------------------------------------------------------------------- */

} /* end namespace tasking */

} /* end namespace scai */
