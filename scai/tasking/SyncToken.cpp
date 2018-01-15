/**
 * @file SyncToken.cpp
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
 * @brief Implementation of non-pure methods of abstract class SyncToken.
 * @author Thomas Brandes, Jiri Kraus
 * @date 22.03.2011
 */

// hpp
#include <scai/tasking/SyncToken.hpp>

// internal scai libraries

#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>

#include <memory>
#include <functional>

using std::shared_ptr;

namespace scai
{

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

void SyncToken::pushRoutine( std::function<void()> routine )
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
    if ( currentSyncToken != NULL )
    {
        SCAI_LOG_ERROR( logger, "setCurrent: " << *this << ", but current is: " << *currentSyncToken )
    }

    currentSyncToken = this;
}

void SyncToken::unsetCurrent()
{
    if ( currentSyncToken == NULL )
    {
        SCAI_LOG_WARN( logger, "unset current sync token, not available" )
    }
    else
    {
        SCAI_LOG_INFO( logger, "no more current sync token " << *currentSyncToken )
        currentSyncToken = NULL;
    }
}

SyncToken* SyncToken::getCurrentSyncToken()
{
    return currentSyncToken;
}

// we can rely on the fact that thread-private variable is initialized with NULL

thread_local SyncToken* SyncToken::currentSyncToken = NULL;


/* ----------------------------------------------------------------------- */

} /* end namespace tasking */

} /* end namespace scai */
