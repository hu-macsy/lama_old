/**
 * @file MICSyncToken.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief MICSyncToken.cpp
 * @author Thomas Brandes
 * @date 03.07.2013
 */

// hpp
#include <scai/tasking/mic/MICSyncToken.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace tasking
{

MICSyncToken::MICSyncToken( int device ) : mSignal( -1 ), mDevice( device )
{
    SCAI_LOG_DEBUG( logger, "MICSyncToken generated" )
}

MICSyncToken::~MICSyncToken()
{
    wait();
}

void MICSyncToken::wait()
{
    if ( isSynchronized() )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger, "wait for offload computation by signal" )

    if ( mSignal >= 0 )
    {
        // SCAI_REGION( "MIC.offloadSynchronize" )
#pragma offload target( mic:mDevice ) wait(mSignal)
        {
        }
    }

    setSynchronized();
}

void MICSyncToken::setSignal( int signal )
{
    SCAI_ASSERT_ERROR( !isSynchronized(), "cannot set signal as SyncToken is already synchronized" )
    SCAI_ASSERT_ERROR( mSignal < 0, "signal already set, cannot handle multiple signals" )
    // this signal will be used to wait for synchronization
    mSignal = signal;
}

bool MICSyncToken::probe() const
{
    // No idea what to do so wait

    // wait();
    if ( isSynchronized() )
    {
        return true;
    }

    return false;
}

MICSyncToken* MICSyncToken::getCurrentSyncToken()
{
    SyncToken* syncToken = SyncToken::getCurrentSyncToken();

    if ( syncToken == NULL )
    {
        return NULL;
    }

    // make a dynamic CAST
    MICSyncToken* micSyncToken = dynamic_cast<MICSyncToken*>( syncToken );

    // If the current sync token is not a CUDA stream token it is very likely an error

    if ( micSyncToken == NULL )
    {
        SCAI_LOG_ERROR( logger, "Current sync token = " << *syncToken << " not MICSyncToken as expected" )
    }

    // But might not be too serious so return NULL that forces synchronous execution
    return micSyncToken;
}

} /* end namespace tasking */

} /* end namespace scai */
