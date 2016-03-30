/**
 * @file MICSyncToken.cpp
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
 * @brief MICSyncToken.cpp
 * @author Thomas Brandes
 * @date 03.07.2013
 * @since 1.1.0
 */

// hpp
#include <scai/tasking/mic/MICSyncToken.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace tasking
{

MICSyncToken::MICSyncToken() : mSignal( -1 )
{
    SCAI_LOG_DEBUG( logger, "MICSyncToken generated" )
}

MICSyncToken::~MICSyncToken()
{
    wait();
}

void MICSyncToken::wait()
{
    if( isSynchronized() )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger, "wait for offload computation by signal" )

    if( mSignal >= 0 )
    {
        // SCAI_REGION( "MIC.offloadSynchronize" )

        // // finally called functions might also need the context, e.g. unbindTexture

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

    if( isSynchronized() )
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
