/**
 * @file MICSyncToken.cpp
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
 * @brief MICSyncToken.cpp
 * @author Thomas Brandes
 * @date 03.07.2013
 * @since 1.1.0
 */

// hpp
#include <lama/mic/MICSyncToken.hpp>

// others
#include <lama/tracing.hpp>
#include <lama/exception/LAMAAssert.hpp>

namespace lama
{

MICSyncToken::MICSyncToken( MICContextPtr micContext )
    : mMICContext( micContext ),
      mSignal( -1 )
{
    LAMA_LOG_DEBUG( logger, "MICSyncToken for " << *micContext << " generated" )
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

    LAMA_LOG_DEBUG( logger, "wait for offload computation by signal" )

    if ( mSignal >= 0 )
    {
        LAMA_REGION( "MIC.offloadSynchronize" )

        #pragma offload target( mic : 0 ), wait( mSignal )
        {
        }

        // // finally called functions might also need the context, e.g. unbindTexture

    }

    setSynchronized();
}

void MICSyncToken::setSignal( int signal )
{
    LAMA_ASSERT_ERROR( !isSynchronized(), "cannot set signal as SyncToken is already synchronized" )

    LAMA_ASSERT_ERROR( mSignal < 0, "signal already set, cannot handle multiple signals" )

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

} // namespace

