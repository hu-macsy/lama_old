/**
 * @file SyncToken.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

// hpp
#include <lama/SyncToken.hpp>

// others
#include <lama/BaseAccess.hpp>
#include <lama/LAMAArray.hpp>

#include <lama/exception/LAMAAssert.hpp>
#include <lama/exception/Exception.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( SyncToken::logger, "SyncToken" )

SyncToken::SyncToken()
    : mSynchronized( false )
{
    LAMA_LOG_DEBUG( logger, "SyncToken constructed" )
}

SyncToken::~SyncToken()
{
    LAMA_LOG_DEBUG( logger, "~SyncToken" )

    if ( !mSynchronized )
    {
        LAMA_LOG_WARN( logger, "no synchronization called on SyncToken" )
    }

    LAMA_LOG_DEBUG( logger, "still delete " << mAccesses.size() << " accesses and " << mArrays.size() << " arrays" )

    while ( !mAccesses.empty() )
    {
        delete mAccesses.back();
        mAccesses.pop_back();
    }

    while ( !mArrays.empty() )
    {
        delete mArrays.back();
        mArrays.pop_back();
    }

    while ( !mChilds.empty() )
    {
        mChilds.back()->wait();
        delete mChilds.back();
        mChilds.pop_back();
    }
}

void SyncToken::writeAt( std::ostream& stream ) const
{
    stream << "SyncToken( synchronized = " << mSynchronized << " )";
}

void SyncToken::pushAccess( std::auto_ptr<BaseAccess> access )
{
    LAMA_ASSERT_ERROR( access.get(), "NULL access cannot be pushed for synchronization." )

    if ( mSynchronized )
    {
        LAMA_LOG_DEBUG( logger, *this << ": push access not done, already synchronized" )

        // delete the access, do not push it anymore

        access.reset();
    }
    else
    {
        LAMA_LOG_DEBUG( logger, *this << ": push access, will be freed at synchronization" )

        // take ownership of the pointer, have to be deleted at synchronization

        mAccesses.push_back( access.release() );
    }
}

void SyncToken::pushArray( std::auto_ptr<_LAMAArray> array )
{
    LAMA_ASSERT_ERROR( array.get(), "NULL array cannot be pushed for synchronization." )

    if ( mSynchronized )
    {
        LAMA_LOG_DEBUG( logger, *this << ": push array not done, already synchronized" )

        // delete the pointer, do not push it anymore

        array.reset();
    }
    else
    {
        LAMA_LOG_DEBUG( logger, *this << ": push array, will be freed at synchronization" )

        // take ownership of the pointer, have to be deleted at synchronization

        mArrays.push_back( array.release() );
    }
}

void SyncToken::pushSyncToken( std::auto_ptr<SyncToken> syncToken )
{
    LAMA_ASSERT_ERROR( syncToken.get(), "NULL SyncToken cannot be pushed for synchronization." )

    if ( mSynchronized )
    {
        LAMA_LOG_DEBUG( logger, *this << ": push SyncToken not done, already synchronized" )

        // delete the pointer, do not push it anymore

        syncToken.reset();
    }
    else
    {
        LAMA_LOG_DEBUG( logger, *this << ": push SyncToken, will be synchronized after synchronization" )

        // take ownership of the pointer, have to be deleted at synchronization

        mChilds.push_back( syncToken.release() );
    }

}

bool SyncToken::isSynchronized() const
{
    return mSynchronized;
}

void SyncToken::setSynchronized()
{
    if ( mSynchronized )
    {
        LAMA_THROWEXCEPTION( *this << " is already synchronized" )
    }

    mSynchronized = true;

    // after synchronization we can unlock used LAMA arrays

    LAMA_LOG_DEBUG( logger,
                    *this << ": delete " << mAccesses.size() << " accesses, " << mArrays.size() << " arrays and synchronizing " << mChilds.size() << " childs " )

    while ( !mAccesses.empty() )
    {
        BaseAccess* lastAccess = mAccesses.back();
        LAMA_LOG_INFO( logger, "delete " << lastAccess )

        delete lastAccess;
        mAccesses.pop_back();
    }

    while ( !mArrays.empty() )
    {
        delete mArrays.back();
        mArrays.pop_back();
    }

    while ( !mChilds.empty() )
    {
        mChilds.back()->wait();
        delete mChilds.back();
        mChilds.pop_back();
    }
}

}
