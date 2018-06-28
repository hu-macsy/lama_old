/**
 * @file SyncToken.hpp
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
 * @brief Definition of a base class for synchronization of computations and communications.
 * @author Thomas Brandes, Jiri Kraus
 * @date 22.03.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>
#include <scai/common/Printable.hpp>

#include <scai/logging.hpp>

// std
#include <vector>
#include <memory>
#include <functional>

namespace scai
{

/** Namespace for all data structures used in the tasking library. */

namespace tasking
{

/** Simple base class where classes might derived from to become SyncTokenMember.
 *
 *  A SyncToken can take over ownership of shared pointer objects derived from this class.
 */

class SyncTokenMember
{
public:

    SyncTokenMember()
    {
    }

    virtual ~SyncTokenMember()
    {
    }
};

/** Abstract class that defines tokens for asynchronous operations.
 *
 * Communication and computations can be executed asynchronously. A
 * token is needed mainly to wait on the completion of the operation.
 *
 * This class also supports the possibility to push LAMA array accesses and LAMA arrays
 * to a token. After successful synchronization, the accesses/arrays are release and
 * the arrays can be accesses for other purposes.
 *
 * All started asynchronous operations in LAMA must be synchronized. This is
 * absolutely mandatory and can be done in the following ways:
 *
 * \code
 *    std::unique_ptr<SyncToken> token ( new XXXSyncToken( ...) )
 *    ! synchronization is alway done when object will be deleted at the end of the scope
 *
 *    token->wait();     // explicit wait
 *
 *    !  This is not recommened but works
 *    SyncToken* token = new XXXSyncToken( ... )
 *       ....
 *    delete token;   // easy to forget, token will never be synchronized
 * \endcode
 *
 */

class COMMON_DLL_IMPORTEXPORT SyncToken: public common::Printable, private common::NonCopyable
{
public:

    /** Destructor.
     *
     *  The destructor will not wait here for completion of the action
     *  this must be done in each derived class.
     */

    virtual ~SyncToken();

    /** Query whether token has already been synchronized. */

    bool isSynchronized() const;

    /** Method to wait on the completion of an operation. */

    virtual void wait() = 0;

    /** Predicate to ask if an asynchronous operation is already completed. */

    virtual bool probe() const = 0;

    /** Base class provides a default implementation for the virtual method of Printable.
     *
     * @see Printable for more details.
     */
    virtual void writeAt( std::ostream& stream ) const;

    /** Add a shared pointer to this SyncToken to take ownership
     *
     *  @param member shared pointer to an object that can be SyncTokenMember
     */

    void pushToken( std::shared_ptr<SyncTokenMember> member );

    /** Add a routine to be called after synchronization. */

    void pushRoutine( std::function<void()> routine );

    /**
     *  Set this SyncToken as the current one of this thread.
     *
     *  Only one SyncToken can be the current one.
     *
     *  Global access to the current sync token makes design easier
     *  as it can be decided locally where and how to start an
     *  asynchronous operation.
     *
     *  Be careful: SyncToken is not current for the executing thread
     *              but for the parent thread that will run it.
     */

    void setCurrent();

    /**
     *  Current SyncToken will be no more current one.
     */
    void unsetCurrent();

    /**
     *  Get the current sync token of this thread.
     */
    static SyncToken* getCurrentSyncToken();

protected:

    /** Default constructor can only be called by derived classes. */

    SyncToken();

    /** This method should be called by base classes after a successful wait. */

    void setSynchronized();

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /** Helper class for a static object on which should be act at termination. */

    class CGuard
    {
    public:

        /** Constructor. */

        CGuard();

        /** Destructor. */

        ~CGuard();
    };

    // counts allocated - freed sync tokens, verify for 0 at then end

    static int countSyncToken;

    static CGuard cguard;//!< required to call routine at its destructor

    /** Each thread can set globally (thread-private) a SyncToken */

    static thread_local SyncToken* currentSyncToken;

    /** Vector of shared pointers  that will be released after completion. */

    std::vector< std::shared_ptr<SyncTokenMember > > mTokens;

    bool mSynchronized;  //!< if true the token has already been synchronized.

    std::vector< std::function<void()> > mSynchronizedFunctions;

public:

    class ScopedAsynchronous
    {
    public:

        ScopedAsynchronous( SyncToken& token ) : mToken( &token )
        {
            mToken->setCurrent();
        }

        ScopedAsynchronous( SyncToken* token ) : mToken( token )
        {
            if ( mToken != NULL )
            {
                mToken->setCurrent();
            }
        }

        ~ScopedAsynchronous()
        {
            if ( mToken != NULL )
            {
                mToken->unsetCurrent();
            }
        }

    private:

        SyncToken* mToken;
    };
};

} /* end namespace tasking */

} /* end namespace scai */

#define SCAI_ASYNCHRONOUS( token ) scai::tasking::SyncToken::ScopedAsynchronous _SCAIAsyncScope( token );
