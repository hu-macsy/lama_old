/**
 * @file SyncToken.hpp
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
 * @brief Definition of a base class for synchronization of computations and communications.
 * @author Thomas Brandes, Jiri Kraus
 * @date 22.03.2011
 */

#pragma once

// for dll_import
#include <common/config.hpp>

// base classes
#include <common/NonCopyable.hpp>
#include <common/Printable.hpp>

// logging
#include <logging/logging.hpp>

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

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
 *    auto_ptr<SyncToken> token = new XXXSyncToken( ... )
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

class COMMON_DLL_IMPORTEXPORT SyncToken: public Printable, private common::NonCopyable
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

    void pushToken( boost::shared_ptr<SyncTokenMember> member );

    /** Add a routine to be called after synchronization. */

    void pushRoutine( boost::function<void()> function );

protected:

    /** Default constructor can only be called by derived classes. */

    SyncToken();

    /** This method should be called by base classes after a successful wait. */

    void setSynchronized();

    /** Logger for this class. */

    LAMA_LOG_DECL_STATIC_LOGGER(logger)

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

    /** Vector of shared pointers  that will be released after completion. */

    std::vector< boost::shared_ptr<SyncTokenMember > > mTokens;

    bool mSynchronized;  //!< if true the token has already been synchronized.

    std::vector< boost::function<void()> > mSynchronizedFunctions;
};

}