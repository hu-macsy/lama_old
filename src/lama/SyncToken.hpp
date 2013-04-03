/**
 * @file SyncToken.hpp
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
 * @brief Definition of a base class for synchronization of computations and communications.
 * @author Thomas Brandes, Jiri Kraus
 * @date 22.03.2011
 * $Id$
 */

#ifndef LAMA_SYNC_TOKEN_HPP_
#define LAMA_SYNC_TOKEN_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>
#include <lama/Printable.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/function.hpp>

#include <vector>
#include <memory>

namespace lama
{

class BaseAccess;
class _LAMAArray;

/** Abstract class that defines tokens for asynchronous operations.
 *
 * Communication and computations can be executed asynchronously. A
 * token is needed mainly to wait on the completion of the operation.
 *
 * This class also supports the possibility to push LAMA array accesses and LAMA arrays
 * to a token. After successful synchronization, the accesses are destroyed and
 * the arrays can be accesses for other purposes. LAMA arrays used for temporary purpose
 * can also be freed at synchronization.
 */

class LAMA_DLL_IMPORTEXPORT SyncToken: public Printable, private NonCopyable
{
public:

    /** Default constructor */

    SyncToken();

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

    /** Add a read/write access to the token so that LAMA arrays will be locked until synchronization. */

    void pushAccess( std::auto_ptr<BaseAccess> access );

    /** Add a LAMA array that will be free after synchronization */

    void pushArray( std::auto_ptr<_LAMAArray> array );

    /** Add a Synctoken that will be synchronized after synchronization */

    void pushSyncToken( std::auto_ptr<SyncToken> syncToken );

protected:

    /** This method should be called by base classes after a successful wait. */

    void setSynchronized();

    LAMA_LOG_DECL_STATIC_LOGGER(logger) //!< logger for this class

private:

    /** Vector of accesses that will be freed after completion. */

    std::vector< BaseAccess* > mAccesses;

    std::vector< _LAMAArray* > mArrays;

    std::vector< SyncToken* > mChilds;

    bool mSynchronized;
};

}

#endif // LAMA_SYNC_TOKEN_HPP_
