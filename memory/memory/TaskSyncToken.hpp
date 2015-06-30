/**
 * @file TaskSyncToken.hpp
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
 * @brief Definition of a class that starts new Task.
 * @author Thomas Brandes
 * @date 04.05.2011
 */

#pragma once

// for dll_import
#include <common/config.hpp>

#include <tasking/Task.hpp>

// base classes
#include <memory/SyncToken.hpp>

// others

namespace memory
{

/** Class to run a function asynchronously as a asynchronous Task. */

class COMMON_DLL_IMPORTEXPORT TaskSyncToken: public SyncToken
{
public:

    /** This constructor starts a function as a Task and
     *  gives the SyncToken to wait on the completion.
     */
    TaskSyncToken( boost::function<void()> function, int numOmpThreads = 0 );

    /** This constructor starts a function as a Task and
     *  gives the SyncToken to wait on the completion.
     */
    TaskSyncToken();

    virtual void run( boost::function<void()> function, int numOmpThreads = 0 );

    virtual ~TaskSyncToken();

    /** Wait for the completion of the Task. */

    virtual void wait();

    /** Test for the completion of the Task. */

    virtual bool probe() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    boost::shared_ptr<tasking::Task> mTask;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

}

