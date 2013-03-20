/**
 * @file TaskSyncToken.cpp
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
 * @brief TaskSyncToken.cpp
 * @author Thomas Brandes
 * @date 04.05.2011
 * $Id$
 */

// hpp
#include <lama/task/TaskSyncToken.hpp>

// tracing
#include <lama/tracing.hpp>

// boost
#include <boost/thread/xtime.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/utility.hpp>

using namespace boost;

namespace lama
{

LAMA_LOG_DEF_LOGGER( TaskSyncToken::logger, "SyncToken.TaskSyncToken" );

TaskSyncToken::TaskSyncToken( boost::function<void()> function, int numOmpThreads /* = 0 */)
    : mTask( new Task( function, numOmpThreads ) )
{
    LAMA_LOG_DEBUG( logger, "Thread " << *mTask << " with function " << function << " started." );
}

TaskSyncToken::TaskSyncToken()
{
    // empty task token
}

void TaskSyncToken::run( boost::function<void()> function, int numOmpThreads /* = 0 */)
{
    mTask = boost::shared_ptr<Task>( new Task( function, numOmpThreads ) );

    LAMA_LOG_DEBUG( logger, "Thread " << *mTask << " with function " << function << " started." );
}

TaskSyncToken::~TaskSyncToken()
{
    wait();
}

void TaskSyncToken::wait()
{
    LAMA_REGION( "TaskSyncToken.wait" );

    if( isSynchronized() )
    {
        return;
    }

    // set synchronized already here as synchronize might throw exception

    setSynchronized();

    if( mTask ) // might be running task
    {
        LAMA_LOG_DEBUG( logger, "Waiting for thread " << mTask );

        mTask->synchronize(); // not really needed as destructor will be called

        mTask = boost::shared_ptr<Task>(); // implies call destructor of mTask
    }
}

bool TaskSyncToken::probe() const
{
    if( !mTask )
    {
        return true;
    }

    return mTask->probe();
}

void TaskSyncToken::writeAt( std::ostream& stream ) const
{
    stream << "TaskSyncToken( ";
    if( mTask )
    {
        stream << *mTask;
    }
    else
    {
        stream << "NULL";
    }
    stream << ", synchronized = " << isSynchronized() << ")";
}

}
