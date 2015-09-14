/**
 * @file TaskSyncToken.cpp
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
 * @brief Implementation of methods for class TaskSyncToken.
 * @author Thomas Brandes
 * @date 14.07.2015
 */

// hpp
#include <scai/tasking/TaskSyncToken.hpp>

// internal scai libraries
#include <scai/common/bind.hpp>
#include <scai/common/exception/Exception.hpp>

using scai::tasking::Task;
using scai::common::function;
using scai::common::shared_ptr;

namespace scai
{

namespace tasking
{

SCAI_LOG_DEF_LOGGER( TaskSyncToken::logger, "SyncToken.TaskSyncToken" )

TaskSyncToken::TaskSyncToken( function<void()> routine, int numOmpThreads ) :

    mTask( new Task( routine, numOmpThreads ) )

{
    SCAI_LOG_DEBUG( logger, "Thread " << *mTask << " with routine started." )
}

TaskSyncToken::TaskSyncToken()
{
    SCAI_LOG_DEBUG( logger, "TaskSyncToken(): no function set" )
    // empty task token
}

void TaskSyncToken::run( function<void()> routine, int numOmpThreads /* = 0 */ )
{
    mTask = shared_ptr<Task>( new Task( routine, numOmpThreads ) );

    SCAI_LOG_DEBUG( logger, "Thread " << *mTask << " with routine started." )
}

TaskSyncToken::~TaskSyncToken()
{
    SCAI_LOG_DEBUG( logger, "~TaskSyncToken: wait" )

    wait();

    SCAI_LOG_DEBUG( logger, "~TaskSyncToken: wait done" )
}

void TaskSyncToken::wait()
{
    SCAI_LOG_DEBUG( logger, "wait" )

    if ( isSynchronized() )
    {
        return;
    }

    if ( mTask ) // might be running task
    {
        SCAI_LOG_DEBUG( logger, "Waiting for thread " << mTask )

        try
        {
            mTask->synchronize(); // not really needed as destructor will be called
        }
        catch ( common::Exception& ex )
        {
            SCAI_LOG_ERROR( logger, "Task caught exception: " << ex.what() );
        }

        mTask = shared_ptr<Task>(); // implies call destructor of mTask

        mTask.reset();  // implies call destructor of mTask
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "no task to synchronize." )
    }

    // Now we can free tokens and call clean functions, do not call it before mTask is synchronized

    setSynchronized();
}

bool TaskSyncToken::probe() const
{
    if ( !mTask )
    {
        return true;
    }

    return mTask->probe();
}

void TaskSyncToken::writeAt( std::ostream& stream ) const
{
    stream << "TaskSyncToken( ";

    if ( mTask )
    {
        stream << *mTask;
    }
    else
    {
        stream << "NULL";
    }

    stream << ", synchronized = " << isSynchronized() << ")";
}

} /* end namespace tasking */

} /* end namespace scai */
