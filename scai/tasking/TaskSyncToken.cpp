/**
 * @file TaskSyncToken.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Implementation of methods for class TaskSyncToken.
 * @author Thomas Brandes
 * @date 14.07.2015
 */

// internal scai libraries

#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/common/macros/throw.hpp>

#include <memory>
#include <functional>

using std::shared_ptr;

namespace scai
{

using tasking::Task;

namespace tasking
{

SCAI_LOG_DEF_LOGGER( TaskSyncToken::logger, "SyncToken.TaskSyncToken" )

TaskSyncToken::TaskSyncToken( std::function<void()> routine, int numOmpThreads ) :

    mTask( new Task( routine, numOmpThreads ) )

{
    SCAI_LOG_DEBUG( logger, "Thread " << *mTask << " with routine started." )
}

TaskSyncToken::TaskSyncToken()
{
    SCAI_LOG_DEBUG( logger, "TaskSyncToken(): no function set" )
    // empty task token
}

void TaskSyncToken::run( std::function<void()> routine, int numOmpThreads /* = 0 */ )
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

TaskSyncToken* TaskSyncToken::getCurrentSyncToken()
{
    SyncToken* syncToken = SyncToken::getCurrentSyncToken();

    if ( syncToken == NULL )
    {
        return NULL;
    }

    // make a dynamic CAST
    TaskSyncToken* taskSyncToken = dynamic_cast<TaskSyncToken*>( syncToken );

    // If the current sync token is not a Task token it is very likely an error

    if ( taskSyncToken == NULL )
    {
        SCAI_LOG_ERROR( logger, "Current sync token = " << *syncToken << " not TaskSyncToken as expected" )
    }

    // But might not be too serious so probably NULL results in synchronous execution
    return taskSyncToken;
}

} /* end namespace tasking */

} /* end namespace scai */
