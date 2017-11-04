/**
 * @file Task.cpp
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
 * @brief Definition of the class Task for asynchronous function execution.
 * @author Thomas Brandes
 * @date 21.11.2011
 */

// hpp
#include <scai/tasking/Task.hpp>

// internal scai libraries
#include <scai/common/macros/throw.hpp>
#include <scai/common/Settings.hpp>

// std
#include <limits>
#include <functional>

namespace scai
{

namespace tasking
{

/* ------------------------------------------------------------------------- */

std::shared_ptr<ThreadPool> Task::theThreadPool;

/* ------------------------------------------------------------------------- */

ThreadPool& Task::getThreadPool()
{
    int poolSize = 1;
    common::Settings::getEnvironment( poolSize, "SCAI_THREADPOOL_SIZE" );

    if ( !theThreadPool )
    {
        theThreadPool.reset( new ThreadPool( poolSize ) );
    }

    return *theThreadPool;
}

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( Task::logger, "Task" )

/* ------------------------------------------------------------------------- */

Task::Task( std::function<void()> taskFunction, int numOmpThreads /* = 0 */ )

{
    SCAI_LOG_DEBUG( logger, "Creating Task" )
    mTask = getThreadPool().schedule( taskFunction, numOmpThreads );
    SCAI_LOG_DEBUG( logger, "Task created" )
}

/* ------------------------------------------------------------------------- */

Task::~Task()
{
    SCAI_LOG_INFO( logger, "~Task" )
    synchronize();
}

/* ------------------------------------------------------------------------- */

void Task::synchronize()
{
    if ( !mTask )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger, "Waiting for task = " << mTask->mTaskId << ", state = " << mTask->mState )
    getThreadPool().wait( mTask );
    SCAI_LOG_DEBUG( logger, "Task = " << mTask->mTaskId << " finished" << ", exception = " << mTask->mException )

    if ( mTask->mException )
    {
        COMMON_THROWEXCEPTION( "Task got exception" )
    }

    mTask.reset();
}

bool Task::probe() const
{
    return ThreadPoolTask::FINISHED == mTask->mState;
}

/* ------------------------------------------------------------------------- */

void Task::writeAt( std::ostream& stream ) const
{
    stream << "Task( id = " << getId() << " )";
}

/* ------------------------------------------------------------------------- */

Task::Id Task::getId() const
{
    if ( mTask )
    {
        return mTask->mTaskId;
    }
    else
    {
        return std::numeric_limits<Task::Id>::max();
    }
}

} /* end namespace tasking */

} /* end namespace scai */
