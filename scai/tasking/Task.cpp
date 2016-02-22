/**
 * @file Task.cpp
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
 * @brief Definition of the class Task for asynchronous function execution.
 * @author Thomas Brandes
 * @date 21.11.2011
 */

// hpp
#include <scai/tasking/Task.hpp>

// internal scai libraries
#include <scai/common/macros/throw.hpp>
#include <scai/common/bind.hpp>

// std
#include <limits>

namespace scai
{

namespace tasking
{

/* ------------------------------------------------------------------------- */

common::shared_ptr<ThreadPool> Task::theThreadPool;

/* ------------------------------------------------------------------------- */

ThreadPool& Task::getThreadPool()
{
    if ( !theThreadPool )
    {
        theThreadPool.reset( new ThreadPool( 1 ) );
    }

    return *theThreadPool;
}

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( Task::logger, "Task" )

/* ------------------------------------------------------------------------- */

Task::Task( common::function<void()> taskFunction, int numOmpThreads /* = 0 */ )

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
    return ThreadTask::FINISHED == mTask->mState;
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
