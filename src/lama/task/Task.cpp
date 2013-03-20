/**
 * @file Task.cpp
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
 * @brief Definition of the class Task for asynchronous function execution.
 * @author Thomas Brandes
 * @date 21.11.2011
 * $Id$
 */

// hpp
#include <lama/task/Task.hpp>

// others
#include <lama/exception/Exception.hpp>

// boost
#include <boost/bind.hpp>

#include <limits>

namespace lama
{

/* ------------------------------------------------------------------------- */

boost::shared_ptr<LAMAThreadPool> Task::theThreadPool;

/* ------------------------------------------------------------------------- */

LAMAThreadPool& Task::getThreadPool()
{
    if( !theThreadPool )
    {
        theThreadPool.reset( new LAMAThreadPool( 1 ) );
    }

    return *theThreadPool;
}

/* ------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( Task::logger, "Task" );

/* ------------------------------------------------------------------------- */

Task::Task( boost::function<void()> taskFunction, int numOmpThreads /* = 0 */)

{
    LAMA_LOG_DEBUG( logger, "Creating Task" );

    mTask = getThreadPool().schedule( taskFunction, numOmpThreads );

    LAMA_LOG_DEBUG( logger, "Task created" );
}

/* ------------------------------------------------------------------------- */

Task::~Task()
{
    LAMA_LOG_INFO( logger, "~Task" );
}

/* ------------------------------------------------------------------------- */

void Task::synchronize()
{
    if( !mTask )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "Waiting for task = " << mTask->mTaskId << ", state = " << mTask->mState );

    getThreadPool().wait( mTask );

    LAMA_LOG_DEBUG( logger, "Task = " << mTask->mTaskId << " finished" << ", exception = " << mTask->mException );

    if( mTask->mException )
    {
        LAMA_THROWEXCEPTION( "Task got exception" );
    }

    mTask.reset();
}

bool Task::probe() const
{
    return LAMAThreadTask::FINISHED == mTask->mState;
}

/* ------------------------------------------------------------------------- */

void Task::writeAt( std::ostream& stream ) const
{
    stream << "Task( id = " << getId() << " )";
}

/* ------------------------------------------------------------------------- */

Task::Id Task::getId() const
{
    if( mTask )
    {
        return mTask->mTaskId;
    }
    else
    {
        return std::numeric_limits<Task::Id>::max();
    }
}

} // namespace lama
