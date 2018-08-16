/**
 * @file ThreadPool.cpp
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
 * @brief Implementation of class for a thread pool
 * @author Thomas Brandes
 * @date 27.12.2012
 */

// hpp
#include <scai/tasking/ThreadPool.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/OpenMP.hpp>

#include <memory>

namespace scai
{

using std::shared_ptr;
using std::function;

namespace tasking
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( ThreadPool::logger, "ThreadPool" )

/* ------------------------------------------------------------------------- */

shared_ptr<ThreadPoolTask> ThreadPoolTask::create(
    function<void()> work,
    unsigned int taskId,
    int numOmpThreads /* = 0 */ )
{
    shared_ptr<ThreadPoolTask> task = shared_ptr<ThreadPoolTask>( new ThreadPoolTask() );
    task->mWork = work;

    if ( numOmpThreads == 0 )
    {
        task->ompThreads = omp_get_max_threads();
    }
    else
    {
        task->ompThreads = numOmpThreads;
    }

    task->mState = DEFINED;
    task->mTaskId = taskId;
    task->mException = false;
    return task;
}

/* ------------------------------------------------------------------------- */

static void threadRoutine( ThreadPool::ThreadData& args )
{
    SCAI_LOG_THREAD( "ThreadPoolWorker_" << args.i )
    args.pool->worker( args.i );
}

ThreadPool::ThreadPool( int size )
{
    SCAI_LOG_INFO( logger, "Construct thread pool with " << size << " threads" )
    mMaxSize = size;
    mThreads.reset( new std::thread[ mMaxSize ] );
    mThreadArgs.reset( new ThreadData[ mMaxSize ] );
    mTaskId = 0; // Initialize counter for task ids

    // Create all threads just from the beginning, on demand might be possible later

    for ( int i = 0; i < mMaxSize; i++ )
    {
        mThreadArgs[i].pool = this;
        mThreadArgs[i].i    = i;
        mThreads[i] = std::thread( threadRoutine, std::ref( mThreadArgs[i] ) );
    }
}

/* ------------------------------------------------------------------------- */

shared_ptr<ThreadPoolTask> ThreadPool::schedule( std::function<void()> work, int numOmpThreads /* = 0 */ )
{
    SCAI_REGION( "ThreadPool::schedule" )
    std::thread::id thisThread = std::this_thread::get_id();
    bool isRecursiveTask = false;

    for ( int i = 0; i < mMaxSize; ++i )
    {
        if ( mThreads[i].get_id() == thisThread )
        {
            isRecursiveTask = true;
            break;
        }
    }

    shared_ptr<ThreadPoolTask> task = ThreadPoolTask::create( work, mTaskId++, numOmpThreads );

    if ( isRecursiveTask )
    {
        SCAI_LOG_WARN( logger, "Executing recursive taks synchronously to avoid deadlocks." )
        work();
        task->mState = ThreadPoolTask::FINISHED;
        return task;
    }

    std::unique_lock<std::mutex> lock( mTaskQueueMutex );
    mTaskQueue.push( task );
    SCAI_LOG_DEBUG( logger, "Added task " << task->mTaskId << " to task queue" )
    task->mState = ThreadPoolTask::QUEUED;
    //  notifiy one waiting worker
    mNotifyTask.notify_one();
    return task;
}

/* ------------------------------------------------------------------------- */

void ThreadPool::wait( shared_ptr<ThreadPoolTask> task )
{
    if ( !task )
    {
        COMMON_THROWEXCEPTION( "NULL pointer for task" )
    }

    SCAI_LOG_DEBUG( logger, "wait on task id = " << task->mTaskId << ", state = " << task->mState )

    while ( task->mState != ThreadPoolTask::FINISHED )
    {
        std::unique_lock<std::mutex> lock( mNotifyFinishMutex );

        if ( task->mState != ThreadPoolTask::FINISHED )
        {
            // wait on signal for a finishing thread
            // Attention: do not output here, as worker thread might finish and notify before wait
            mNotifyFinished.wait( lock );
        }

        SCAI_LOG_DEBUG( logger, "notified, task state = " << task->mState )
    }
}

/* ------------------------------------------------------------------------- */

void ThreadPool::worker( int id )
{
    // This method will be executed by all threads in the pool
    // Each thread picks up a task if available in the task queue
    // No busy wait, waits on signal mNotifyTask
    SCAI_LOG_INFO( logger, "worker thread " << id << " starts" )
    mWorkerState = WORKING;
    int ompThreads = -1;
    // wait for a new task

    while ( true )
    {
        shared_ptr<ThreadPoolTask> task;
        // pick up a new task if available
        {
            std::unique_lock<std::mutex> lock( mTaskQueueMutex );

            if ( mTaskQueue.empty() )
            {
                // Instead of busy wait this thread waits on notification
                SCAI_LOG_DEBUG( logger, "worker thread " << id << " waits on notify for new task" )
                mWorkerState = WAITING;
                mNotifyTask.wait( lock );
                mWorkerState = WORKING;
                SCAI_LOG_DEBUG( logger, "worker thread " << id << " notified about a new task" )
            }
            else
            {
                task = mTaskQueue.front();
                mTaskQueue.pop();

                if ( !task )
                {
                    SCAI_LOG_DEBUG( logger, "worker thread " << id << " picked shutdown task" )
                    // this is the shutdown task
                    break;
                }
            }
        }

        if ( task )
        {
            SCAI_LOG_DEBUG( logger,
                            "worker thread " << id << " runs task " << task->mTaskId << " with " << task->ompThreads << " OMP threads" )
            task->mState = ThreadPoolTask::RUNNING;

            if ( task->ompThreads != ompThreads )
            {
                omp_set_num_threads( task->ompThreads );
                ompThreads = task->ompThreads;
            }

            try
            {
                task->mWork();
            }
            catch ( common::Exception& ex )
            {
                SCAI_LOG_WARN( logger, "worker thread got exception, has been caught: " << ex.what() )
                task->mException = true;
            }
            catch ( ... )
            {
                SCAI_LOG_WARN( logger, "worker thread got exception, has been caught" )
                task->mException = true;
            }

            std::unique_lock<std::mutex> lock( mNotifyFinishMutex );
            task->mState = ThreadPoolTask::FINISHED;
            SCAI_LOG_DEBUG( logger, "worker thread " << id << " finished task " << task->mTaskId )
            // notify threads waiting on a finished task
            mNotifyFinished.notify_all();
        }
    }

    // worker is finished
    SCAI_LOG_INFO( logger, "worker thread " << id << " finishes" )
}

/* ------------------------------------------------------------------------- */

void ThreadPool::shutdown()
{
    SCAI_LOG_INFO( logger, "shut down " << mMaxSize << " threads, "
                   << mTaskQueue.size() << " tasks in queue" )
    shared_ptr<ThreadPoolTask> shutdownTask; // NULL pointer
    {
        // lock access to the task queue before adding shutdown tasks
        std::unique_lock<std::mutex> lock( mTaskQueueMutex );

        for ( int i = 0; i < mMaxSize; i++ )
        {
            mTaskQueue.push( shutdownTask );
        }

        // notifiy all waiting worker threads about new task
        mNotifyTask.notify_all();
    }
    SCAI_LOG_DEBUG( logger, "added " << mMaxSize << " shutdown tasks" )

    // and now wait for completion of all worker threads and delete them

    for ( int i = 0; i < mMaxSize; ++i )
    {
        SCAI_LOG_DEBUG( logger, "wait for worker thread " << i )
        mThreads[i].join();
        SCAI_LOG_DEBUG( logger, "worker thread " << i << " terminated (joined)" )
    }
}

/* ------------------------------------------------------------------------- */

ThreadPool::~ThreadPool()
{
    shutdown();
}

} /* end namespace tasking */

} /* end namespace scai */
