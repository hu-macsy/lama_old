/**
 * @file ThreadPool.cpp
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

namespace scai
{

using common::Thread;
using common::shared_ptr;
using common::function;

namespace tasking
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( ThreadPool::logger, "ThreadPool" )

/* ------------------------------------------------------------------------- */

shared_ptr<ThreadTask> ThreadTask::create(
    function<void()> work,
    unsigned int taskId,
    int numOmpThreads /* = 0 */ )
{
    shared_ptr<ThreadTask> task = shared_ptr<ThreadTask>( new ThreadTask() );

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

    mThreads.reset( new common::Thread[ mMaxSize ] );
    mThreadArgs.reset( new ThreadData[ mMaxSize ] );

    mTaskId = 0; // Initialize counter for task ids

    // Create all threads just from the beginning, on demand might be possible later

    for ( int i = 0; i < mMaxSize; i++ )
    {
        mThreadArgs[i].pool = this;
        mThreadArgs[i].i    = i;

        mThreads[i].run( &threadRoutine, mThreadArgs[i] );
    }
}

/* ------------------------------------------------------------------------- */

shared_ptr<ThreadTask> ThreadPool::schedule( function<void()> work, int numOmpThreads /* = 0 */ )
{
    SCAI_REGION( "ThreadPool::schedule" )
    Thread::Id thisThread = Thread::getSelf();
    bool isRecursiveTask = false;

    for ( int i = 0; i < mMaxSize; ++i )
    {
        if ( mThreads[i].getId() == thisThread )
        {
            isRecursiveTask = true;
            break;
        }
    }

    shared_ptr<ThreadTask> task = ThreadTask::create( work, mTaskId++, numOmpThreads );

    if ( isRecursiveTask )
    {
        SCAI_LOG_WARN( logger, "Executing recursive taks synchronously to avoid deadlocks." )
        work();
        task->mState = ThreadTask::FINISHED;
        return task;
    }

    Thread::ScopedLock lock( mTaskQueueMutex );

    mTaskQueue.push( task );
    SCAI_LOG_DEBUG( logger, "Added task " << task->mTaskId << " to task queue" )
    task->mState = ThreadTask::QUEUED;

    //  notifiy one waiting worker

    mNotifyTask.notifyOne();

    return task;
}

/* ------------------------------------------------------------------------- */

void ThreadPool::wait( shared_ptr<ThreadTask> task )
{
    if ( !task )
    {
        COMMON_THROWEXCEPTION( "NULL pointer for task" )
    }

    SCAI_LOG_DEBUG( logger, "wait on task id = " << task->mTaskId << ", state = " << task->mState )

    while ( task->mState != ThreadTask::FINISHED )
    {
        Thread::ScopedLock lock( mNotifyFinishMutex );

        if ( task->mState != ThreadTask::FINISHED )
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
        shared_ptr<ThreadTask> task;

        // pick up a new task if available

        {
            Thread::ScopedLock lock( mTaskQueueMutex );

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

            task->mState = ThreadTask::RUNNING;

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
                SCAI_LOG_INFO( logger, "worker thread got exception, has been caught: " << ex.what() )
                task->mException = true;
            }
            catch ( ... )
            {
                SCAI_LOG_INFO( logger, "worker thread got exception, has been caught" )
                task->mException = true;
            }

            Thread::ScopedLock lock( mNotifyFinishMutex );

            task->mState = ThreadTask::FINISHED;

            SCAI_LOG_DEBUG( logger, "worker thread " << id << " finished task " << task->mTaskId )

            // notify threads waiting on a finished task

            mNotifyFinished.notifyAll();
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

    shared_ptr<ThreadTask> shutdownTask; // NULL pointer

    {
        // lock access to the task queue before adding shutdown tasks

        Thread::ScopedLock lock( mTaskQueueMutex );

        for ( int i = 0; i < mMaxSize; i++ )
        {
            mTaskQueue.push( shutdownTask );
        }

        // notifiy all waiting worker threads about new task

        mNotifyTask.notifyAll();
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
