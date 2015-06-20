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
 * @since 1.0.0
 */

// hpp
#include <tasking/ThreadPool.hpp>

// tracing
#include <tracing/tracing.hpp>

// assert
#include <common/Exception.hpp>

// boost
#include <boost/bind.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_set_num_threads( x )
#endif

using common::Thread;

namespace tasking
{

/* ------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( ThreadPool::logger, "ThreadPool" )

/* ------------------------------------------------------------------------- */

boost::shared_ptr<ThreadTask> ThreadTask::create(
    boost::function<void()> work,
    unsigned int taskId,
    int numOmpThreads /* = 0 */)
{
    boost::shared_ptr<ThreadTask> task = boost::shared_ptr<ThreadTask>( new ThreadTask() );

    task->mWork = work;

    if( numOmpThreads == 0 )
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

static void* threadRoutine( void* p )
{
    ThreadPool::ThreadData* args = ( ThreadPool::ThreadData* ) p;

    LAMA_LOG_THREAD( "ThreadPoolWorker_" << args->i )

    args->pool->worker( args->i );
}

ThreadPool::ThreadPool( int size )
{
    LAMA_LOG_INFO( logger, "Construct thread pool with " << size << " threads" )
    mMaxSize = size;
    mThreads.reserve( mMaxSize );
    mThreadArgs.reserve( mMaxSize );

    mTaskId = 0; // Initialize counter for task ids

    // Create all threads just from the beginning, on demand might be possible later

    for( int i = 0; i < mMaxSize; i++ )
    {
        pthread_t id;
        
        mThreadArgs[i].pool = this;
        mThreadArgs[i].i    = i;

        int rc = pthread_create( &id, NULL, &threadRoutine, &mThreadArgs[i] );

        mThreads.push_back( id );
    }
}

/* ------------------------------------------------------------------------- */

boost::shared_ptr<ThreadTask> ThreadPool::schedule( boost::function<void()> work, int numOmpThreads /* = 0 */)
{
    LAMA_REGION( "ThreadPool::schedule" )
    Thread::Id thisThread = Thread::getSelf();
    bool isRecursiveTask = false;
    std::vector<Thread::Id>::const_iterator end = mThreads.end();

    for( std::vector<pthread_t>::const_iterator it = mThreads.begin(); it != end; ++it )
    {
        if( *it == thisThread )
        {
            isRecursiveTask = true;
            break;
        }
    }

    boost::shared_ptr<ThreadTask> task = ThreadTask::create( work, mTaskId++, numOmpThreads );

    if ( isRecursiveTask )
    {
        LAMA_LOG_WARN( logger, "Executing recursive taks synchronously to avoid deadlocks." )
        work();
        task->mState = ThreadTask::FINISHED;
        return task;
    }

    Thread::ScopedLock lock( mTaskQueueMutex );

    mTaskQueue.push( task );
    LAMA_LOG_DEBUG( logger, "Added task " << task->mTaskId << " to task queue" )
    task->mState = ThreadTask::QUEUED;

    //  notifiy one waiting worker

    mNotifyTask.notify_one();

    return task;
}

/* ------------------------------------------------------------------------- */

void ThreadPool::wait( boost::shared_ptr<ThreadTask> task )
{
    if ( !task )
    {
        COMMON_THROWEXCEPTION( "NULL pointer for task" )
    }

    LAMA_LOG_DEBUG( logger, "wait on task id = " << task->mTaskId << ", state = " << task->mState )

    while( task->mState != ThreadTask::FINISHED )
    {
        Thread::ScopedLock lock( mNotifyFinishMutex );

        if( task->mState != ThreadTask::FINISHED )
        {
            // wait on signal for a finishing thread
            // Attention: do not output here, as worker thread might finish and notify before wait

            mNotifyFinished.wait( lock );
        }

        LAMA_LOG_DEBUG( logger, "notified, task state = " << task->mState )
    }
}

/* ------------------------------------------------------------------------- */

void ThreadPool::worker( int id )
{
    // This method will be executed by all threads in the pool
    // Each thread picks up a task if available in the task queue
    // No busy wait, waits on signal mNotifyTask

    LAMA_LOG_INFO( logger, "worker thread " << id << " starts" )

    mWorkerState = WORKING;

    int ompThreads = -1;
    // wait for a new task

    while( true )
    {
        boost::shared_ptr<ThreadTask> task;

        // pick up a new task if available

        {
            Thread::ScopedLock lock( mTaskQueueMutex );

            if( mTaskQueue.empty() )
            {
                // Instead of busy wait this thread waits on notification

                LAMA_LOG_DEBUG( logger, "worker thread " << id << " waits on notify for new task" )

                mWorkerState = WAITING;

                mNotifyTask.wait( lock );

                mWorkerState = WORKING;

                LAMA_LOG_DEBUG( logger, "worker thread " << id << " notified about a new task" )
            }
            else
            {
                task = mTaskQueue.front();
                mTaskQueue.pop();

                if( !task )
                {
                    LAMA_LOG_DEBUG( logger, "worker thread " << id << " picked shutdown task" )
                    // this is the shutdown task
                    break;
                }
            }
        }

        if( task )
        {
            LAMA_LOG_DEBUG( logger,
                            "worker thread " << id << " runs task " << task->mTaskId << " with " << task->ompThreads << " OMP threads" )

            task->mState = ThreadTask::RUNNING;

            if( task->ompThreads != ompThreads )
            {
                omp_set_num_threads( task->ompThreads );
                ompThreads = task->ompThreads;
            }

            try
            {
                task->mWork();
            }
            catch( ... )
            {
                LAMA_LOG_INFO( logger, "worker thread got exception, has been caught" )
                task->mException = true;
            }

            Thread::ScopedLock lock( mNotifyFinishMutex );

            task->mState = ThreadTask::FINISHED;

            LAMA_LOG_DEBUG( logger, "worker thread " << id << " finished task " << task->mTaskId )

            // notify threads waiting on a finished task

            mNotifyFinished.notify_all();
        }
    }

    // worker is finished

    LAMA_LOG_INFO( logger, "worker thread " << id << " finishes" )
}

/* ------------------------------------------------------------------------- */

ThreadPool::~ThreadPool()
{
    LAMA_LOG_INFO( logger, "~ThreadPool: shut down " << mThreads.size() << " threads" )

    boost::shared_ptr<ThreadTask> shutdownTask; // NULL pointer

    {
        // lock access to the task queue before adding shutdown tasks

        Thread::ScopedLock lock( mTaskQueueMutex );

        for( size_t i = 0; i < mThreads.size(); i++ )
        {
            mTaskQueue.push( shutdownTask );
        }

        // notifiy all waiting worker threads about new task

        mNotifyTask.notify_all();
    }

    LAMA_LOG_DEBUG( logger, "added " << mThreads.size() << " shutdown tasks" )

    // and now wait for completion of all worker threads and delete them

    for( size_t i = 0; i < mThreads.size(); ++i )
    {
        LAMA_LOG_DEBUG( logger, "wait for worker thread " << i )
        pthread_join( mThreads[i], NULL );
        LAMA_LOG_DEBUG( logger, "worker thread " << i << " terminated (joined)" )
    }
}

} // namespace lama
