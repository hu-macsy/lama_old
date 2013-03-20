/**
 * @file LAMAThreadPool.cpp
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
 * @brief Implementation of class for a thread pool
 * @author Thomas Brandes
 * @date 27.12.2012
 * $Id$
 */

// hpp
#include <lama/task/LAMAThreadPool.hpp>

// tracing
#include <lama/tracing.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// boost
#include <boost/bind.hpp>

#include <omp.h>

namespace lama
{

/* ------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( LAMAThreadPool::logger, "LAMAThreadPool" );

/* ------------------------------------------------------------------------- */

boost::shared_ptr<LAMAThreadTask> LAMAThreadTask::create(
    boost::function<void()> work,
    unsigned int taskId,
    int numOmpThreads /* = 0 */)
{
    boost::shared_ptr<LAMAThreadTask> task = boost::shared_ptr<LAMAThreadTask>( new LAMAThreadTask() );

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

LAMAThreadPool::LAMAThreadPool( int size )
{
    LAMA_LOG_INFO( logger, "Construct thread pool with " << size << " threads" );
    mMaxSize = size;
    mThreads.reserve( mMaxSize );

    mTaskId = 0; // Initialize counter for task ids

    // Create all threads just from the beginning, on demand might be possible later

    for( int i = 0; i < mMaxSize; i++ )
    {
        mThreads.push_back( new boost::thread( boost::bind( &LAMAThreadPool::worker, this, i ) ) );
    }
}

/* ------------------------------------------------------------------------- */

boost::shared_ptr<LAMAThreadTask> LAMAThreadPool::schedule( boost::function<void()> work, int numOmpThreads /* = 0 */)
{
    LAMA_REGION( "LAMAThreadPool::schedule" );
    boost::thread thisThread;
    bool isRecursiveTask = false;
    std::vector<boost::thread*>::const_iterator end = mThreads.end();
    for( std::vector<boost::thread*>::const_iterator it = mThreads.begin(); it != end; ++it )
    {
        if( *( *it ) == thisThread )
        {
            isRecursiveTask = true;
            break;
        }
    }

    boost::shared_ptr<LAMAThreadTask> task = LAMAThreadTask::create( work, mTaskId++, numOmpThreads );

    if( isRecursiveTask )
    {
        LAMA_LOG_WARN( logger, "Executing recursive taks synchronously to avoid deadlocks." );
        work();
        task->mState = LAMAThreadTask::FINISHED;
        return task;
    }
    boost::mutex::scoped_lock lock( mTaskQueueMutex );

    mTaskQueue.push( task );
    LAMA_LOG_DEBUG( logger, "Added task " << task->mTaskId << " to task queue" );
    task->mState = LAMAThreadTask::QUEUED;

    //  notifiy one waiting worker

    mNotifyTask.notify_one();

    return task;
}

/* ------------------------------------------------------------------------- */

void LAMAThreadPool::wait( boost::shared_ptr<LAMAThreadTask> task )
{
    LAMA_ASSERT_ERROR( task, "NULL pointer for task" );

    LAMA_LOG_DEBUG( logger, "wait on task id = " << task->mTaskId << ", state = " << task->mState );

    while( task->mState != LAMAThreadTask::FINISHED )
    {
        boost::mutex::scoped_lock lock( mNotifyFinishMutex );

        if( task->mState != LAMAThreadTask::FINISHED )
        {
            // wait on signal for a finishing thread
            // Attention: do not output here, as worker thread might finish and notify before wait

            mNotifyFinished.wait( lock );
        }

        LAMA_LOG_DEBUG( logger, "notified, task state = " << task->mState );
    }
}

/* ------------------------------------------------------------------------- */

void LAMAThreadPool::worker( int id )
{
    // This method will be executed by all threads in the pool
    // Each thread picks up a task if available in the task queue
    // No busy wait, waits on signal mNotifyTask

    LAMA_LOG_INFO( logger, "worker thread " << id << " starts" );

    mWorkerState = WORKING;

    int ompThreads = -1;
    // wait for a new task

    while( true )
    {
        boost::shared_ptr<LAMAThreadTask> task;

        // pick up a new task if available

        {
            boost::mutex::scoped_lock lock( mTaskQueueMutex );

            if( mTaskQueue.empty() )
            {
                // Instead of busy wait this thread waits on notification

                LAMA_LOG_DEBUG( logger, "worker thread " << id << " waits on notify for new task" );

                mWorkerState = WAITING;

                mNotifyTask.wait( lock );

                mWorkerState = WORKING;

                LAMA_LOG_DEBUG( logger, "worker thread " << id << " notified about a new task" );
            }
            else
            {
                task = mTaskQueue.front();
                mTaskQueue.pop();

                if( !task )
                {
                    LAMA_LOG_DEBUG( logger, "worker thread " << id << " picked shutdown task" );
                    // this is the shutdown task
                    break;
                }
            }
        }

        if( task )
        {
            LAMA_LOG_DEBUG( logger,
                            "worker thread " << id << " runs task " << task->mTaskId << " with " << task->ompThreads << " OMP threads" );

            task->mState = LAMAThreadTask::RUNNING;
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
                LAMA_LOG_INFO( logger, "worker thread got exception, has been caught" );
                task->mException = true;
            }

            boost::mutex::scoped_lock lock( mNotifyFinishMutex );

            task->mState = LAMAThreadTask::FINISHED;

            LAMA_LOG_DEBUG( logger, "worker thread " << id << " finished task " << task->mTaskId );

            // notify threads waiting on a finished task

            mNotifyFinished.notify_all();
        }
    }

    // worker is finished

    LAMA_LOG_INFO( logger, "worker thread " << id << " finishes" );
}

/* ------------------------------------------------------------------------- */

LAMAThreadPool::~LAMAThreadPool()
{
    LAMA_LOG_INFO( logger, "~LAMAThreadPool: shut down " << mThreads.size() << " threads" );

    boost::shared_ptr<LAMAThreadTask> shutdownTask; // NULL pointer

    {   // lock access to the task queue before adding shutdown tasks

        boost::mutex::scoped_lock lock( mTaskQueueMutex );

        for( size_t i = 0; i < mThreads.size(); i++ )
        {
            mTaskQueue.push( shutdownTask );
        }

        // notifiy all waiting worker threads about new task

        mNotifyTask.notify_all();
    }

    LAMA_LOG_DEBUG( logger, "added " << mThreads.size() << " shutdown tasks" );

    // and now wait for completion of all worker threads and delete them

    for( size_t i = 0; i < mThreads.size(); ++i )
    {
        LAMA_LOG_DEBUG( logger, "wait for worker thread " << i );
        mThreads[i]->join();
        LAMA_LOG_DEBUG( logger, "worker thread " << i << " terminated (joined)" );
        delete mThreads[i];
        mThreads[i] = NULL;
    }
}

} // namespace lama
