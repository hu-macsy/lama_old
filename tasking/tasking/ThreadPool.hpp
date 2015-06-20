/**
 * @file ThreadPool.hpp
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
 * @brief Definition of class for a thread pool
 * @author Thomas Brandes
 * @date 27.12.2012
 * @since 1.0.0
 */

#pragma once

#include <common/config.hpp>
#include <common/Thread.hpp>

#include <logging/logging.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/version.hpp>

#include <pthread.h>

#include <climits>
#include <queue>

namespace tasking
{

/** Task identifiers are given by an int value. */

typedef unsigned int TaskId;

/** Structure that specifies worker task and its state.
 *
 *  As this structure might be used either by a running thread or by
 *  the thread that issued the task it should always be used via a
 *  shared pointer so that ownership never matters.
 */

struct COMMON_DLL_IMPORTEXPORT ThreadTask
{
    enum TaskState
    {
        DEFINED,    //!< The task structure has been filled
        QUEUED,     //!< Task is queued in the pool, ready for execution
        RUNNING,    //!< Task is executed by free thread
        FINISHED    //!< Task is terminated, structure still exists
    };

    boost::function<void()> mWork;  //!< task function to be executed

    volatile TaskState mState; //!< current state of the task

    bool mException; //!< true if task got an exception

    unsigned int mTaskId; //!< unique task id

    int ompThreads; //!< number of OpenMP threads to be used

    /** Create a new task as a shared pointer */

    static boost::shared_ptr<ThreadTask> create(
        boost::function<void()> work,
        unsigned int taskId,
        int numOmpThreads = 0 );
};

/** Class that defines a pool of threads that will execute scheduled tasks
 (see ThreadTask).

 The main advantage of this class is that creation of running threads is
 only done once. Non-busy threads will wait on notification for new tasks.

 Note: tasks will be handled by shared pointers so task data might be freed
 either here or at its creation.
 */

class COMMON_DLL_IMPORTEXPORT ThreadPool
{

public:

    /** Create a new thread pool with a certain number of threads. */

    ThreadPool( int size );

    /** Schedules a new function that will be executed by available free thread.
     *
     *  @param[in] work      is function to be executed by available thread
     *  @param[in] numOmpThreads number of openmp threads the task should use
     *               (if numOmpThreads == 0 omp_get_max_threads() is used)
     *  @return shared pointer for the task
     */

    boost::shared_ptr<ThreadTask> schedule( boost::function<void()> work, int numOmpThreads = 0 );

    /** Wait on completion of a task. */

    void wait( boost::shared_ptr<ThreadTask> task );

    /** Destructor shuts down and deletes all threads of the pool. */

    ~ThreadPool();

    /** Arguments that will be passed to each worker thread */

    struct ThreadData
    {
        ThreadPool* pool;
        int         i;
    };

    /** Method that will be run by the worker threads. */

    void worker( int id );

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    unsigned    int mTaskId; // last given task id

    int mMaxSize;// number of worker threads

    std::vector<pthread_t> mThreads;     // worker threads of this pool
    std::vector<ThreadData> mThreadArgs; // arguments for each worker thread

    std::queue<boost::shared_ptr<ThreadTask> > mTaskQueue;

    common::Thread::Condition mNotifyFinished;// notify about finished tasks
    common::Thread::Condition mNotifyTask;// notify about new task

    common::Thread::Mutex mTaskQueueMutex;// make access to taskqueue thread-safe
    common::Thread::Mutex mNotifyFinishMutex;// used for wait on mNotifyFinished

    enum WorkerState
    {   
        WORKING, 
        WAITING
    };

    WorkerState mWorkerState; //used to indicate the state of the worker
};

} // namespace 

