/**
 * @file ThreadPool.hpp
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
 * @brief Definition of class for a thread pool
 * @author Thomas Brandes
 * @date 27.12.2012
 */

#pragma once

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/config.hpp>
#include <scai/common/NonCopyable.hpp>

// std
#include <climits>
#include <queue>
#include <memory>
#include <functional>

namespace scai
{

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

struct COMMON_DLL_IMPORTEXPORT ThreadPoolTask
{
    enum TaskState
    {
        DEFINED,    //!< The task structure has been filled
        QUEUED,     //!< Task is queued in the pool, ready for execution
        RUNNING,    //!< Task is executed by free thread
        FINISHED    //!< Task is terminated, structure still exists
    };

    std::function<void()> mWork;  //!< task function to be executed

    volatile TaskState mState; //!< current state of the task

    bool mException; //!< true if task got an exception

    unsigned int mTaskId; //!< unique task id

    int ompThreads; //!< number of OpenMP threads to be used

    /** Create a new task as a shared pointer */

    static std::shared_ptr<ThreadPoolTask> create(
        std::function<void()> work,
        unsigned int taskId,
        int numOmpThreads = 0 );
};

/** Class that defines a pool of threads that will execute scheduled tasks
 (see ThreadPoolTask).

 The main advantage of this class is that creation of running threads is
 only done once. Non-busy threads will wait on notification for new tasks.

 Note: tasks will be handled by shared pointers so task data might be freed
 either here or at its creation.
 */

class COMMON_DLL_IMPORTEXPORT ThreadPool : private common::NonCopyable
{

public:

    /** Create a new thread pool with a certain number of threads. */

    ThreadPool( int size );

    int size() const
    {
        return mMaxSize;
    }

    /** Schedules a new function that will be executed by available free thread.
     *
     *  @param[in] work      is function to be executed by available thread
     *  @param[in] numOmpThreads number of openmp threads the task should use
     *               (if numOmpThreads == 0 omp_get_max_threads() is used)
     *  @return shared pointer for the task
     */

    std::shared_ptr<ThreadPoolTask> schedule( std::function<void()> work, int numOmpThreads = 0 );

    /** Wait on completion of a task. */

    void wait( std::shared_ptr<ThreadPoolTask> task );

    /** Wait on completion of all scheduled tasks */

    void shutdown();

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

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    unsigned int mTaskId; // last given task id

    int mMaxSize;         // number of worker threads

    // use scoped array instead of vector as no copy constructor is available

    std::unique_ptr<std::thread[]> mThreads;    // worker threads of this pool
    std::unique_ptr<ThreadData[]> mThreadArgs;     // arguments for each worker thread

    std::queue<std::shared_ptr<ThreadPoolTask> > mTaskQueue;

    std::condition_variable_any mNotifyFinished;// notify about finished tasks
    std::condition_variable_any mNotifyTask;// notify about new task

    std::mutex mTaskQueueMutex;      // make access to taskqueue thread-safe
    std::mutex mNotifyFinishMutex;   // used for wait on mNotifyFinished

    enum WorkerState
    {
        WORKING,
        WAITING
    };

    WorkerState mWorkerState; //used to indicate the state of the worker
};

} /* end namespace tasking */

} /* end namespace scai */
