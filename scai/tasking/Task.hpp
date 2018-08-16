/**
 * @file Task.hpp
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
 * @brief Declaration of the class Task for asynchronous function execution.
 * @author Thomas Brandes
 * @date 21.11.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/tasking/ThreadPool.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <memory>
#include <functional>

namespace scai
{

namespace tasking
{

/** @brief Class for asynchronous tasks.
 *
 * Task represents a asynchronous task, which allows synchronization. To start an
 * asynchronous task one simply have to construct a Task object and pass a pointer to
 * the function to execute to its constructor. If the signature of the function does not fit
 * std::bind can be used to build an appropriate function pointer. Be careful with reference because they
 * need special wrapper objects. E.g.
 * \code
 *
 * void foo( Vector& y, const Matrix& A, Scalar alpha, const Vector& x);
 *
 * ...
 *
 * void bar( Vector& y, const Matrix& A, const Scalar alpha, const Vector& x)
 * {
 *     //construct function pointer
 *     std::function<void()> function = std::bind( foo, std::ref(y), alpha, std::cref(A), std::cref(x) );
 *
 *     //start asynchronous tasks
 *     Task task( function );
 *
 *     //do other work
 *
 *     //wait for asynchronous task to finish
 *     task.synchronize();
 * }
 * \endcode
 */
class COMMON_DLL_IMPORTEXPORT Task: public common::Printable
{

public:

    typedef unsigned int Id;

    /** Construct a new task that will execute a void function, no arguments.
     *
     * Its safe to start a recursive Task, e.g. to start a Task from a Task. The
     * recursive Task might execute synchronously to avoid deadlock situations.
     *
     *  @param[in] function      is void() function executed by the new created task.
     *  @param[in] numOmpThreads number of openmp threads the task should use
     *                           (if numOmpThreads == 0 omp_get_max_threads() is used)
     *
     *  By using std::function, it is possible to bind arguments via std::bind.
     */
    Task( std::function<void()> function, int numOmpThreads = 0 );

    /**
     *  Destructor. Implicitly synchronizes the Task.
     */
    ~Task();

    /**
     *  Method synchronize waits for finalization of the task.
     */
    void synchronize();

    /**
     * Method test if this Task is already finished.
     *
     * @return if this Task is already finished.
     */
    bool probe() const;

    /**
     *  Query a task for its id.
     */
    Id getId() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Thread pool will be allocated on need, shared pointer guarantees deallocation. */

    static std::shared_ptr<ThreadPool> theThreadPool;

    static ThreadPool& getThreadPool();

    std::shared_ptr<ThreadPoolTask> mTask;
};

} /* end namespace tasking */

} /* end namespace scai */
