/**
 * @file Task.hpp
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
 * @brief Declaration of the class Task for asynchronous function execution.
 * @author Thomas Brandes
 * @date 21.11.2011
 * $Id$
 */
#ifndef LAMA_TASK_HPP_
#define LAMA_TASK_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Printable.hpp>

// others
#include <lama/task/LAMAThreadPool.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

namespace lama
{

/** @brief Class for asynchronous tasks.
 *
 * Task represents a asynchronous task, which allows synchronization. To start an
 * asynchronous task one simply have to construct a Task object and pass a pointer to
 * the function to execute to its constructor. If the signature of the function does not fit
 * boost::bind can be used to build an appropriate function pointer. Be careful with reference because they
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
 *     boost::function<void()> function = boost::bind( foo, boost::ref(y), alpha ,boost::cref(A), boost::cref(x) );
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
class LAMA_DLL_IMPORTEXPORT Task: public Printable
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
     *  By using boost::function, it is possible to bind arguments via boost::bind.
     */
    Task( boost::function<void()> function, int numOmpThreads = 0 );

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

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    /** Thread pool will be allocated on need, shared pointer guarantees deallocation. */

    static boost::shared_ptr<LAMAThreadPool> theThreadPool;

    static LAMAThreadPool& getThreadPool();

    boost::shared_ptr<LAMAThreadTask> mTask;
};

} // namespace lama

#endif // LAMA_TASK_HPP_
