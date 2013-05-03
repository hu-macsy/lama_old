/**
 *  * @file ThreadPoolTest.hpp
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
 * @brief Contains the implementation of the class ThreadPoolTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 02.02.2012
 * $Id$
 */

#include <boost/test/unit_test.hpp>

#include <lama/CommunicatorFactory.hpp>
#include <lama/task/LAMAThreadPool.hpp>
#include <lama/tracing.hpp>

#include <memory>
#include <boost/bind.hpp>
#include <boost/scoped_array.hpp>

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ThreadPoolTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.ThreadPoolTest" );

lama::CommunicatorPtr mComm; //!< communicator used for distribution, parallel execution

/** Maximal number of threads in the pool. */

#define LAMA_THREADS { 1, 2, 3 }

/** Number of tasks that will be scheduled. */

#define LAMA_TASKS { 1, 8, 15 }

/** Factor used for workload of one task. */

const int LAMA_WORKLOAD = 1000000;
/* ----------------------------------------------------------------------- */

void work( const int in, int& out )
{
    out = in;

    int factor = in % 4 + 1;

    std::ostringstream regionName;

    regionName << "work_" << factor;
    LAMA_REGION( regionName.str().c_str() )

    // just do some stupid work, workload depends on in

    for ( int i = 0; i < LAMA_WORKLOAD * factor; i++ )
    {
        int dir = i & 1;
        if ( dir )
        {
            out += 13;
        }
        else
        {
            out -= 13;
        }
    }
    out = in;
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( runTest )
{
    LAMA_LOG_INFO( logger, "runTest" );

    LAMA_REGION( "runTest" )

    IndexType thread_sizes[] = LAMA_THREADS;
    IndexType thread_configs = sizeof( thread_sizes ) / sizeof(IndexType);

    IndexType task_sizes[] = LAMA_TASKS;
    IndexType task_configs = sizeof( task_sizes ) / sizeof(IndexType);

    for ( IndexType i = 0; i < task_configs; ++i )
    {
        IndexType ntasks = task_sizes[i];

        for ( IndexType j = 0; j < thread_configs; ++j )
        {
            LAMA_REGION_N( "PoolRun", thread_sizes[j] * 100 + ntasks )

            scoped_array<IndexType> x( new IndexType[ntasks] );

            {
                LAMAThreadPool pool( thread_sizes[j] );

                // Just issue the tasks, do not keep references

                for ( IndexType i = 0; i < ntasks; i++ )
                {
                    x[i] = -1;
                    pool.schedule( boost::bind( &ThreadPoolTest::work, i, ref( x[i] ) ) );
                }

                // end of scope: pool waits for all tasks/work to be finished
            }

            for ( IndexType i = 0; i < ntasks; i++ )
            {
                BOOST_CHECK_EQUAL( i, x[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( waitTest )
{
    LAMA_LOG_INFO( logger, "waitTest" );

    LAMA_REGION( "waitTest" )

    IndexType thread_sizes[] = LAMA_THREADS;
    IndexType thread_configs = sizeof( thread_sizes ) / sizeof(IndexType);

    IndexType task_sizes[] = LAMA_TASKS;
    IndexType task_configs = sizeof( task_sizes ) / sizeof(IndexType);

    for ( IndexType i = 0; i < task_configs; ++i )
    {
        int ntasks = task_sizes[i];

        for ( IndexType j = 0; j < thread_configs; ++j )
        {
            LAMA_REGION_N( "PoolWait", thread_sizes[j] * 100 + ntasks )

            scoped_array<IndexType> x( new IndexType[ntasks] ); // array with result for each task

            scoped_array<shared_ptr<LAMAThreadTask> > tasks( new shared_ptr<LAMAThreadTask> [ntasks] );

            LAMAThreadPool pool( thread_sizes[j] );

            for ( IndexType i = 0; i < ntasks; i++ )
            {
                x[i] = -1;
                tasks[i] = pool.schedule( boost::bind( &ThreadPoolTest::work, i, ref( x[i] ) ) );
            }

            for ( IndexType i = 0; i < ntasks; i++ )
            {
                pool.wait( tasks[i] );
                BOOST_CHECK_EQUAL( i, x[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( singleTest )
{
    LAMA_REGION( "singleTest" )

    LAMA_LOG_INFO( logger, "singleTest" );

    // Extensive test of a thread pool with one thread
    // Should verify that master never misses a notify at wait

    LAMAThreadPool pool( 1 );

    IndexType rnd = 15;

    for ( IndexType i = 0; i < 100; ++i )
    {
        IndexType resultThread;
        IndexType resultMaster;

        boost::shared_ptr<LAMAThreadTask> task = pool.schedule(
                    boost::bind( &ThreadPoolTest::work, i, ref( resultThread ) ) );

        // Master thread does something and then waits

        rnd = ( rnd + 19 ) % 17;

        work( i + rnd, ref( resultMaster ) );

        BOOST_CHECK_EQUAL( i + rnd, resultMaster );

        pool.wait( task );

        BOOST_CHECK_EQUAL( i, resultThread );
    }
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
