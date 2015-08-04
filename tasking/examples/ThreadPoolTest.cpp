/**
 * @file ThreadPoolTest.hpp
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
 * @brief Extensive test program for ThreadPool
 * @author Thomas Brandes
 * @date 13.07.2015
 */

#include <tracing.hpp>
#include <tasking/ThreadPool.hpp>

#include <memory>
#include <common/bind.hpp>

using namespace common;
using namespace tasking;

LAMA_LOG_DEF_LOGGER( logger, "Test.ThreadPoolTest" )

/** Maximal number of threads in the pool. */

#define POOL_SIZES { 1, 2, 3 }

/** Number of tasks that will be scheduled. */

#define TASK_SIZES { 1, 8, 15 }

/** Factor used for workload of one task. */

static const int WORKLOAD = 1000000;

/* ----------------------------------------------------------------------- */

void work( const int in, int& out )
{
    out = in;
    int factor = in % 4 + 1;
    std::ostringstream regionName;
    regionName << "work_" << factor;
    LAMA_REGION( regionName.str().c_str() )

    // just do some stupid work, workload depends on in

    for ( int i = 0; i < WORKLOAD * factor; i++ )
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

void runTest()
{
    LAMA_LOG_INFO( logger, "runTest" );
    LAMA_REGION( "runTest" )
    int thread_sizes[] = POOL_SIZES;
    int thread_configs = sizeof( thread_sizes ) / sizeof( int );
    int task_sizes[] = TASK_SIZES;
    int task_configs = sizeof( task_sizes ) / sizeof( int );

    for ( int i = 0; i < task_configs; ++i )
    {
        int ntasks = task_sizes[i];

        for ( int j = 0; j < thread_configs; ++j )
        {
            LAMA_REGION_N( "PoolRun", thread_sizes[j] * 100 + ntasks )

            std::vector<int> x ( ntasks );

            {
                ThreadPool pool( thread_sizes[j] );

                // Just issue the tasks, do not keep references

                for ( int i = 0; i < ntasks; i++ )
                {
                    x[i] = -1;
                    pool.schedule( bind( &work, i, ref( x[i] ) ) );
                }

                // end of scope: pool waits for all tasks/work to be finished
            }

            for ( int i = 0; i < ntasks; i++ )
            {
                // BOOST_CHECK_EQUAL( i, x[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

void waitTest()
{
    LAMA_LOG_INFO( logger, "waitTest" );
    LAMA_REGION( "waitTest" )
    int thread_sizes[] = POOL_SIZES;
    int thread_configs = sizeof( thread_sizes ) / sizeof( int );
    int task_sizes[] = TASK_SIZES;
    int task_configs = sizeof( task_sizes ) / sizeof( int );

    for ( int i = 0; i < task_configs; ++i )
    {
        int ntasks = task_sizes[i];

        for ( int j = 0; j < thread_configs; ++j )
        {
            LAMA_REGION_N( "PoolWait", thread_sizes[j] * 100 + ntasks )
            std::vector<int> x( ntasks );  // array with result for each task
            std::vector<shared_ptr<ThreadTask> > tasks( ntasks );
            ThreadPool pool( thread_sizes[j] );

            for ( int i = 0; i < ntasks; i++ )
            {
                x[i] = -1;
                tasks[i] = pool.schedule( bind( &work, i, ref( x[i] ) ) );
            }

            for ( int i = 0; i < ntasks; i++ )
            {
                pool.wait( tasks[i] );
                // BOOST_CHECK_EQUAL( i, x[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

void singleTest()
{
    LAMA_LOG_THREAD( "main:singleTest" )

    LAMA_REGION( "singleTest" )
    LAMA_LOG_INFO( logger, "singleTest" );
    // Extensive test of a thread pool with one thread
    // Should verify that master never misses a notify at wait
    ThreadPool pool( 1 );
    int rnd = 15;

    const int NTIMES = 100;

    for ( int i = 0; i < NTIMES; ++i )
    {
        int resultThread;
        int resultMaster;

        shared_ptr<ThreadTask> task = pool.schedule( bind( &work, i, ref( resultThread ) ) );

        // Master thread does something and then waits
        rnd = ( rnd + 19 ) % 17;
        work( i + rnd, ref( resultMaster ) );
        // BOOST_CHECK_EQUAL( i + rnd, resultMaster );
        pool.wait( task );
        // BOOST_CHECK_EQUAL( i, resultThread );
    }
}

/* --------------------------------------------------------------------- */

int main()
{
    singleTest();
    runTest();
    waitTest();
}
