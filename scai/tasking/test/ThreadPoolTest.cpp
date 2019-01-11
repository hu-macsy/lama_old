/**
 * @file ThreadPoolTest.cpp
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
 * @brief Extensive test program for ThreadPool
 * @author Thomas Brandes
 * @date 13.07.2015
 */

#include <boost/test/unit_test.hpp>

#include <scai/tracing.hpp>
#include <scai/tasking/ThreadPool.hpp>

#include <memory>
#include <functional>

using namespace scai::common;
using namespace scai::tasking;

using std::shared_ptr;
using std::bind;
using std::ref;

/** Maximal number of threads in the pool. */

#define POOL_SIZES { 1, 2, 3 }

/** Number of tasks that will be scheduled. */

#define TASK_SIZES { 1, 8, 15 }

/** Factor used for workload of one task. */

static const int WORKLOAD = 1000000;

/* ----------------------------------------------------------------------- */

static void work( const int in, int& out )
{
    out = in;
    int factor = in % 4 + 1;

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

BOOST_AUTO_TEST_SUITE( ThreadPoolTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.ThreadPoolTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( runTest )
{
    SCAI_LOG_INFO( logger, "runTest" );
    int thread_sizes[] = POOL_SIZES;
    int thread_configs = sizeof( thread_sizes ) / sizeof( int );
    int task_sizes[] = TASK_SIZES;
    int task_configs = sizeof( task_sizes ) / sizeof( int );

    for ( int i = 0; i < task_configs; ++i )
    {
        int ntasks = task_sizes[i];

        for ( int j = 0; j < thread_configs; ++j )
        {
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
                BOOST_CHECK_EQUAL( i, x[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( waitTest )
{
    SCAI_LOG_INFO( logger, "waitTest" );
    int thread_sizes[] = POOL_SIZES;
    int thread_configs = sizeof( thread_sizes ) / sizeof( int );
    int task_sizes[] = TASK_SIZES;
    int task_configs = sizeof( task_sizes ) / sizeof( int );

    for ( int i = 0; i < task_configs; ++i )
    {
        int ntasks = task_sizes[i];

        for ( int j = 0; j < thread_configs; ++j )
        {
            std::vector<int> x( ntasks );  // array with result for each task
            std::vector<shared_ptr<ThreadPoolTask> > tasks( ntasks );
            ThreadPool pool( thread_sizes[j] );

            for ( int i = 0; i < ntasks; i++ )
            {
                x[i] = -1;
                tasks[i] = pool.schedule( bind( &work, i, ref( x[i] ) ) );
            }

            for ( int i = 0; i < ntasks; i++ )
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
    SCAI_LOG_THREAD( "main:singleTest" )
    SCAI_LOG_INFO( logger, "singleTest" );
    // Extensive test of a thread pool with one thread
    // Should verify that master never misses a notify at wait
    ThreadPool pool( 1 );
    int rnd = 15;
    const int NTIMES = 100;

    for ( int i = 0; i < NTIMES; ++i )
    {
        int resultThread;
        int resultMaster;
        shared_ptr<ThreadPoolTask> task = pool.schedule( bind( &work, i, ref( resultThread ) ) );
        // Master thread does something and then waits
        rnd = ( rnd + 19 ) % 17;
        work( i + rnd, ref( resultMaster ) );
        BOOST_CHECK_EQUAL( i + rnd, resultMaster );
        pool.wait( task );
        BOOST_CHECK_EQUAL( i, resultThread );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
