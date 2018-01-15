/**
 * @file common/test/ThreadTest.cpp
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
 * @brief Test program for naming of threads
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/thread.hpp>
#include <scai/common/macros/assert.hpp>

#include <iostream>
#include <cstdlib>

using namespace scai;
using namespace common;
using namespace thread;

BOOST_AUTO_TEST_SUITE( ThreadTest )

static const int N_THREADS   = 31;  // number of threads

std::mutex barrierMutex;
std::condition_variable_any barrierCondition;

static int thread_cnt = 0;

static void barrier()
{
    std::unique_lock<std::mutex> lock( barrierMutex );
    thread_cnt ++;

    if ( thread_cnt != N_THREADS + 1 )
    {
        // Some others not at barrier so wait
        barrierCondition.wait( lock );
    }
    else
    {
        // Now all threads have reached
        thread_cnt = 0;
        barrierCondition.notify_all();
    }
}

/** Function that is that is executed by each thread */

static void runRoutine( const int i )
{
    std::string name = "Thread_" + std::to_string( i );
    thread::defineCurrentThreadName( name.c_str() );
    barrier();
    std::string name1 = getCurrentThreadName();
    // BOOST_CHECK_EQUAL on separate threads can produce strange outputs
    SCAI_ASSERT_EQ_ERROR( name, name1, "serious error" );
}

BOOST_AUTO_TEST_CASE( nameThreadTest )
{
    // run a number of threads and give them names

    std::thread threads[N_THREADS];

    for ( int i = 0; i < N_THREADS; ++i )
    {
        threads[i] = std::thread( runRoutine, i );
    }

    barrier();

    for ( int i = 0; i < N_THREADS; ++i )
    {
        std::string expected = "Thread_" + std::to_string( i );
        std::string name = getThreadName( threads[i].get_id() );
        BOOST_CHECK_EQUAL( expected, name );
    }

    for ( int i = 0; i < N_THREADS; ++i )
    {
        threads[i].join();
    }
}

BOOST_AUTO_TEST_SUITE_END()
