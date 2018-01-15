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

static void runRoutine( bool& okay, const int i )
{
    std::string name = "Thread_" + std::to_string( i );
    thread::defineCurrentThreadName( name.c_str() );
    barrier();  // Barrier 1 : all threads have defined their name
    std::string name1 = *getCurrentThreadName();
    // BOOST_CHECK_EQUAL is not thread-safe
    okay = name == name1;
    barrier();  // Barrier 2 : checks are done
    name = "NewName_" + std::to_string( i );
    thread::defineCurrentThreadName( name.c_str() );
    barrier();  // Barrier 3 : names are redefined
}

BOOST_AUTO_TEST_CASE( nameThreadTest )
{
    // run a number of threads and give them names

    std::thread threads[N_THREADS];
    bool results[N_THREADS];
    std::shared_ptr<std::string> names[N_THREADS];

    for ( int i = 0; i < N_THREADS; ++i )
    {
        results[i] = false; 
        threads[i] = std::thread( runRoutine, std::ref( results[i] ), i );
    }

    barrier();  

    for ( int i = 0; i < N_THREADS; ++i )
    {
        std::string expected = "Thread_" + std::to_string( i );
        names[i] = getThreadName( threads[i].get_id() );
        BOOST_CHECK_EQUAL( expected, *names[i] );
    }

    barrier(); // checks done
    barrier(); // thread names are now redefined

    for ( int i = 0; i < N_THREADS; ++i )
    {
        std::string expected = "Thread_" + std::to_string( i );
        BOOST_CHECK_EQUAL( expected, *names[i] );

        expected = "NewName_" + std::to_string( i );
        std::string name = *getThreadName( threads[i].get_id() );
        BOOST_CHECK_EQUAL( expected, name );
    }

    for ( int i = 0; i < N_THREADS; ++i )
    {
        threads[i].join();
        BOOST_CHECK( results[i] );   // get correct results
    }
}

BOOST_AUTO_TEST_SUITE_END()
