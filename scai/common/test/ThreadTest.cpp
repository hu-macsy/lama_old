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
 * @brief Test program for threads
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Thread.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>

using namespace scai;
using namespace common;

Thread::Mutex barrierMutex;
Thread::Condition barrierCondition;

static const int NB_THREADS   = 16;

// Define routine that is executed by one thread

static int thread_cnt = 0;

static void barrier()
{
    Thread::ScopedLock lock( barrierMutex );
    thread_cnt ++;

    if ( thread_cnt != NB_THREADS )
    {
        // Some others not at barrier so wait
        barrierCondition.wait( lock );
    }
    else
    {
        // Now all threads have reached
        thread_cnt = 0;
        barrierCondition.notifyAll();
    }
}

// Shared array for all threads

static int sharedArray[ NB_THREADS ];

static void barrierRoutine( int& arg )
{
    sharedArray[arg] = arg;
    barrier();
    int sum = 0;

    for ( int i = 0; i < NB_THREADS; ++i )
    {
        sum += sharedArray[i];
    }

    int expected_sum = NB_THREADS * ( NB_THREADS - 1 ) / 2;
    BOOST_CHECK_EQUAL( sum, expected_sum );
}

BOOST_AUTO_TEST_CASE( barrierTest )
{
    Thread threads[NB_THREADS];
    int threadArgs[NB_THREADS];

    for ( int i = 0; i < NB_THREADS; ++i )
    {
        threadArgs[i] = i;
        threads[i].run( barrierRoutine, threadArgs[i] );
    }

    for ( int i = 0; i < NB_THREADS; ++i )
    {
        threads[i].join();
    }
}

static const int SLEEP_TIME  = 1;  // in seconds
static const int C_THREADS   = 4;

// Define routine that is executed by one thread

Thread::RecursiveMutex critMutex;    // recursive mutex needed here

static void criticalRoutine( int& n )
{
    std::ostringstream nstream;
    nstream << "Thread_" << n;
    Thread::defineCurrentThreadName( nstream.str().c_str() );
    Thread::ScopedLock lock( critMutex );
    Thread::ScopedLock lock1( critMutex );   // second lock by same thread is okay for recursive mutex
    Walltime::sleep( SLEEP_TIME * 1000 );
    BOOST_CHECK_EQUAL( nstream.str(), Thread::getCurrentThreadName() );
}

BOOST_AUTO_TEST_CASE( criticalRegionTest )
{
    // macro to give the current thread a name that appears in further logs
    Thread threads[C_THREADS];
    int threadArgs[C_THREADS];
    double time = Walltime::get();

    for ( int i = 0; i < C_THREADS; ++i )
    {
        threadArgs[i] = i;
        threads[i].run( criticalRoutine, threadArgs[i] );
    }

    for ( int i = 0; i < C_THREADS; ++i )
    {
        threads[i].join();
    }

    time = Walltime::get() - time;
    // If critical region is implemented correctly, time must be > ( #threds * sleep_time )
    BOOST_CHECK( C_THREADS* SLEEP_TIME <= time );
}

// Define routine that is executed by one thread

static void runRoutine( int& )
{
    Walltime::sleep( SLEEP_TIME * 1000 );
}

BOOST_AUTO_TEST_CASE( concurrentTest )
{
    // macro to give the current thread a name that appears in further logs
    Thread threads[C_THREADS];
    int threadArgs[C_THREADS];
    double time = Walltime::get();

    for ( int i = 0; i < C_THREADS; ++i )
    {
        threadArgs[i] = i;
        threads[i].run( runRoutine, threadArgs[i] );
    }

    for ( int i = 0; i < C_THREADS; ++i )
    {
        threads[i].join();
    }

    time = Walltime::get() - time;
    // If threads are implemented correctly, time must be ~ sleep_time
    BOOST_CHECK_CLOSE( time, double( SLEEP_TIME ), 20 );
}

