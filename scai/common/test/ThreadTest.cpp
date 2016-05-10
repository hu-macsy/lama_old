/**
 * @file common/test/ThreadTest.cpp
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
 * @brief Test program for threads
 *
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Thread.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <unistd.h>

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

    int expected_sum = NB_THREADS * ( NB_THREADS-1 ) / 2;

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

static const int SLEEP_TIME  = 1;
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

    sleep( SLEEP_TIME );

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

    BOOST_CHECK( C_THREADS * SLEEP_TIME <= time );
}

// Define routine that is executed by one thread

static void runRoutine( int& )
{
    sleep( SLEEP_TIME );
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

