/**
 * @file common/examples/Barrier.cpp
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
 * @brief Example with pthreads using Condition of the common library to implement barrier.
 *
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/common/Thread.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/Exception.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <unistd.h>

using namespace std;
using namespace scai::common;

Thread::Mutex barrierMutex;
Thread::Condition barrierCondition;

Thread::Mutex printMutex;

static const int N_THREADS   = 16;

// Define routine that is executed by one thread

static int thread_cnt = 0;

static void barrier()
{
    Thread::ScopedLock lock( barrierMutex );

    thread_cnt ++;
     
    if ( thread_cnt != N_THREADS )
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

static int sharedArray[ N_THREADS ];

static void threadRoutine( int& arg )
{
    sharedArray[arg] = arg;

    barrier();

    int sum = 0;

    for ( int i = 0; i < N_THREADS; ++i )
    {
        sum += sharedArray[i];
    }

    int expected_sum = N_THREADS * ( N_THREADS-1 ) / 2;

    COMMON_ASSERT_EQUAL( sum, expected_sum, "Wrong value after thread barrier" )

    {
        Thread::ScopedLock lock( printMutex );
        std::cout << "Thread " << arg << " has correct sum = " << sum << std::endl;
    }
}

int main( int, char** )
{
    Thread threads[N_THREADS];
    int threadArgs[N_THREADS];

    for ( int i = 0; i < N_THREADS; ++i )
    {
        threadArgs[i] = i;
        threads[i].run( threadRoutine, threadArgs[i] );
    }

    for ( int i = 0; i < N_THREADS; ++i )
    {
        threads[i].join();
    }

    std::cout << "All threads are terminated correctly." << std::endl;
}
