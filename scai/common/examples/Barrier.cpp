/**
 * @file common/examples/Barrier.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Example with pthreads using Condition of the common library to implement barrier.
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/common/Thread.hpp>
#include <scai/common/macros/assert.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <unistd.h>

using scai::common::Thread;

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

    SCAI_ASSERT_EQUAL( sum, expected_sum, "Wrong value after thread barrier" )

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
