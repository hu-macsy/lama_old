/**
 * @file common/examples/CriticalRegion.cpp
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
 * @brief Example with pthreads and a critical region using the common library 
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

Thread::RecursiveMutex mutex;    // recursive mutex needed here

static int SLEEP_TIME  = 2;
static int N_THREADS   = 4;

// Define routine that is executed by one thread

void* threadRoutine( void* )
{
    Thread::Id self = Thread::getSelf();

    cout << "Thread " << self << " starts" << endl;

    Thread::ScopedLock lock( mutex );
    Thread::ScopedLock lock1( mutex );   // second lock by same thread is okay for recursive mutex

    cout << "Thread " << self << " enters critical region" << endl;
    sleep( SLEEP_TIME );
    cout << "Thread " << self << " leaves critical region" << endl;

    return NULL;
}

int main( int, char** )
{
    // macro to give the current thread a name that appears in further logs

    std::vector<pthread_t> threads;

    threads.resize( N_THREADS );

    double time = Walltime::get();

    for ( int i = 0; i < N_THREADS; ++i )
    {
        void* arg = NULL;

        int rc = pthread_create( &threads[i], NULL, &threadRoutine, arg );

        if ( rc != 0 )
        {
            COMMON_THROWEXCEPTION( "Could not create pthread " << i << " of " << N_THREADS << ", rc = " << rc )
        }
    }

    for ( int i = 0; i < N_THREADS; ++i )
    {
        pthread_join( threads[i], NULL );
    }

    time = Walltime::get() - time;

    cout << "Termination of " << N_THREADS << " threads after " << time << " seconds" << endl;

    // If critical region is implemented correctly, time must be > ( #threds * sleep_time )

    SCAI_ASSERT_LT( N_THREADS * SLEEP_TIME, time, 
                      "ERROR: " << N_THREADS << " threads seem to enter critial region at same time" )
}
