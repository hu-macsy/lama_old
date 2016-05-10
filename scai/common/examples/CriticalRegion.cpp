/**
 * @file common/examples/CriticalRegion.cpp
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
 * @brief Example with pthreads and a critical region using the common library 
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/common/Thread.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <unistd.h>

using namespace std;
using namespace scai::common;

Thread::RecursiveMutex threadRecursiveMutex;    // recursive threadRecursiveMutex needed here

static int SLEEP_TIME  = 2;
static int N_THREADS   = 4;

// Define routine that is executed by one thread

void* threadRoutine( void* )
{
    Thread::Id self = Thread::getSelf();

    cout << "Thread " << self << " starts" << endl;

    Thread::ScopedLock lock( threadRecursiveMutex );
    Thread::ScopedLock lock1( threadRecursiveMutex );   // second lock by same thread is okay for recursive threadRecursiveMutex

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
