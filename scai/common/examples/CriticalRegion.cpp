/**
 * @file common/examples/CriticalRegion.cpp
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
 * @brief Example with thhreads and a critical region using the common library
 * @author Thomas Brandes
 * @date 19.06.2015
 */

#include <scai/common/Walltime.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <thread>
#include <mutex>

using namespace std;
using namespace scai::common;

std::recursive_mutex threadRecursiveMutex;    // recursive threadRecursiveMutex needed here

static const int SLEEP_TIME  = 2;
static const int N_THREADS   = 4;

// Define routine that is executed by one thread

static void threadRoutine()
{
    std::thread::id self = std::this_thread::get_id();
    cout << "Thread " << self << " starts" << endl;
    std::unique_lock<std::recursive_mutex> lock( threadRecursiveMutex );
    std::unique_lock<std::recursive_mutex> lock1( threadRecursiveMutex );   
    // second lock by same thread is okay for recursive threadRecursiveMutex
    cout << "Thread " << self << " enters critical region" << endl;
    Walltime::sleep( SLEEP_TIME * 1000 );
    cout << "Thread " << self << " leaves critical region" << endl;
}

int main( int, char** )
{
    // macro to give the current thread a name that appears in further logs

    double time = Walltime::get();

    std::thread threads[N_THREADS];

    for ( int i = 0; i < N_THREADS; ++i )
    {
        threads[i] = std::thread( threadRoutine );
    }

    for ( int i = 0; i < N_THREADS; ++i )
    {
        threads[i].join();
    }

    time = Walltime::get() - time;
    cout << "Termination of " << N_THREADS << " threads after " << time << " seconds" << endl;
    // If critical region is implemented correctly, time must be > ( #threds * sleep_time )
    SCAI_ASSERT_LT( N_THREADS * SLEEP_TIME, time,
                    "ERROR: " << N_THREADS << " threads seem to enter critial region at same time" )
}
