/**
 * @file BenchTasking.cpp
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
 * @brief Compare performance of Task vs. Thread
 * @author Thomas Brandes
 * @date 02.02.2012
 */

#include <scai/tasking/Task.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/macros/throw.hpp>

#include <memory>
#include <thread>
#include <functional>

using namespace std;
using namespace scai;
using namespace scai::tasking;

static const int WORKLOAD = 2000;

/* ----------------------------------------------------------------------- */

void work( int& out )
{
    int in = out;
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
/* ----------------------------------------------------------------- */

void doTasking( int N )
{
    std::unique_ptr<int[]> arg( new int[N] );
    std::unique_ptr<Task*[]> tasks( new Task*[N] );

    for ( int i = 0; i < N; ++i )
    {
        arg[i] = 1;
        int omp_threads = 1;
        tasks[i] = new Task( std::bind( &work, std::ref( arg[i] ) ), omp_threads );
    }

    for ( int i = 0; i < N; ++i )
    {
        delete tasks[i];
    }
}

/* --------------------------------------------------------------------- */

void doThreading( int N )
{
    static int MAX_THREADS = 256;
    std::unique_ptr<int[]> arg( new int[N] );
    std::unique_ptr<std::thread[]> threads( new std::thread[N] );

    for ( int i = 0; i < N; ++i )
    {
        arg[i] = 1;

        threads[i] = std::thread( work, std::ref( arg[i] ) );

        // if we have more than maximal number of threads we wait for previous ones

        if ( i > MAX_THREADS )
        {
            threads[i - MAX_THREADS].join();
        }
    }
}

/* --------------------------------------------------------------------- */

void doSelf( int N )
{
    std::unique_ptr<int[]> arg( new int[N] );

    for ( int i = 0; i < N; ++i )
    {
        arg[i] = 1;
        work( arg[i] );
    }
}

/* --------------------------------------------------------------------- */

int main()
{
    // static int N = 100000;
    static int N = 5000;
    double time0 = common::Walltime::get();
    doSelf( N );
    time0 = common::Walltime::get() - time0;
    double time1 = common::Walltime::get();
    doThreading( N );
    time1 = common::Walltime::get() - time1;
    double time2 = common::Walltime::get();
    doTasking( N );
    time2 = common::Walltime::get() - time2;
    cout << "Execution of " << N << " work routines." << endl;
    cout << "Time for self (all work by this thread)        = " << time0 << endl;
    cout << "Time for threads (one thread for each work)    = " << time1 << endl;
    cout << "Time for tasks (work scheduled by thread pool) = " << time2 << endl;
}
