/**
 * @file BenchTasking.cpp
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
 * @brief Compare performance of Task vs. Thread
 * @author Thomas Brandes
 * @date 02.02.2012
 * @since 1.0.0
 */

#include <scai/tasking/Task.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Thread.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/throw.hpp>

#include <scai/common/bind.hpp>

using namespace std;
using namespace scai;
using namespace scai::tasking;

static const int WORKLOAD = 200;

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
    common::scoped_array<int> arg( new int[N] );
    common::scoped_array<Task*> tasks( new Task*[N] );

    for ( int i = 0; i < N; ++i )
    {
        arg[i] = 1;
        int omp_threads = 1;

        tasks[i] = new Task( common::bind( &work, common::ref( arg[i] )), omp_threads );
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

    using common::Thread;

    common::scoped_array<int> arg( new int[N] );
    common::scoped_array<Thread*> threads( new Thread*[N] );

    for ( int i = 0; i < N; ++i )
    {
        arg[i] = 1;
        threads[i] = new Thread( &work, arg[i] );

        if ( i > MAX_THREADS )
        {
            threads[i - MAX_THREADS]->join();
        }
    }

    for ( int i = 0; i < N; ++i )
    {
        delete threads[i];
    }
}

/* --------------------------------------------------------------------- */

void doSelf( int N )
{
    common::scoped_array<int> arg( new int[N] );

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
    static int N = 100000;

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
    cout << "Time for self    = " << time0 << endl;
    cout << "Time for threads = " << time1 << endl;
    cout << "Time for tasks   = " << time2 << endl;
}
