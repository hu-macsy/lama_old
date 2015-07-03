/**
 * @file BenchArray.cpp
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
 * @brief Benchmarking of operations on LAMA arrays
 * @author: Thomas Brandes
 * @date 03.07.2015
 **/

#include <memory/LAMAArray.hpp>
#include <memory/HostReadAccess.hpp>

#include <common/Walltime.hpp>

#include <iostream>
#include <logging/logging.hpp>

using namespace std;
using namespace memory;

LAMA_LOG_DEF_LOGGER( logger, "BenchArray" )

void init( double X[], IndexType size, double val )
{
    for ( IndexType i = 0; i < size; ++i )
    {
        X[i] = val;
    }
}

void sub( double& res, const double data[], IndexType size )
{
    for ( IndexType i = 0; i < size; ++i )
    {
        res += data[i];
    }
}

void routineLAMA( double& res, IndexType n )
{
    LAMA_LOG_TRACE( logger, "routineLAMA, n = " << n )
    LAMAArray<double> X( n, 1.0 );
    HostReadAccess<double> read( X );
    sub( res, read.get(), n );
}

void routineLAMA1( double& res, IndexType n )
{
    LAMAArray<double> X( n );
}

void routineSimple( double& res, IndexType n )
{
    LAMA_LOG_TRACE( logger, "routineSimple, n = " << n )
    double* X = new double[n];
    init ( X, n, 1.0 );
    sub ( res, X, n );
    delete [] X;
}

int main()
{
    static int ITER_VEC[] = { 1000000, 100000, 10000, 1000, 100, 10, 1 };
    static int N_VEC[]    = { 1, 10, 100, 1000, 10000, 100000, 1000000 };

    int NCASES = sizeof( ITER_VEC ) / sizeof( int );

    for ( int k = 0; k < NCASES; ++k )
    {
        int ITER = ITER_VEC[k];
        int N    = N_VEC[k];
        double res = 0.0;  // avoids dead code elimination
        double time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSimple( res, N );
        }

        double time1 = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA1( res, N );
        }

        double time2 = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA( res, N );
        }

        double time3 = common::Walltime::get();

        double ts = ( time1 - time ) * 1000.0;
        double tl1 = ( time2 - time1 ) * 1000.0;
        double tl = ( time3 - time2 ) * 1000.0;

        cout << "Case " << k << ": N = " << N << ", ITER = " << ITER << endl;
        cout << "res = " << res << endl;
        cout << "routineSimple: " << ts << " ms "  << endl;
        cout << "routineLAMA1: " << tl1 << " ms "  << endl;
        cout << "routineLAMA: " << tl << " ms "  << endl;
        cout << endl;
    }
}
