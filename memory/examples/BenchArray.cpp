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

#include <memory/memory.hpp>

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

    LAMAArray<double> X;
    {
        WriteOnlyAccess<double> write( X, n );
        init( write.get(), n, 1.0 );
    }
    {
        ReadAccess<double> read( X );
        sub( res, read.get(), n );
    }
}

void routineLAMA_1( double& res )
{
    LAMAArray<double> X;
    res = 0.0;
}

void routineLAMA_2( double& res, IndexType n )
{
    LAMAArray<double> X;
    {
        WriteOnlyAccess<double> write( X, n );
    }
    res = 0.0;
}

void routineLAMA_3( double& res, IndexType n )
{
    LAMAArray<double> X;
    {
        WriteOnlyAccess<double> write( X, n );
    }
    {
        ReadAccess<double> read( X );
    }
    res = 0.0;
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
    ContextPtr host = Context::getContextPtr( context::Host );
    MemoryPtr hostMem = host->getMemoryPtr();

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

        double ts = ( common::Walltime::get() - time ) * 1000.0;

        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA_1( res );
        }

        double tl1 = ( common::Walltime::get() - time ) * 1000.0;


        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA_2( res, N );
        }

        double tl2 = ( common::Walltime::get() - time ) * 1000.0;

        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA_3( res, N );
        }

        double tl3 = ( common::Walltime::get() - time ) * 1000.0;

        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA( res, N );
        }

        double tl = ( common::Walltime::get() - time ) * 1000.0;


        cout << "Case " << k << ": N = " << N << ", ITER = " << ITER << endl;
        cout << "res = " << res << endl;
        cout << "routineSimple: " << ts << " ms "  << endl;
        cout << "routineLAMA: " << tl << " ms "  << endl;
        cout << "routineLAMA_1: " << tl1 << " ms "  << endl;
        cout << "routineLAMA_2: " << tl2 << " ms "  << endl;
        cout << "routineLAMA_3: " << tl3 << " ms "  << endl;
        cout << endl;
    }
}
