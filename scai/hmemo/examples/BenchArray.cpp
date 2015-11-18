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

#include <scai/hmemo.hpp>

#include <scai/common/Walltime.hpp>

#include <iostream>
#include <scai/logging.hpp>

using namespace std;
using namespace scai::hmemo;

SCAI_LOG_DEF_LOGGER( logger, "BenchArray" )

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
    SCAI_LOG_TRACE( logger, "routineLAMA, n = " << n )

    ContextPtr contextPtr = Context::getHostPtr();

    LAMAArray<double> X( contextPtr );

    {
        WriteOnlyAccess<double> write( X, contextPtr, n );
        init( write.get(), n, 1.0 );
    }
    {
        ReadAccess<double> read( X, contextPtr );
        sub( res, read.get(), n );
    }
}

void routineSCAI_1( double& res )
{
    LAMAArray<double> X;
    res = 0.0;
}

void routineSCAI_2( double& res, IndexType n )
{
    ContextPtr contextPtr = Context::getHostPtr();

    LAMAArray<double> X;

    {
        WriteOnlyAccess<double> write( X, contextPtr, n );
    }
    res = 0.0;
}

void routineSCAI_3( double& res, IndexType n )
{
    ContextPtr contextPtr = Context::getHostPtr();

    LAMAArray<double> X;

    {
        WriteOnlyAccess<double> write( X, contextPtr, n );
    }
    {
        ReadAccess<double> read( X, contextPtr );
    }
    res = 0.0;
}

void routineSimple( double& res, IndexType n )
{
    SCAI_LOG_TRACE( logger, "routineSimple, n = " << n )
    double* X = new double[n];
    init ( X, n, 1.0 );
    sub ( res, X, n );
    delete [] X;
}

int main()
{
    ContextPtr host = Context::getHostPtr();
    MemoryPtr hostMem = host->getMemoryPtr();

    static int ITER_VEC[] = { 1000000, 100000, 10000, 1000, 100, 10, 1 };
    static int N_VEC[]    = { 1, 10, 100, 1000, 10000, 100000, 1000000 };

    int NCASES = sizeof( ITER_VEC ) / sizeof( int );

    for ( int k = 0; k < NCASES; ++k )
    {
        int ITER = ITER_VEC[k];
        int N    = N_VEC[k];
        double res = 0.0;  // avoids dead code elimination
        double time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSimple( res, N );
        }

        double ts = ( scai::common::Walltime::get() - time ) * 1000.0;

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_1( res );
        }

        double tl1 = ( scai::common::Walltime::get() - time ) * 1000.0;


        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_2( res, N );
        }

        double tl2 = ( scai::common::Walltime::get() - time ) * 1000.0;

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_3( res, N );
        }

        double tl3 = ( scai::common::Walltime::get() - time ) * 1000.0;

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA( res, N );
        }

        double tl = ( scai::common::Walltime::get() - time ) * 1000.0;


        cout << "Case " << k << ": N = " << N << ", ITER = " << ITER << endl;
        cout << "res = " << res << endl;
        cout << "routineSimple (malloc/write/read)                : " << ts << " ms "  << endl;
        cout << "routineLAMA (LAMA array write/read)              : " << tl << " ms "  << endl;
        cout << "routineSCAI_1 (constructor LAMA array, size = 0) : " << tl1 << " ms "  << endl;
        cout << "routineSCAI_2 (constructor LAMA array, size = N) : " << tl2 << " ms "  << endl;
        cout << "routineSCAI_3 (as before, but also read access)  : " << tl3 << " ms "  << endl;
        cout << "LAMA overhead : " << ( ( tl - ts ) / ITER * 1000.0 ) << " us "  << endl;
        cout << endl;
    }
}
