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
#include <vector>
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

void sumIt( double& res, const double data[], IndexType size )
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

    HArray<double> X( contextPtr );

    {
        WriteOnlyAccess<double> write( X, contextPtr, n );
        init( write.get(), n, 1.0 );
    }
    {
        ReadAccess<double> read( X, contextPtr );
        sumIt( res, read.get(), n );
    }
}

void routineSCAI_1( double& res )
{
    HArray<double> X;
    res = 0.0;
}

void routineSCAI_2( double& res, IndexType n )
{
    ContextPtr contextPtr = Context::getHostPtr();

    HArray<double> X;

    {
        WriteOnlyAccess<double> write( X, contextPtr, n );
    }
    res = 0.0;
}

void routineSCAI_3( double& res, IndexType n )
{
    ContextPtr contextPtr = Context::getHostPtr();

    HArray<double> X;

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
    double* X = new double[n];
    init ( X, n, 1.0 );
    sumIt ( res, X, n );
    delete [] X;
}

void routineVector( double& res, IndexType n )
{
    std::vector<double> X( n );
    init ( &X[0], n, 1.0 );
    sumIt ( res, &X[0], n );
}

int main()
{
    ContextPtr host = Context::getHostPtr();
    MemoryPtr hostMem = host->getMemoryPtr();

    static int ITER_VEC[] = { 1000000, 100000, 20000, 5000,  1000,    300,     100 };
    static int N_VEC[]    = {       1,     10,   100, 1000, 10000, 100000, 1000000 };

    int NCASES = sizeof( ITER_VEC ) / sizeof( int );

    for ( int k = 0; k < NCASES; ++k )
    {
        int ITER = ITER_VEC[k];
        int N    = N_VEC[k];
        double res = 0.0;  // avoids dead code elimination

        routineLAMA( res, N ); // warm up
        routineLAMA( res, N ); // warm up

        double time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA( res, N );
        }

        double tl = ( scai::common::Walltime::get() - time );
        tl *= 1000.0 * 1000.0 / ITER;

        routineVector( res, N );  // warm up
        routineVector( res, N );  // warm up

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineVector( res, N );
        }

        double tv = ( scai::common::Walltime::get() - time );
        tv *= 1000.0 * 1000.0 / ITER;

        routineSimple( res, N );   // warm up
        routineSimple( res, N );   // warm up

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSimple( res, N );
        }

        double ts = ( scai::common::Walltime::get() - time );
        ts *= 1000.0 * 1000.0 / ITER;

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_1( res );
        }

        double tl1 = ( scai::common::Walltime::get() - time );
        tl1 *= 1000.0 * 1000.0 / ITER;

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_2( res, N );
        }

        double tl2 = ( scai::common::Walltime::get() - time );
        tl2 *= 1000.0 * 1000.0 / ITER;

        time = scai::common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_3( res, N );
        }

        double tl3 = ( scai::common::Walltime::get() - time );
        tl3 *= 1000.0 * 1000.0 / ITER;

        cout << "Case " << k << ": N = " << N << ", ITER = " << ITER << endl;
        cout << "res = " << res << endl;
        cout << "routineSimple (malloc/write/read)                : " << ts << " µs "  << endl;
        cout << "routineVector (vector/write/read)                : " << tv << " µs "  << endl;
        cout << "routineLAMA (LAMA array write/read)              : " << tl << " µs "  << endl;
        cout << "routineSCAI_1 (constructor LAMA array, size = 0) : " << tl1 << " µs "  << endl;
        cout << "routineSCAI_2 (constructor LAMA array, size = N) : " << tl2 << " µs "  << endl;
        cout << "routineSCAI_3 (as before, but also read access)  : " << tl3 << " µs "  << endl;
        cout << "LAMA overhead : " << ( ( tl - ts ) / ITER * 1000.0 ) << " us "  << endl;
        cout << endl;
    }
}
