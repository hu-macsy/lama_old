/**
 * @file BenchArray.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Benchmarking of operations on LAMA arrays
 * @author Thomas Brandes
 * @date 03.07.2015
 */

#include <scai/hmemo.hpp>

#include <scai/common/Walltime.hpp>

#include <iostream>
#include <vector>
#include <scai/logging.hpp>

using namespace scai;
using namespace hmemo;

using std::cout;
using std::endl;

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
        double time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineLAMA( res, N );
        }

        double tl = ( common::Walltime::get() - time );
        tl *= 1000.0 * 1000.0 / ITER;
        routineVector( res, N );  // warm up
        routineVector( res, N );  // warm up
        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineVector( res, N );
        }

        double tv = ( common::Walltime::get() - time );
        tv *= 1000.0 * 1000.0 / ITER;
        routineSimple( res, N );   // warm up
        routineSimple( res, N );   // warm up
        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSimple( res, N );
        }

        double ts = ( common::Walltime::get() - time );
        ts *= 1000.0 * 1000.0 / ITER;
        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_1( res );
        }

        double tl1 = ( common::Walltime::get() - time );
        tl1 *= 1000.0 * 1000.0 / ITER;
        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_2( res, N );
        }

        double tl2 = ( common::Walltime::get() - time );
        tl2 *= 1000.0 * 1000.0 / ITER;
        time = common::Walltime::get();

        for ( int i = 0; i < ITER; ++i )
        {
            routineSCAI_3( res, N );
        }

        double tl3 = ( common::Walltime::get() - time );
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
