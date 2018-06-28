/**
 * @file matmul.cpp
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
 * @brief Benchmark of matrix multiplication on Host and CUDA
 * @author Thomas Brandes
 * @date 03.06.2013
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// _Matrix & vector related includes
#include <scai/lama/matrix/all.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/common/Walltime.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace std;
using scai::IndexType;
using scai::common::Walltime;
using scai::common::ContextType;

static bool verboseFlag = false;

template<typename ValueType>
static void bench( IndexType size, float fillRate )
{
    ContextPtr host = Context::getHostPtr();

    auto a = zero<CSRSparseMatrix<ValueType>>( size, size );
    auto b = zero<CSRSparseMatrix<ValueType>>( size, size );
    auto c = zero<CSRSparseMatrix<ValueType>>( size, size );

    MatrixCreator::fillRandom( a, fillRate );
    MatrixCreator::fillRandom( b, fillRate );

    a.setContextPtr( host );
    b.setContextPtr( host );
    c.setContextPtr( host );
    a.prefetch();
    b.prefetch();
    a.wait();
    b.wait();
    double timeHost = Walltime::get();
    c = a * b;
    timeHost = Walltime::get() - timeHost;
    ContextPtr gpu = Context::getContextPtr( ContextType::CUDA );
    a.setContextPtr( gpu );
    b.setContextPtr( gpu );

    auto c1 = zero<CSRSparseMatrix<ValueType>>( size, size );

    a.prefetch();
    b.prefetch();
    a.wait();
    b.wait();
    double timeGPU = Walltime::get();
    c1 = a * b;
    timeGPU = Walltime::get() - timeGPU;
    // check maxDiff
    auto maxDiff = c.maxDiffNorm( c1 );
    cout << "max diff Host/GPU matrix = " << maxDiff << endl;

    if ( verboseFlag )
    {
        cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        cout << "c = " << c << endl;
        cout << "time <" << scai::common::getScalarType<ValueType>() << "> on " << *gpu
             << ", size = " << size << ", rate = " << fillRate << endl;
    }

    const int precision = 1;

    cout << "Size = " << size << ", rate = " << ( fillRate * 100 )
         << "%, type = " << scai::common::getScalarType<ValueType>() << endl;

    cout << "===================================" << endl;

    cout << setiosflags( std::ios::fixed ) << std::setprecision( precision );

    cout << "time host = " << setw( 6 ) << timeHost * 1000.0 << endl;

    cout << setiosflags( std::ios::fixed ) << std::setprecision( precision );

    cout << "time cuda = " << setw( 6 ) << timeGPU * 1000.0 << endl;

    cout << setiosflags( std::ios::fixed ) << std::setprecision( precision );

    cout << "speedup   = " << setw( 6 ) << ( timeHost / timeGPU ) << endl;;

    cout << endl;
}

int main()
{
    if ( !Context::hasContext( ContextType::CUDA ) )
    {
        cout << "This examples compares the Host and CUDA implementation. You build without CUDA, so it's skipped." << endl;
        return 0;
    }

    IndexType sizes[] = { 100, 350, 1000, 2500 };
    float fillrates[] = { 0.005, 0.01, 0.02, 0.05 };
    int nsizes = sizeof( sizes ) / sizeof( IndexType );
    int nrates = sizeof( fillrates ) / sizeof( float );

    for ( int i = 0; i < nsizes; ++i )
    {
        for ( int j = 0; j < nrates; ++j )
        {
#define     DO_BENCH( ValueType ) bench<ValueType>( sizes[i], fillrates[j] );
            // do the benchmark for each supported CUDA type
            SCAI_COMMON_LOOP( DO_BENCH, SCAI_NUMERIC_TYPES_CUDA )
        }
    }
}

