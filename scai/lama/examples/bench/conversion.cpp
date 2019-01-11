/**
 * @file conversion.cpp
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
 * @brief Benchmark of matrix format/type conversions on Host and CUDA
 * @author Thomas Brandes
 * @date 06.06.2013
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
using std::cout;
using std::endl;
using std::setw;
using scai::IndexType;
using scai::common::Walltime;
using scai::common::ContextType;

//static bool verboseFlag = false;

static void bench( _Matrix& b, _Matrix& a )
{
    ContextPtr host = Context::getHostPtr();
    a.setContextPtr( host );
    b.setContextPtr( host );
    a.prefetch();
    a.wait();
    const int nrepeat = 10;
    double timeHost = Walltime::get();

    for ( int k = 0; k < nrepeat; ++k )
    {
        b.assign( a );
    }

    timeHost = Walltime::get() - timeHost;
    ContextPtr gpu = Context::getContextPtr( ContextType::CUDA );
    a.setContextPtr( gpu );
    b.setContextPtr( gpu );
    a.prefetch();
    a.wait();
    double timeGPU = Walltime::get();

    for ( int k = 0; k < nrepeat; ++k )
    {
        b.assign( a );
    }

    timeGPU = Walltime::get() - timeGPU;
    // check maxDiff
    const int precision = 1;
    cout << "Size = " << a.getNumRows() << " x " << a.getNumColumns() << ", nnz = " << a.getNumValues() << endl;
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

    IndexType sizes[] = { 10000, 30000 };
    float fillrates[] = { 0.001, 0.002, 0.003 };
    int nsizes = sizeof( sizes ) / sizeof( IndexType );
    int nrates = sizeof( fillrates ) / sizeof( double );

    for ( int i = 0; i < nsizes; ++i )
    {
        IndexType size = sizes[i];

        for ( int j = 0; j < nrates; ++j )
        {
            // a bit tricky here as we need two different value types
            typedef SCAI_COMMON_FIRST_ARG( SCAI_NUMERIC_TYPES_HOST ) ValueType;
#if SCAI_COMMON_COUNT_NARG( SCAI_ARITHTMETIC_HOST ) == 1
            // only one supported type, so take the same value type
            typedef ValueType ValueType1;
#else
            // take the second supported value type
            typedef SCAI_COMMON_FIRST_ARG( SCAI_COMMON_TAIL( SCAI_NUMERIC_TYPES_HOST ) ) ValueType1;
#endif
            float rate = fillrates[j];

            auto a  = zero<CSRSparseMatrix<ValueType>>( size, size );
            auto a1 = zero<CSRSparseMatrix<ValueType>>( size, size );
            auto b  = zero<ELLSparseMatrix<ValueType1>>( size, size );
            auto c  = zero<JDSSparseMatrix<ValueType1>>( size, size );
            auto d  = zero<CSRSparseMatrix<ValueType1>>( size, size );

            MatrixCreator::fillRandom( a, rate );
            cout << "ELL <-- CSR" << endl;
            bench( b, a );
            cout << "CSR <-- ELL" << endl;
            bench( a1, b );
            // test for same
            auto maxDiff = a.maxDiffNorm( a1 );
            cout << "max diff = " << maxDiff << endl;
            cout << "JDS <-- CSR" << endl;
            bench( c, a );
            cout << "CSR <-- JDS" << endl;
            bench( a1, c );
            maxDiff = a.maxDiffNorm( a1 );
            cout << "max diff = " << maxDiff << endl;
            cout << "CSR <-- CSR" << endl;
            bench( d, a );
            cout << "CSR <-- CSR" << endl;
            bench( a1, d );
            maxDiff = a.maxDiffNorm( a1 );
            cout << "max diff = " << maxDiff << endl;
        }
    }
}

