/**
 * @file conversion.cpp
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
 * @brief Benchmark of matrix format/type conversions on Host and CUDA
 * @author Thomas Brandes
 * @date 06.06.2013
 * @since 1.0.1
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/common/Walltime.hpp>

using namespace scai::lama;
using namespace std;
using scai::common::Walltime;

static bool verboseFlag = false;

static void bench( Matrix& b, Matrix& a )
{
    ContextPtr host = Context::getContextPtr( context::Host );

    a.setContext( host );
    b.setContext( host );

    a.prefetch();
    a.wait();

    const int nrepeat = 10;

    double timeHost = Walltime::get();

    for ( int k = 0; k < nrepeat; ++k )
    {
        b = a;
    }

    timeHost = Walltime::get() - timeHost;

    ContextPtr gpu = Context::getContextPtr( context::CUDA );

    a.setContext( gpu );
    b.setContext( gpu );

    a.prefetch();
    a.wait();

    double timeGPU = Walltime::get();

    for ( int k = 0; k < nrepeat; ++k )
    {
        b = a;
    }

    timeGPU = Walltime::get() - timeGPU;

    // check maxDiff

    const int precision = 1;

    cout << "Size = " << a.getNumRows() << " x " << a.getNumColumns() << ", nnz = " << a.getNumValues() << endl;
    cout << "===================================" << endl;
    cout << setiosflags( std::ios::fixed ) << std::setprecision( precision );
    cout << "time host = " << setw(6) << timeHost * 1000.0 << endl;
    cout << setiosflags( std::ios::fixed ) << std::setprecision( precision );
    cout << "time cuda = " << setw(6) << timeGPU * 1000.0 << endl;
    cout << setiosflags( std::ios::fixed ) << std::setprecision( precision );
    cout << "speedup   = " << setw(6) << ( timeHost / timeGPU ) << endl;;
    cout << endl;
}

int main()
{
    IndexType sizes[] = { 10000, 30000 };
    double fillrates[] = { 0.001, 0.002, 0.003 };

    int nsizes = sizeof( sizes ) / sizeof( IndexType );
    int nrates = sizeof( fillrates ) / sizeof( double );

    for ( int i = 0; i < nsizes; ++i )
    {
        IndexType size = sizes[i];

        for ( int j = 0; j < nrates; ++j )
        {
            typedef float ValueType;
            typedef double ValueType1;

            double rate = fillrates[j];
 
            CSRSparseMatrix<ValueType> a( size, size );
            CSRSparseMatrix<ValueType> a1( size, size );

            ELLSparseMatrix<ValueType1> b( size, size );
            JDSSparseMatrix<ValueType1> c( size, size );
            CSRSparseMatrix<ValueType1> d( size, size );

            MatrixCreator<ValueType>::fillRandom( a, rate );

            cout << "ELL <-- CSR" << endl;
            bench( b, a );
            cout << "CSR <-- ELL" << endl;
            bench( a1, b );

            // test for same

            Scalar maxDiff = a.maxDiffNorm( a1 );
            cout << "max diff = " << maxDiff.getValue<ValueType>() << endl;

            cout << "JDS <-- CSR" << endl;
            bench( c, a );
            cout << "CSR <-- JDS" << endl;
            bench( a1, c );

            maxDiff = a.maxDiffNorm( a1 );
            cout << "max diff = " << maxDiff.getValue<ValueType>() << endl;

            cout << "CSR <-- CSR" << endl;
            bench( d, a );
            cout << "CSR <-- CSR" << endl;
            bench( a1, d );

            maxDiff = a.maxDiffNorm( a1 );
            cout << "max diff = " << maxDiff.getValue<ValueType>() << endl;
        }
    }
}

