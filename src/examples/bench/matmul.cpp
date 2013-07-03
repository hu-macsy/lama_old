/**
 * @file matmul.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Benchmark of matrix multiplication on Host and CUDA
 * @author Thomas Brandes
 * @date 03.06.2013
 * @since 1.0.1
 */

#include <iostream>
#include <iomanip>

#include <lama.hpp>

// Matrix & vector related includes
#include <lama/expression/all.hpp>
#include <lama/matrix/all.hpp>

#include <lama/matutils/MatrixCreator.hpp>
#include <lama/Walltime.hpp>

using namespace lama;
using namespace std;

static bool verboseFlag = false;

template<typename ValueType>
static void bench( IndexType size, double fillRate )
{
    ContextPtr host = ContextFactory::getContext( Context::Host );

    CSRSparseMatrix<ValueType> a( size, size );
    CSRSparseMatrix<ValueType> b( size, size );
    CSRSparseMatrix<ValueType> c( size, size );

    MatrixCreator<ValueType>::fillRandom( a, fillRate );
    MatrixCreator<ValueType>::fillRandom( b, fillRate );
 
    a.setContext( host );
    b.setContext( host );
    c.setContext( host );

    a.prefetch();
    b.prefetch();
    a.wait();
    b.wait();

    double timeHost = Walltime::get();

    c = a * b;

    timeHost = Walltime::get() - timeHost;

    ContextPtr gpu = ContextFactory::getContext( Context::CUDA );

    a.setContext( gpu );
    b.setContext( gpu );

    CSRSparseMatrix<ValueType> c1( size, size );

    a.prefetch();
    b.prefetch();
    a.wait();
    b.wait();

    double timeGPU = Walltime::get();

    c1 = a * b;

    timeGPU = Walltime::get() - timeGPU;

    // check maxDiff

    Scalar maxDiff = c.maxDiffNorm( c1 );

    cout << "max diff Host/GPU matrix = " << maxDiff.getValue<ValueType>() << endl;

    if ( verboseFlag )
    {
        cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        cout << "c = " << c << endl;

        cout << "time <" << Scalar::getType<ValueType>() << "> on " << *gpu
                 << ", size = " << size << ", rate = " << fillRate << " = " << time << endl;
    }

    const int precision = 1;

    cout << "Size = " << size << ", rate = " << ( fillRate * 100 )
         << "%, type = " << Scalar::getType<ValueType>() << endl;
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
    IndexType sizes[] = { 100, 350, 1000, 2500 };
    double fillrates[] = { 0.005, 0.01, 0.02, 0.05 };

    int nsizes = sizeof( sizes ) / sizeof( IndexType );
    int nrates = sizeof( fillrates ) / sizeof( double );

    for ( int i = 0; i < nsizes; ++i )
    {
        IndexType size = sizes[i];

        for ( int j = 0; j < nrates; ++j )
        {
            double rate = fillrates[j];
 
            bench<float>( size, rate );
            bench<double>( size, rate );

        }
    }
}

