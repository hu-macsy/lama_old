/**
 * @file lama/examples/bench/matvecmul.cpp
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
 * @brief Benchmark of matrix-vector multiplication on Host and CUDA
 * @author Thomas Brandes
 * @date 03.06.2013
 */

#include <iostream>
#include <iomanip>
#include <memory>

#include <scai/lama.hpp>

// _Matrix & vector related includes

#include <scai/lama/matrix/all.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

using namespace scai;
using namespace scai::lama;
using namespace scai::hmemo;
using namespace std;
using scai::common::Walltime;

template<typename ValueType>
static void bench( Matrix<ValueType>& mat )
{
    ContextPtr ctx = Context::getContextPtr();

    DenseVector<ValueType> x( ctx );
    DenseVector<ValueType> y1( ctx );
    DenseVector<ValueType> y2( ctx );

    x.allocate( mat.getRowDistributionPtr() );
    y1.allocate( mat.getRowDistributionPtr() );
    y2.allocate( mat.getRowDistributionPtr() );

    const IndexType size = mat.getNumRows();
    const IndexType bound = 1; 

    x.setRandom( size, bound );

    mat.setCommunicationKind( SyncKind::SYNCHRONOUS );

    mat.setContextPtr( ctx );
    x.setContextPtr( ctx );
    mat.prefetch();
    x.prefetch();
    mat.wait();
    x.wait();

    cout << "x = " << x << endl;

    std::unique_ptr<Matrix<ValueType> > matT( mat.newMatrix() );

    double timeT = Walltime::get();
    {
        SCAI_REGION( "Main.Bench.transpose" )
        *matT = transpose( mat );
    }
    timeT = Walltime::get() - timeT;
    timeT *= 1000.0;   // scale to ms

    cout << "matT = " << *matT << endl;

    double time1 = Walltime::get();

    {
        SCAI_REGION( "Main.Bench.gemv_normal" )
        y1 = mat * x;
    }

    time1 = Walltime::get() - time1;
    time1 *= 1000.0;   // scale to ms

    cout << "y1  = mat * x = " << y1 << endl;

    double time2 = Walltime::get();
    {
        SCAI_REGION( "Main.Bench.gemv_transpose" )
        y2 = transpose( *matT ) * x;
    }

    time2 = Walltime::get() - time2;
    time2 *= 1000.0;   // scale to ms

    cout << "Benchmark results: GEMV, normal vs transpose" << endl;
    cout << "============================================" << endl;
    cout << "gemv normal: " << time1 << " ms, gemv transpose: " << time2 << " ms" << endl;
    cout << "explicit transpose: " << timeT << " ms" << endl;

    // check result

    y1 -= y2;

    cout << "max diff = " << y1.maxNorm() << endl;
}

/** Benchmark driver program to measure performance of op( matrix ) * vector (gemv) with op = NORMAL or TRANSPOSE
 *
 *  The benchmark creates a sparse matrix of size N x N and fills it randomly.
 *  
 */
int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.Bench.main" )

    common::Settings::parseArgs( argc, argv );

    const IndexType N = 10000;

    const float fillRate = 0.1;

    auto C = zero<CSRSparseMatrix<DefaultReal>>( N, N );

    MatrixCreator::fillRandom( C, fillRate );

    cout << "Bench this matrix: " << C << endl;

    {
        SCAI_REGION( "Main.Bench.main" )
        bench( C );
    }
}

