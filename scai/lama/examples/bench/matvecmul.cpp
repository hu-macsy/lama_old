/**
 * @file lama/examples/bench/matvecmul.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Benchmark of matrix-vector multiplication on Host and CUDA
 * @author Thomas Brandes
 * @date 03.06.2013
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// Matrix & vector related includes

#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/all.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Settings.hpp>

using namespace scai;
using namespace scai::lama;
using namespace scai::hmemo;
using namespace std;
using scai::common::Walltime;

typedef common::shared_ptr<_DenseVector> DenseVectorPtr;

static void bench( Matrix& mat )
{
    ContextPtr ctx = Context::getContextPtr();

    DenseVectorPtr x( mat.newDenseVector() );
    DenseVectorPtr y1( mat.newDenseVector() );
    DenseVectorPtr y2( mat.newDenseVector() );

    const IndexType size = mat.getNumRows();

    x->setSequence( Scalar( 0 ), Scalar( 0.1 ), size );

    mat.setCommunicationKind( Matrix::SYNCHRONOUS );

    mat.setContextPtr( ctx );
    x->setContextPtr( ctx );
    mat.prefetch();
    x->prefetch();
    mat.wait();
    x->wait();

    cout << "x = " << *x << endl;

    common::unique_ptr<Matrix> matT( mat.newMatrix() );

    double timeT = Walltime::get();
    {
        SCAI_REGION( "Main.Bench.transpose" )
        matT->assignTranspose( mat );
    }
    timeT = Walltime::get() - timeT;
    timeT *= 1000.0;   // scale to ms

    cout << "matT = " << *matT << endl;

    double time1 = Walltime::get();

    {
        SCAI_REGION( "Main.Bench.gemv" )
        *y1 = mat * *x;
    }

    time1 = Walltime::get() - time1;
    time1 *= 1000.0;   // scale to ms

    cout << "y1  = mat * x = " << *y1 << endl;

    double time2 = Walltime::get();
    {
        SCAI_REGION( "Main.Bench.gevm" )
        *y2 = *x * *matT;
    }

    time2 = Walltime::get() - time2;
    time2 *= 1000.0;   // scale to ms

    cout << "transpose: " << timeT << " ms" << endl;
    cout << "gemv: " << time1 << " ms, gevm: " << time2 << " ms" << endl;

    // check result

    *y1 -= *y2;

    cout << "max diff = " << y1->maxNorm() << endl;
}

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.Bench.main" )

    common::Settings::parseArgs( argc, argv );

    COOSparseMatrix<RealType> C( 10000, 10000 );

    MatrixCreator::fillRandom( C, 0.1 );

    cout << "Bench this matrix: " << C << endl;

    {
        SCAI_REGION( "Main.Bench.main" )
        bench( C );
    }
}

