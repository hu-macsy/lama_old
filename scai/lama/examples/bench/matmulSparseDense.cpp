/**
 * @file matmul.cpp
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
 * @brief Benchmark of matrix multiplication on Host and CUDA
 * @author Thomas Brandes
 * @date 03.06.2013
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

using namespace scai;
using namespace lama;
using utilskernel::HArrayUtils;

int main()
{
    auto A = read<CSRSparseMatrix<double>>( "inputA.mtx" );

    const IndexType m = A.getNumRows();
    const IndexType k = A.getNumColumns();

    const IndexType n = 10;  // number of rhs

    HArray<double> denseData( k * n, 1.0 );
    DenseMatrix<double>B( DenseStorage( k, n, std::move( denseData ) ) );

    B.writeToFile( "B.mtx" );

    // sparse Matrix * dense Matrix

    DenseMatrix<double> M;
    M = A * B;

    M.writeToFile( "AB.mtx" );

    std::cout << "M = " << M << std::endl;

    DenseMatrix<double> M1( m, n );

    for ( IndexType i = 0; i < n; ++i )
    {
        DenseVector<double> x;
        B.getColumn( x, i );
        DenseVector<double> y;
        y = A * x;
        M1.setColumn( y, i, common::BinaryOp::COPY );
    }

    std::cout << "M1 = " << M1 << std::endl;
    M1.redistribute( M.getRowDistributionPtr(), M.getColDistributionPtr() );
    std::cout << "M1 = " << M1 << std::endl;

    M = M - M1;
    std::cout << "max diff = " << M.maxNorm() << std::endl;
}

