/**
 * @file lama/examples/kaczmarz/solver.cpp
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
 * @brief Kaczmarz method with sparse vector
 * @author Thomas Brandes
 * @date 23.02.2017
 */

#include <scai/lama.hpp>

// _Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai::hmemo;
using namespace scai::lama;
using namespace scai::dmemo;

int main()
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    SCAI_REGION( "main.Kacmarz" )

    CSRSparseMatrix<double> matrix( "matrix.frm" );

    SCAI_ASSERT_EQ_ERROR( matrix.getNumRows(), matrix.getNumColumns(), "example only for square matrices" )

    IndexType size = matrix.getNumRows();

    DistributionPtr blockDist( new BlockDistribution( size, comm ) );
    DistributionPtr repDist( new NoDistribution( size ) );

    matrix.redistribute( blockDist, blockDist );

    std::cout << "Matrix = " << matrix << std::endl;

    DenseVector<double> b( blockDist, 1.0 );
    DenseVector<double> x( blockDist, 0.0 );
    DenseVector<double> rowDotP( size, 0.0 );

    DenseVector<double> row;

    {
        SCAI_REGION( "main.RowDotp" )

        for ( IndexType i = 0; i < size; ++i )
        {
            matrix.getRow( row, i );
            Scalar dotP = row.dotProduct( row );
            rowDotP[i]  = dotP.getValue<double>();
        }
    }

    std::cout << "built dotp of each row, now start" << std::endl;

    // x = x + ( b(i) - < matrix(i,:) * x(:)> / 

    const IndexType maxIter = 1;

    Scalar s1, s2, bi, alpha;

    for ( IndexType iter = 0; iter < maxIter; ++iter )
    {
        {
            SCAI_REGION( "main.FullIter" )

            for ( IndexType i = 0; i < size; ++i )
            {
                matrix.getRow( row, i );
                s1 = row.dotProduct( x );
                s2 = rowDotP[i];
                bi = b[i];
                alpha = ( bi - s1 ) / s2;
                if ( i > 650 && i < 700 )
                {
                    std::cout << "Iter = " << i << " update: alpha = " << alpha << ", bi = " << bi 
                              << ", s1 = " << s1 << ", s2 = " << s2 << std::endl;
                }
                x += alpha * row;
            }
        }

        std::cout << "x norm = " << x.l2Norm() << std::endl;

        // if ( ( iter + 1 ) % 50 == 0 )
        {
            SCAI_REGION( "main.Residual" )

            DenseVector<double> y( matrix * x );
            DenseVector<double> res( y - b );
            std::cout << "Iter = " << ( iter + 1 ) << ", res = " << res.l2Norm() << std::endl;
        } 
    }
}
