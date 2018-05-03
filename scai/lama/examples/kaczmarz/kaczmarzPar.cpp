/**
 * @file kaczmarzPar.cpp
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
 * @brief Kaczmarz method with full matrix vector operations
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
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace hmemo;
using namespace lama;
using namespace dmemo;

int main( int argc, const char* argv[] )
{
    if ( argc != 3 )
    {
        std::cout << "Correct call: " << argv[0] << " <filename> <n_iter>" << std::endl;
        return -1;
    }

    std::string filename = argv[1];
    IndexType maxIter = atoi( argv[2] );

    std::cout << "Run " << argv[0] << " " << filename << ", #iter = " << maxIter << std::endl;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    SCAI_REGION( "main.Kacmarz" )

    CSRSparseMatrix<double> a( filename );

    SCAI_ASSERT_EQ_ERROR( a.getNumRows(), a.getNumColumns(), "example only for square matrices" )

    IndexType size = a.getNumRows();

    DistributionPtr blockDist( new BlockDistribution( size, comm ) );
    DistributionPtr repDist( new NoDistribution( size ) );

    a.redistribute( blockDist, repDist );

    std::cout << "Matrix = " << a << std::endl;

    DenseVector<double> rowDotP( blockDist, 0.0 );

    SparseVector<double> row;

    {
        SCAI_REGION( "main.RowDotp" )

        double time = common::Walltime::get();

        const Distribution& rowDist = a.getRowDistribution();

        HArray<double>& localRowDot = rowDotP.getLocalValues();

        for ( IndexType i = 0; i < size; ++i )
        {
            IndexType localI = rowDist.global2local( i );

            if ( localI == invalidIndex )
            {
                continue;
            }

            a.getRowLocal( row, localI );

            Scalar dotP = row.dotProduct( row );

            utilskernel::HArrayUtils::setVal( localRowDot, localI, 1.0 / dotP.getValue<double>() );
        }

        time = common::Walltime::get() - time;

        std::cout << "Row norm calculation took " << time << " seconds." << std::endl;
    }
 
    std::cout << "Norms: min = " << rowDotP.min() << ", max = " << rowDotP.max() << std::endl;

    a.redistribute( blockDist, blockDist );

    DenseVector<double> b( a.getRowDistributionPtr(), 1.0 );
    DenseVector<double> x( a.getColDistributionPtr(), 0.0 );

    std::cout << "built dotp of each row, now start" << std::endl;

    // x = x + ( b(i) - < a(i,:) * x(:)> / 

    const IndexType printIter = 10;

    DenseVector<double> res;

    double omega = 1.0;  // underrelaxation required

    for ( IndexType iter = 0; iter < maxIter; ++iter )
    {
        {
            SCAI_REGION( "main.FullIter" )

            res = b - a * x;
            res *= rowDotP;
            x += omega * res * a;

        }

        if ( ( iter + 1 ) % printIter == 0 )
        {
            SCAI_REGION( "main.Residual" )

            DenseVector<double> res( a * x - b );

            Scalar norm = res.l2Norm();

            if ( comm->getRank() == 0 )
            {
                std::cout << "Iter = " << ( iter + 1 ) << ", res = " << norm << std::endl;
            }
        } 
    }
}
