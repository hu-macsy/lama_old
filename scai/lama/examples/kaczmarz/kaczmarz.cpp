/**
 * @file kaczmarz.cpp
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
 * @brief Kaczmarz method with sparse vector
 * @author Thomas Brandes
 * @date 23.02.2017
 */

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

// import common 
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <iostream>
#include <stdlib.h>

using namespace scai;
using namespace hmemo;
using namespace lama;
using namespace dmemo;


int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    if ( argc != 3 )
    {
        std::cout << "Correct call: " << argv[0] << " <filename> <n_iter>" << std::endl;
        return -1;
    }

    std::string filename = argv[1];
    IndexType maxIter = atoi( argv[2] );

    std::cout << "Run " << argv[0] << " " << filename << ", #iter = " << maxIter << std::endl;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    ContextPtr ctx = Context::getContextPtr();

    // planned: SCAI_CONTEXT_ACCESS( ctx );   // becomes default context for all created objects

    SCAI_REGION( "main.Kacmarz" )

    CSRSparseMatrix<double> matrix( filename );

    SCAI_ASSERT_EQ_ERROR( matrix.getNumRows(), matrix.getNumColumns(), "example only for square matrices" )

    IndexType size = matrix.getNumRows();

    DistributionPtr blockDist( new BlockDistribution( size, comm ) );
    DistributionPtr repDist( new NoDistribution( size ) );

    matrix.setContextPtr( ctx );
    matrix.redistribute( blockDist, blockDist );

    std::cout << "Matrix = " << matrix << std::endl;

    DenseVector<double> b( blockDist, 1.0, ctx );
    DenseVector<double> x( blockDist, 0.0, ctx );
    DenseVector<double> rowDotP( size, 0.0, ctx );

    SparseVector<double> row( ctx );

    {
        SCAI_REGION( "main.RowDotp" )

        double time = common::Walltime::get();

        for ( IndexType i = 0; i < size; ++i )
        {
            matrix.getRow( row, i );
            Scalar dotP = row.dotProduct( row );
            rowDotP[i]  = dotP.getValue<double>();
        }

        time = common::Walltime::get() - time;

        std::cout << "Row norm calculation took " << time << " seconds." << std::endl;
    }

    std::cout << "built dotp of each row, now start" << std::endl;

    // x = x + ( b(i) - < matrix(i,:) * x(:)> / 

    const IndexType printIter = 1;

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
                x += alpha * row;
            }
        }

        if ( ( iter + 1 ) % printIter == 0 )
        {
            SCAI_REGION( "main.Residual" )

            DenseVector<double> y( matrix * x );
            DenseVector<double> res( y - b );
            std::cout << "Iter = " << ( iter + 1 ) << ", res = " << res.l2Norm() << std::endl;
        } 
    }
}
