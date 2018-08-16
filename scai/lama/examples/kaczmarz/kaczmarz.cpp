/**
 * @file kaczmarz.cpp
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

    CSRSparseMatrix<double> a( filename );

    SCAI_ASSERT_EQ_ERROR( a.getNumRows(), a.getNumColumns(), "example only for square matrices" )

    IndexType size = a.getNumRows();

    DistributionPtr blockDist( new BlockDistribution( size, comm ) );
    DistributionPtr repDist( new NoDistribution( size ) );

    a.setContextPtr( ctx );
    a.redistribute( blockDist, blockDist );

    std::cout << "Matrix = " << a << std::endl;

    DenseVector<double> b( blockDist, 1.0, ctx );
    DenseVector<double> x( blockDist, 0.0, ctx );
    DenseVector<double> rowNorm( size, 0.0, ctx );

    SparseVector<double> ai( ctx );

    {
        SCAI_REGION( "main.RowDotp" )

        double time = common::Walltime::get();

        for ( IndexType i = 0; i < size; ++i )
        {
            a.getRow( ai, i );
            Scalar dotP = ai.dotProduct( ai );
            rowNorm[i]  = dotP.getValue<double>();
        }

        time = common::Walltime::get() - time;

        std::cout << "Row norm calculation took " << time << " seconds." << std::endl;
    }

    std::cout << "built dotp of each row, now start" << std::endl;

    // x = x + ( b(i) - < a(i,:) * x(:)> / 

    const IndexType printIter = 10;

    Scalar s1, s2, bi, alpha;

    for ( IndexType iter = 0; iter < maxIter; ++iter )
    {
        {
            SCAI_REGION( "main.FullIter" )

            for ( IndexType i = 0; i < size; ++i )
            {
                a.getRow( ai, i );
                s1 = ai.dotProduct( x );
                s2 = rowNorm[i];
                bi = b[i];
                alpha = ( bi - s1 ) / s2;
                x += alpha * ai;
            }
        }

        if ( ( iter + 1 ) % printIter == 0 )
        {
            SCAI_REGION( "main.Residual" )

            DenseVector<double> y( a * x );
            DenseVector<double> res( y - b );
            std::cout << "Iter = " << ( iter + 1 ) << ", res = " << res.l2Norm() << std::endl;
        } 
    }
}
