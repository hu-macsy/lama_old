/**
 * @file lsbcJoined.cpp
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
 * @brief Example of least square problem with boundary conditions
 * @author Thomas Brandes, Andreas Borgen Langva
 * @date 21.07.2017
 */

#include <scai/ipbcls/ConstrainedLeastSquares.hpp>

#include <scai/tracing.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/lama/matrix/MatrixWithT.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/common/Settings.hpp>
#include <scai/partitioning/Partitioning.hpp>

#include "JoinedMatrix.hpp"
#include "JoinedVector.hpp"

#include <cstdlib>

using namespace scai;
using namespace lama;
using namespace ipbcls;

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
 
    SCAI_ASSERT_EQ_ERROR( comm->getSize(), 1, "This example program does not run parallel yet" )

    SCAI_REGION( "Main.driver" )

    common::Settings::parseArgs( argc, argv );

    if ( argc < 7 )
    {
        if ( comm->getRank() == 0 )
        {
            std::cout << "Please call " << argv[0]
                      << " <tolerance> <filename_D> <filename_L> <filename_T> <filename_lb> <filename_ub> [ <filename_x> ]" << std::endl
                      << " where <tolerance> is a floating point number between 0 and 1. " << std::endl;
        }

        return -1;
    }

    double tol = std::strtod( argv[1], NULL );

    if ( tol == 0.0 )
    {
        std::cout << "Invalid tolerance supplied." << std::endl;
        return -1;
    }

    auto D  = read<CSRSparseMatrix<ValueType>>( argv[2] );
    auto L  = read<CSRSparseMatrix<ValueType>>( argv[3] );

    // Alternatively if geometry is known.
    // StencilMatrix<ValueType> L( common::Grid3D( 73, 32, 121 ), common::Stencil3D<ValueType>( 7 ) );

    auto T  = read<DenseVector<ValueType>>( argv[4] );
    auto lb = read<DenseVector<ValueType>>( argv[5] );
    auto ub = read<DenseVector<ValueType>>( argv[6] );

    D.setCommunicationKind( SyncKind::SYNCHRONOUS );
    L.setCommunicationKind( SyncKind::SYNCHRONOUS );

    auto colDist = L.getColDistributionPtr();

    auto Zero = denseVectorFill<ValueType>( colDist, 0 );

    MatrixWithT<ValueType> Lopt( L, L );
    MatrixWithT<ValueType> Dopt( D );

    JoinedMatrix<ValueType> A( Dopt, Lopt );
    JoinedVector<ValueType> b( T, Zero );

    lb.redistribute( colDist );
    ub.redistribute( colDist );

    std::cout << "A = " << A << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "lb = " << lb << std::endl;
    std::cout << "ub = " << ub << std::endl;

    DenseVector<ValueType> x( colDist, 0 );

    ConstrainedLeastSquares<ValueType> lsq( A );

    // lsq.setInnerSolverType( InnerSolverType::StandardCG );
    lsq.setInnerSolverType( InnerSolverType::NewtonStepCG );

    // lsq.useTranspose();       // will matrixTimesVector instead ov vectorTimesMatrix

    lsq.setObjectiveTolerance( static_cast<ValueType>( tol ) );
    lsq.setMaxIter( 500 );

    try
    {
        lsq.solve( x, b, lb, ub );
    }
    catch ( common::Exception& ex )
    {
        std::cout << "Caught exception: " << ex.what() << std::endl;
        std::cout << "Stop execution." << std::endl;
        return 1;
    }

    std::unique_ptr<Vector<ValueType>> residual( A.newTargetVector() );
    *residual =  A * x - b;

    std::cout << "res norm = " << residual->l2Norm() << std::endl;

    if ( argc > 7 )
    {
        x.writeToFile( argv[7] );
        std::cout << "written solution to file " << argv[6] << std::endl;
    }
}
