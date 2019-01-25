/**
 * @file doAll.cpp
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

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>
#include <scai/lama/matrix/MatrixWithT.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>

#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/ipbcls.hpp>

#include "JoinedMatrix.hpp"
#include "JoinedVector.hpp"

using namespace scai;
using namespace lama;

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    SCAI_REGION( "Main.driver" )
    
    common::Settings::parseArgs( argc, argv );

    if ( argc <= 4 )
    {
        std::cout << argv[0] << " <input_D> <input_T> <input_So> <input_hrz> " << std::endl;
        std::cout << argv[0] << " [ <output_x] " << std::endl;
        std::cout << "            [ --SCAI_STRENGTH=<strength> ]" << std::endl;
        std::cout << "            [ --SCAI_VARIATION=<variation> ]" << std::endl;
        return -1;
    }

    std::cout << "Read D from "  << argv[1] << ", T from " << argv[2] 
              << ", So from " << argv[3] << ", hrz from " << argv[4] << std::endl;

    auto D   = read<CSRSparseMatrix<ValueType>>( argv[1] );
    auto T   = read<DenseVector<ValueType>>( argv[2] );
    auto So  = read<DenseVector<ValueType>>( argv[3]  );
    auto hrz = read<DenseVector<ValueType>>( argv[4] );

    std::cout << "D = " << D << std::endl;
    std::cout << "hrz = " << hrz << std::endl;
    std::cout << "So = " << So << std::endl;
    std::cout << "T = " << T << std::endl;

    IndexType ny = hrz.size();
    IndexType n  = So.size();
    IndexType nz = n / ny ;

    SCAI_ASSERT_EQ_ERROR( ny * nz, n , "Illegal factors ny = " << ny << ", nz = " << nz )

    IndexType nray = D.getNumRows();

    SCAI_ASSERT_EQ_ERROR( D.getNumColumns(), n, "D must have #colums equal to problem size " << ny << " x " << nz )
    SCAI_ASSERT_EQ_ERROR( T.size(), D.getNumRows(), "T cannot be rhs for D" )

    int strength = 10;
    int variation = 2;

    common::Settings::getEnvironment( strength, "SCAI_STRENGTH" );
    common::Settings::getEnvironment( variation, "SCAI_VARIATION" );

    std::cout << "Use strength = " << strength << ", variation = " << variation << std::endl;

    common::Stencil2D<ValueType> stencil( 5 ); stencil.scale( - strength / 4.0 );

    common::Grid2D grid( ny, nz );

    grid.setBorderType( 0, common::Grid::BORDER_ABSORBING );
    grid.setBorderType( 1, common::Grid::BORDER_ABSORBING );

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
 
    dmemo::DistributionPtr gridDist( new dmemo::GridDistribution( grid, comm ) );
    dmemo::DistributionPtr rayDist( new dmemo::BlockDistribution( nray, comm ) );
 
    StencilMatrix<ValueType> L( gridDist, stencil );

    MatrixWithT<ValueType> Lopt( L, L );   // tricky stuff for symmetric matrix

    DenseVector<ValueType> Zero( L.getNumRows(), 0 );

    So.redistribute( gridDist );

    DenseVector<ValueType> lb( So );
    DenseVector<ValueType> ub( So );

    lb *= ( 1. - variation / 100. ); //for testing: 0.01;
    ub  *= ( 1. + variation / 100. ); //for testing: 100.0;
 
    for ( IndexType i = 0; i < ny; ++i )
    {
        for ( IndexType j = 0; j < nz; ++j )
        {
            if ( hrz[i] > j + 1 )
            {
                IndexType pos = i * nz + j;
                lb[pos] = So[pos] * 0.9999;
                ub[pos] = So[pos] * 1.0001;
            }
        }
    }


    // take context as specified by SCAI_CONTEXT

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    D.setContextPtr( ctx );
    D.setCommunicationKind( SyncKind::SYNCHRONOUS );
    ub.setContextPtr( ctx );
    lb.setContextPtr( ctx );

    D.redistribute( rayDist, gridDist );
    T.redistribute( rayDist );
    Zero.redistribute( gridDist );

    auto x = denseVectorFill<ValueType>( gridDist, 0 );

    MatrixWithT<ValueType> Dopt( D );

    JoinedMatrix<ValueType> A( Dopt, Lopt );
    JoinedVector<ValueType> T_ext( T, Zero );

    std::cout << "construct lsq." << std::endl;

    ipbcls::ConstrainedLeastSquares<ValueType> lsq( A );

    // lsq.useTranspose();       // will matrixTimesVector instead ov vectorTimesMatrix

    lsq.setObjectiveTolerance( 0.01 );
    lsq.setMaxIter( 50 );

    try
    {
        std::cout << "solve lsq with boundary cond" << std::endl;
        lsq.solve( x, T_ext, lb, ub );
    }
    catch ( common::Exception& ex )
    {
        std::cout << "Caught exception: " << ex.what() << std::endl;
        std::cout << "Stop execution." << std::endl;
        return 1;
    }

    VectorPtr<ValueType> residual( T_ext.newVector() );

    *residual = A * x - T_ext;

    std::cout << "res norm = " << residual->l2Norm() << std::endl;

    if ( argc > 5 )
    {
        x.writeToFile( argv[5] );
        std::cout << "written solution to file " << argv[5] << std::endl;
    }
}
