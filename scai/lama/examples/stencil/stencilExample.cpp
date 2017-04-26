/**
 * @file stencilExample.cpp
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
 * @brief Demo of stencil class and generation of Stencil matrix
 * @author Thomas Brandes
 * @date 23.02.2017
 */

#include <scai/lama.hpp>

#include "Stencil.hpp"
#include "StencilStorage.hpp"
#include "StencilMatrix.hpp"

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/NoCommunicator.hpp>

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

    Stencil3D<double> stencil3( 27 ); // 27-point stencil
    Stencil2D<double> stencil2( 5 );  // 5-point stencil, two dims
    Stencil1D<double> stencil1( 3 );  // 3-point stencil, one dimension

    Stencil1D<double> stencilFD8;

    stencilFD8.reserve( 8 );   // just for convenience, not mandatory

    stencilFD8.addPoint( -3, -5.0/7168.0 );
    stencilFD8.addPoint( -2, 49.0/5120.0 );
    stencilFD8.addPoint( -1, -245.0/3072.0 );
    stencilFD8.addPoint( 0, 1225.0/1024.0 );
    stencilFD8.addPoint( 1, -1225.0/1024.0 );
    stencilFD8.addPoint( 2, 245.0/3072.0 ) ;
    stencilFD8.addPoint( 3, -49.0/5120.0 );
    stencilFD8.addPoint( 4, 5.0/7168.0 );

    Stencil1D<double> stencilDummy( 1 );

    // 1-dimensional stencils can be combined

    Stencil3D<double> stencilX( stencilFD8, stencilDummy, stencilDummy );
    Stencil3D<double> stencilY( stencilDummy, stencilFD8, stencilDummy );
    Stencil3D<double> stencilZ( stencilDummy, stencilDummy, stencilFD8 );

    Stencil2D<double> stencil2_5( stencil1, stencil1 /*, 1 */ ); 
    Stencil2D<double> stencil2_9( 9 );

    Stencil3D<double> stencil3_7( stencil1, stencil1, stencil1 );
    Stencil3D<double> stencil3_19( 19 );
    Stencil3D<double> stencil3_27( 27 );

    for ( IndexType i = 0; i < stencil3_27.nPoints(); ++i )
    {
        int i1, i2, i3;  // relative position
        double val;
        stencil3_27.getPoint( i1, i2, i3, val, i );
    }
 
    // stencil2 == stencil2X 

    dmemo::CommunicatorPtr comm( new dmemo::NoCommunicator() );

    const IndexType N = 300;

    common::Grid3D grid( N, N, N );

    dmemo::DistributionPtr gridDistribution( new GridDistribution( grid, comm ) );

    // StencilStorage<double> stencilStorage( grid, stencil3_7 );

    // std::cout << "stencilStorage " << stencilStorage << std::endl;

    // CSRStorage<double> csrStorage( stencilStorage );

    StencilMatrix<double> stencilMatrix( grid, stencil3_7 );

    std::cout << "stencilMatrix " << stencilMatrix << std::endl;

    stencilMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );

    // stencilMatrix.writeToFile( "stencil.txt" );

    CSRSparseMatrix<double> csrMatrix;

    MatrixCreator::buildPoisson( csrMatrix, 3, 7, N, N, N );

    csrMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );

    std::cout << "csrStencilMatrix " << csrMatrix << std::endl;

    // csrMatrix.writeToFile( "csr.txt" );

    // Scalar s = csrMatrix.maxDiffNorm( stencilMatrix );
    // std::cout << "max diff = " << s << std::endl;

    DenseVector<double> x( stencilMatrix.getColDistributionPtr() );
    x = 1.0;

    DenseVector<double> y1 ( stencilMatrix * x );
    // y1.writeToFile( "y1.txt" );
    DenseVector<double> y2 ( csrMatrix * x );
    // y2.writeToFile( "y2.txt" );

    std::cout << "diff = " << y1.maxDiffNorm( y2 ) << std::endl;
}
