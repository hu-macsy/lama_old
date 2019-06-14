/**
 * @file stencilDist.cpp
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
 * @brief Demo of distributed stencil
 * @author Thomas Brandes
 * @date 27.04.2017
 */

#include <scai/lama.hpp>


// _Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>
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

// DefaultReal: double or float

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device
    //   SCAI_NP      = ...    set the default processor grid for grid distribution

    common::Settings::parseArgs( argc, argv );

    common::Stencil1D<ValueType> stencilFD8;

    stencilFD8.reserve( 8 );   // just for convenience, not mandatory

    stencilFD8.addPoint( -3, -5.0/7168.0 );
    stencilFD8.addPoint( -2, 49.0/5120.0 );
    stencilFD8.addPoint( -1, -245.0/3072.0 );
    stencilFD8.addPoint( 0, 1225.0/1024.0 );
    stencilFD8.addPoint( 1, -1225.0/1024.0 );
    stencilFD8.addPoint( 2, 245.0/3072.0 ) ;
    stencilFD8.addPoint( 3, -49.0/5120.0 );
    stencilFD8.addPoint( 4, 5.0/7168.0 );

    common::Stencil1D<ValueType> stencilDummy( 1 );

    common::Stencil2D<ValueType> stencilX( stencilFD8, stencilDummy );

    common::Stencil3D<ValueType> stencil( 7 );

    const IndexType N = 4; 

    common::Grid3D grid( N, N, N );

    // take default communicator and default processor array for mapping of the grid

    auto gridDistribution = dmemo::gridDistribution( grid );

    StencilMatrix<ValueType> distStencilMatrix( gridDistribution, stencil );
    CSRSparseMatrix<ValueType> distCSRMatrix( distStencilMatrix );

    std::cout << "distributed stencilMatrix " << distStencilMatrix << std::endl;

    StencilMatrix<ValueType> repStencilMatrix( grid, stencil );

    std::cout << "replicated stencilMatrix " << repStencilMatrix << std::endl;

    DenseVector<ValueType> repX;

    repX.setRandom( repStencilMatrix.getColDistributionPtr(), 1 );

    DenseVector<ValueType> distX( repX );
    distX.redistribute( distStencilMatrix.getColDistributionPtr() );

    std::cout << "max diff X = " << repX.maxDiffNorm( distX ) << std::endl;

    // replicated and distributed matrix-vector multiplication

    const auto repY  = denseVectorEval( repStencilMatrix * repX );
    const auto distY = denseVectorEval( distStencilMatrix * distX );

    std::cout << "max diff Y = " << repY.maxDiffNorm( distY ) << std::endl;
}
