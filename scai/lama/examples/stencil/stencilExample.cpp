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

// _Matrix & vector related includes
#include <scai/lama/GridVector.hpp>
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

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

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

    // 1-dimensional stencils can be combined

    common::Stencil3D<ValueType> stencilY( stencilDummy, stencilFD8, stencilDummy );

    common::Grid3D grid( 20, 20, 20 );

    StencilMatrix<ValueType> m( grid, stencilY );

    StencilMatrix<ValueType> copyM1( m );

    StencilMatrix<ValueType> copyM2;
    copyM2 = m;

    m.scale( 0.3 );

    CSRSparseMatrix<ValueType> m1;

    m1.assignTranspose( m );
    m1.writeToFile( "m1.mtx" );

    StencilMatrix<ValueType> m2;

    m2.assignTranspose( m );
    m2.writeToFile( "m2.mtx" );

    GridVector<ValueType> v( grid, 1.0 );

    GridVector<ValueType> v1;
    v1 = m * v;
    GridVector<ValueType> v2;
    v2 = transpose( m2 ) * v;

    std::cout << "max diff = " << v1.maxDiffNorm( v2 ) << std::endl;
}
