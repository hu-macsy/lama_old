/**
 * @file lama/examples/stencil/stencilTranspose.cpp
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
 * @brief Demo of stencil class and generation of Stencil matrix
 * @author Thomas Brandes
 * @date 23.02.2017
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

typedef DefaultReal ValueType;

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    common::Stencil1D<ValueType> stencilFD8;

    stencilFD8.reserve( 4 );   // just for convenience, not mandatory

    stencilFD8.addPoint( -1, -5 );
    stencilFD8.addPoint( 0, 3  );
    stencilFD8.addPoint( 1, -3 );
    stencilFD8.addPoint( 2, 5 ) ;

    // stencilFD8.addPoint( -3, -5.0/7168.0 );
    // stencilFD8.addPoint( -2, 49.0/5120.0 );
    // stencilFD8.addPoint( -1, -245.0/3072.0 );
    // stencilFD8.addPoint( 0, 1225.0/1024.0 );
    // stencilFD8.addPoint( 1, -1225.0/1024.0 );
    // stencilFD8.addPoint( 2, 245.0/3072.0 ) ;
    // stencilFD8.addPoint( 3, -49.0/5120.0 );
    // stencilFD8.addPoint( 4, 5.0/7168.0 );

    common::Stencil1D<ValueType> stencilBD8;
    stencilBD8.transpose( stencilFD8 );

    common::Grid1D grid1( 5 );
    common::Grid1D grid2( 5 );
 
    grid1.setBorderType( 0, common::BorderType::ABSORBING, common::BorderType::ABSORBING );
    grid2.setBorderType( 0, common::BorderType::ABSORBING, common::BorderType::ABSORBING );

    StencilStorage<ValueType> s1( grid1, stencilFD8 );
    StencilStorage<ValueType> s2( grid2, stencilBD8 );

    auto csr1 = convert<CSRStorage<ValueType>>( s1 );
    csr1.writeToFile( "csr1.mtx" );
    auto csr2 = convert<CSRStorage<ValueType>>( s2 );
    csr2.writeToFile( "csr2.mtx" );
    CSRStorage<ValueType> csrT;
    csrT.assignTranspose( csr1 );
    csrT.writeToFile( "csrT.mtx" );

    ValueType diff = csr2.maxDiffNorm( csrT );
    std::cout << "max diff = " << diff << std::endl;
}
