/**
 * @file stencilGEMV.cpp
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
 * @brief Benchmark to compare stencil GEMV with CSR GEMV.
 * @author Thomas Brandes
 * @date 26.04.2017
 */

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/expression/all.hpp>
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


int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    // Take a default stencil

    common::Stencil1D<double> stencil1( 3 );
    common::Stencil4D<double> stencil( stencil1, stencil1, stencil1, stencil1 );

    // Define a grid with same number of dimensions as stencil

    const IndexType N1 = 40;
    const IndexType N2 = 40;
    const IndexType N3 = 50;
    const IndexType N4 = 50;

    common::Grid4D grid( N1, N2, N3, N4 );

    // Distibute grid onto default processor array, can be set by --SCAI_NP=2x3x2

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    ContextPtr ctx = Context::getContextPtr();

    dmemo::DistributionPtr gridDistribution( new GridDistribution( grid, comm ) );

    // The stencil matrix just needs the grid distribution and the stencil

    StencilMatrix<double> stencilMatrix( gridDistribution, stencil );
    stencilMatrix.setContextPtr( ctx );

    std::cout << "stencilMatrix " << stencilMatrix << std::endl;

    stencilMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );

    // Conversion stencil matrix -> CSR matrix is full supported

    CSRSparseMatrix<double> csrMatrix( stencilMatrix );

    csrMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );
    csrMatrix.setContextPtr( ctx );

    std::cout << "csrStencilMatrix " << csrMatrix << std::endl;

    DenseVector<double> x( stencilMatrix.getColDistributionPtr() );

    x = 1.0;

    const IndexType NITER = 20;

    DenseVector<double> y1 ( stencilMatrix.getRowDistributionPtr(), 0.0 );
    DenseVector<double> y2 ( csrMatrix.getRowDistributionPtr(), 0.0 );

    {
        SCAI_REGION( "main.stencilGEMV" )

        for ( IndexType iter = 0; iter < NITER; ++iter )
        {
            y1 += stencilMatrix * x;
        }
    }
    {
        SCAI_REGION( "main.csrGEMV" )
        for ( IndexType iter = 0; iter < NITER; ++iter )
        {
            y2 += csrMatrix * x;
        }
    }

    // stencil and CSR matrix are same, so result vectors y1 and y2 must be same or close

    std::cout << "diff = " << y1.maxDiffNorm( y2 ) << std::endl;
}
