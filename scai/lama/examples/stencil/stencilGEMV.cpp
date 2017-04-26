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

    Stencil3D<double> stencil( 27 ); // 27-point stencil

    const IndexType N1 = 300;
    const IndexType N2 = 300;
    const IndexType N3 = 300;

    common::Grid3D grid( N1, N2, N3 );

    StencilMatrix<double> stencilMatrix( grid, stencil );

    std::cout << "stencilMatrix " << stencilMatrix << std::endl;

    stencilMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );

    CSRSparseMatrix<double> csrMatrix( stencilMatrix );

    csrMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );

    std::cout << "csrStencilMatrix " << csrMatrix << std::endl;

    DenseVector<double> x( stencilMatrix.getColDistributionPtr() );

    x = 1.0;

    DenseVector<double> y1 ( stencilMatrix.getNumRows(), 0.0 );
    DenseVector<double> y2 ( csrMatrix.getNumRows(), 0.0 );

    for ( IndexType iter = 0; iter < 10; ++iter )
    {
        y1 += stencilMatrix * x;
        y2 += csrMatrix * x;
    }

    std::cout << "diff = " << y1.maxDiffNorm( y2 ) << std::endl;
}
