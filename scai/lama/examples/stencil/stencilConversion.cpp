/**
 * @file stencilConversion.cpp
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
 * @brief Benchmarking for conversion of StencilMatrix to CSRMatrix
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

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    common::Stencil2D<DefaultReal> stencil( 5 );

    const IndexType N = 10;

    common::Grid2D grid( N, N );

    StencilMatrix<DefaultReal> stencilMatrix( grid, stencil );

    CSRSparseMatrix<DefaultReal> csrMatrix1;

    {
        SCAI_REGION( "main.buildPoisson" )
        MatrixCreator::buildPoisson( csrMatrix1, 2, 5, N, N, N );
    }

    const CSRStorage<DefaultReal>& s = csrMatrix1.getLocalStorage();

    FileIO::write( s.getIA(), "ia.txt" );
    FileIO::write( s.getJA(), "ja.txt" );

    CSRSparseMatrix<DefaultReal> csrMatrix2;

    {
        SCAI_REGION( "main.convertStencil" )
        csrMatrix2 = stencilMatrix;
    }

    SCAI_ASSERT_EQUAL( csrMatrix1.getNumRows(), csrMatrix2.getNumRows(), "serious mismatch" )

    csrMatrix1.writeToFile( "poisson.mtx" );
    csrMatrix2.writeToFile( "stencil.mtx" );

    std::cout << "diff = " << csrMatrix1.maxDiffNorm( csrMatrix2 ) << std::endl;
}
