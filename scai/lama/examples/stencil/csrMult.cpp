/**
 * @file lama/examples/stencil/csrMult.cpp
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
 * @brief Benchmark to compare stencil GEMV with CSR GEMV.
 * @author Thomas Brandes
 * @date 26.04.2017
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

static const IndexType NITER = 10;

int main( int argc, const char* argv[] )
{
    typedef DefaultReal ValueType;

    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    if ( argc < 2 )
    {
        std::cout << "call: " << argv[0] << " matrix" << std::endl;
        return -1;
    }

    // read in a CSR matrix from file, name specified by command line argument

    auto csrMatrix = read<CSRSparseMatrix<ValueType>>( argv[1] );

    auto x = fillDenseVector<ValueType>( csrMatrix.getColDistributionPtr(), 1 );
    auto y = fillDenseVector<ValueType>( csrMatrix.getRowDistributionPtr(), 0 );

    for ( IndexType i = 0; i < NITER; ++i )
    {
        y += csrMatrix * x;
    }

    std::cout << "y norm = " << y.l2Norm() << std::endl;
}
