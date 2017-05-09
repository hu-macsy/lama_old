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

static const IndexType NITER = 10;

template<typename ValueType> 
void bench( const common::Grid& grid, const common::Stencil<ValueType>& stencil )
{
    std::cout << "Benchmark: grid = " << grid << ", stencil = " << stencil << std::endl;

    // Distibute grid onto default processor array, can be set by --SCAI_NP=2x3x2

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    ContextPtr ctx = Context::getContextPtr();

    dmemo::DistributionPtr gridDistribution( new GridDistribution( grid, comm ) );

    // The stencil matrix just needs the grid distribution and the stencil

    StencilMatrix<ValueType> stencilMatrix( gridDistribution, stencil );
    stencilMatrix.setContextPtr( ctx );

    std::cout << "stencilMatrix " << stencilMatrix << std::endl;

    stencilMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );

    // Conversion stencil matrix -> CSR matrix is full supported

    CSRSparseMatrix<ValueType> csrMatrix( stencilMatrix );

    csrMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );
    csrMatrix.setContextPtr( ctx );

    std::cout << "csrStencilMatrix " << csrMatrix << std::endl;

    DenseVector<ValueType> x( stencilMatrix.getColDistributionPtr() );

    x = 1.0;

    DenseVector<ValueType> y1 ( stencilMatrix.getRowDistributionPtr(), 0.0 );
    DenseVector<ValueType> y2 ( csrMatrix.getRowDistributionPtr(), 0.0 );

    {
        SCAI_REGION( "bench.stencilGEMV" )

        for ( IndexType iter = 0; iter < NITER; ++iter )
        {
            y1 += stencilMatrix * x;
        }
    }
    {
        SCAI_REGION( "bench.csrGEMV" )
        for ( IndexType iter = 0; iter < NITER; ++iter )
        {
            y2 += csrMatrix * x;
        }
    }

    // stencil and CSR matrix are same, so result vectors y1 and y2 must be same or close

    std::cout << "diff = " << y1.maxDiffNorm( y2 ) << std::endl;
}

int main( int argc, const char* argv[] )
{
    typedef float ValueType;

    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    if ( argc < 3 )
    {
        std::cout << "call: " << argv[0] << " <ndims> <stencilType>" << std::endl;
        return -1;
    }

    ContextPtr ctx = Context::getContextPtr();

    // on GPUs we can only run smaller problems

    bool isGPU = ctx->getType() == Context::CUDA;

    IndexType nDims   = atoi( argv[1] );
    IndexType nPoints = atoi( argv[2] );

    std::cout << "Bench " << argv[0] << " nDims = " << nDims << ", nPoints = " << nPoints << std::endl;

    if ( nDims == 3 )
    { 
        common::Stencil3D<ValueType> stencil( nPoints );

        // Define a grid with same number of dimensions as stencil

        const IndexType N1 = isGPU ? 200 : 400;
        const IndexType N2 = 400;
        const IndexType N3 = 500;
    
        common::Grid3D grid( N1, N2, N3 );

        bench( grid, stencil );
    }
    else if ( nDims == 2 )
    { 
        common::Stencil2D<ValueType> stencil( nPoints );

        // Define a grid with same number of dimensions as stencil

        const IndexType N1 = isGPU ? 4000 : 8000;
        const IndexType N2 = 10000;
    
        common::Grid2D grid( N1, N2 );

        bench( grid, stencil );
    }
    else if ( nDims == 1 )
    { 
        common::Stencil1D<ValueType> stencil( nPoints );

        // Define a grid with same number of dimensions as stencil

        const IndexType N1 = isGPU ? 40000000 : 80000000;
    
        common::Grid1D grid( N1 );

        bench( grid, stencil );
    }
    else
    {
        COMMON_THROWEXCEPTION( "ndims = " << nDims << " not supported yet" );
    }
}
