/**
 * @file stencilGEMV.cpp
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

#include <scai/tracing.hpp>

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

    // Distribute grid onto default processor array, can be set by --SCAI_NP=2x3x2

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    ContextPtr ctx = Context::getContextPtr();

    dmemo::DistributionPtr gridDistribution( new GridDistribution( grid, comm ) );

    std::cout << *comm << ": distribution = " << *gridDistribution << std::endl;

    lama::SyncKind syncKind = lama::SyncKind::SYNCHRONOUS;

    // by default do matrix-vector operations synchronously, but can be set via environment

    bool isSet;

    if ( common::Settings::getEnvironment( isSet, "SCAI_ASYNCHRONOUS" ) )
    {
        if ( isSet )
        {
            syncKind = scai::lama::SyncKind::ASYNC_LOCAL;
        }
    }

    // The stencil matrix just needs the grid distribution and the stencil

    StencilMatrix<ValueType> stencilMatrix( gridDistribution, stencil );
    stencilMatrix.setContextPtr( ctx );
    stencilMatrix.setCommunicationKind( syncKind );

    std::cout << *comm << ": stencilMatrix " << stencilMatrix << std::endl;

    // Conversion stencil matrix -> CSR matrix is full supported

    CSRSparseMatrix<ValueType> csrMatrix( stencilMatrix );

    csrMatrix.setContextPtr( ctx );
    csrMatrix.setCommunicationKind( syncKind );

    std::cout << "csrStencilMatrix " << csrMatrix << std::endl;

    auto x = fill<DenseVector<ValueType>>( stencilMatrix.getColDistributionPtr(), 1 );

    auto y1 = fill<DenseVector<ValueType>>( stencilMatrix.getRowDistributionPtr(), 0 );
    auto y2 = fill<DenseVector<ValueType>>( csrMatrix.getRowDistributionPtr(), 0 );

    double timeStencil;
    double timeCSR;

    {
        SCAI_REGION( "bench.stencilGEMV" )

        double time = common::Walltime::get();

        for ( IndexType iter = 0; iter < NITER; ++iter )
        {
            y1 += stencilMatrix * x;
        }

        timeStencil = common::Walltime::get() - time;

        // compute time for one stencil in ms

        timeStencil = ( timeStencil / NITER) * 1000.0;
    }
    {
        SCAI_REGION( "bench.csrGEMV" )

        double time = common::Walltime::get();

        for ( IndexType iter = 0; iter < NITER; ++iter )
        {
            y2 += csrMatrix * x;
        }

        timeCSR = common::Walltime::get() - time;

        // compute time for one stencil in ms

        timeCSR = ( timeCSR / NITER) * 1000.0;
    }

    // stencil and CSR matrix are same, so result vectors y1 and y2 must be same or close

    std::cout << "diff = " << y1.maxDiffNorm( y2 ) << std::endl;

    std::cout << "Time (ms): stencil = " << timeStencil << ", csr = " << timeCSR << std::endl;
}

int main( int argc, const char* argv[] )
{
    typedef DefaultReal ValueType;

    // relevant SCAI arguments: 
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    if ( argc < 3 )
    {
        std::cout << "call: " << argv[0] << " [options] <ndims> <stencilType>" << std::endl;
        std::cout << "  possible options:" << std::endl;
        std::cout << "    --SCAI_ASYNCHRONOUS=0|1" << std::endl;
        std::cout << "    --SCAI_CONTEXt=Host|CUDA" << std::endl;
        return -1;
    }

    ContextPtr ctx = Context::getContextPtr();

    // on GPUs we can only run smaller problems

    IndexType nDims   = atoi( argv[1] );
    IndexType nPoints = atoi( argv[2] );

    std::cout << "Bench " << argv[0] << " nDims = " << nDims << ", nPoints = " << nPoints << std::endl;

    if ( nDims == 4 )
    { 
        common::Stencil4D<ValueType> stencil( nPoints );

        // Define a grid with same number of dimensions as stencil

        const IndexType N1 = 50;
        const IndexType N2 = 80;
        const IndexType N3 = 100;
        const IndexType N4 = 100;
    
        common::Grid4D grid( N1, N2, N3, N4 );

        bench( grid, stencil );
    }
    else if ( nDims == 3 )
    { 
        common::Stencil3D<ValueType> stencil( nPoints );

        // Define a grid with same number of dimensions as stencil

        const IndexType N1 = 200;
        const IndexType N2 = 400;
        const IndexType N3 = 500;
    
        common::Grid3D grid( N1, N2, N3 );

        bench( grid, stencil );
    }
    else if ( nDims == 2 )
    { 
        common::Stencil2D<ValueType> stencil( nPoints );

        // Define a grid with same number of dimensions as stencil

        const IndexType N1 =  4000;
        const IndexType N2 = 10000;
    
        common::Grid2D grid( N1, N2 );

        // grid.setBorderType( 0, common::Grid::BORDER_REFLECTING, common::Grid::BORDER_REFLECTING );
        // grid.setBorderType( 0, common::Grid::BORDER_PERIODIC, common::Grid::BORDER_PERIODIC );

        bench( grid, stencil );
    }
    else if ( nDims == 1 )
    { 
        common::Stencil1D<ValueType> stencil( nPoints );

        // Define a grid with same number of dimensions as stencil

        const IndexType N1 = 40 * 1000 * 1000;
    
        common::Grid1D grid( N1 );

        // grid.setBorderType( 0, common::Grid::BORDER_PERIODIC );

        bench( grid, stencil );
    }
    else
    {
        COMMON_THROWEXCEPTION( "ndims = " << nDims << " not supported yet" );
    }
}
