/**
 * @file lama/examples/tutorial/tasking2.cpp
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
 * @brief Example of using task parallelism
 * @author Thomas Brandes
 * @date 08.11.2018
 */
#include <scai/lama.hpp>

#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

int main( int argc, const char* argv[] )
{
    // start with communicator for all processors (WORLD)

    CommunicatorPtr commWorld = Communicator::getCommunicatorPtr();

    SCAI_REGION( "main.tasking" )

    // parse command line arguments, here to use SCAI_NP=NP1xNP2 on command line

    {
        SCAI_REGION( "main.parseArgs." )
        common::Settings::parseArgs( argc, argv );
    }

    IndexType NP  = 1000;    // (default) number of problems
    IndexType NV1 = 1000;    // (default) size dim 1 of each problem
    IndexType NV2 = 1000;    // (default) size dim 2 of each problem

    if ( argc >= 2 )
    {
        NP = atoi( argv[1] );
    }

    if ( argc >= 3 )
    {
        NV1 = atoi( argv[2] );
    }

    if ( argc >= 3 )
    {
        NV2 = atoi( argv[3] );
    }

    // build space = #problems x #size1 x #size2

    common::Grid3D space( NP, NV1, NV2 );

    // (block) distribute the problem space onto the available processors

    auto dist = std::make_shared<dmemo::GridDistribution>( space );

    std::cout << *commWorld << ": distributed grid = " << *dist << std::endl;

    // find out the position of this processor in the processor grid

    const common::Grid& procGrid = dist->getProcGrid();

    IndexType procGridRank[3];
    procGrid.gridPos( procGridRank, commWorld->getRank() );

    CommunicatorPtr intraProblem;  // communicator for set of processors that solve one problem
    CommunicatorPtr interProblem;  // communicator between processor sets 

    {
        SCAI_REGION( "main.split" )

        // all processors with same first dimension in processor array build one problem set

        intraProblem = commWorld->split( procGridRank[0] );

        // this communicator is used for reducing the solutions of problems

        interProblem = commWorld->split( intraProblem->getRank() );
    }

    std::cout << *commWorld << ": is ( " << procGridRank[0] << ", " << procGridRank[1] 
                       << ", " << procGridRank[2] << " ) in " << procGrid
                       << ", intraProblem = " << *intraProblem << ", interProblem = " << *interProblem << std::endl;
    double sum = 0;

    // set up grid and distibution for one problem

    common::Grid2D problemSpace( NV1, NV2 );
    common::Grid2D procArray( procGrid.size(1), procGrid.size(2 ) );

    auto distProblem = std::make_shared<dmemo::GridDistribution>( problemSpace, intraProblem, procArray );

    // now loop over my assigned problems that is first dimension of distributed grid

    for ( IndexType i = dist->localLB()[0]; i < dist->localUB()[0]; ++i ) 
    {
        SCAI_REGION( "main.problem" )

        // Note: each problem has a different data set

        auto v = denseVectorLinear<double>( distProblem, double(i), double(i) / ( NV1 * NV2 ) );
        sum += v.sum();  
    }
 
    double allSum = interProblem->sum( sum );

    std::cout << *commWorld << ": local sum = " << sum << ", all = " << allSum << std::endl;
}
