/**
 * @file lama/examples/tutorial/tasking.cpp
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
#include <scai/dmemo/CommunicatorStack.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

typedef DefaultReal ValueType;     // double if enabled, otherwise float

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

    IndexType NP = 10;       // (default) number of problems
    IndexType NV = 1000000;    // (default) size of each problem

    if ( argc >= 2 )
    {
        NP = atoi( argv[1] );
    }

    if ( argc >= 3 )
    {
        NV = atoi( argv[2] );
    }

    // build space = #problems x #problemSize

    common::Grid2D space( NP, NV );

    // (block) distribute the problem space onto the available processors

    auto dist = std::make_shared<dmemo::GridDistribution>( space );

    std::cout << *commWorld << ": distributed grid = " << *dist << std::endl;

    // find out the position of this processor in the processor grid

    const common::Grid& procGrid = dist->getProcGrid();

    IndexType procGridRank[2];
    procGrid.gridPos( procGridRank, commWorld->getRank() );

    CommunicatorPtr intraProblem;  // communicator for set of processors that solve one problem
    CommunicatorPtr interProblem;  // communicator between processor sets 

    // this communicator is used for distributions of data in one problem

    {
        SCAI_REGION( "main.split" )

        intraProblem = commWorld->split( procGridRank[0] );

        // this communicator is used for reducing the solutions of problems

        interProblem = commWorld->split( procGridRank[1] );

    }

    std::cout << *commWorld << ": is ( " << procGridRank[0] << ", " << procGridRank[1] << " ) in " << procGrid
                       << ", intraProblem = " << *intraProblem << ", interProblem = " << *interProblem << std::endl;
    ValueType sum = 0;

    auto result = denseVectorFill<ValueType>( blockDistribution( 100, intraProblem ), 0 );

    // now loop over my assigned problems that is first dimension of distributed grid

    for ( IndexType i = dist->localLB()[0]; i < dist->localUB()[0]; ++i ) 
    {
        SCAI_REGION( "main.problem" )

        SCAI_DMEMO_TASK( intraProblem )

        auto dist = blockDistribution( NV );  // uses now default communicator intraProblem
        auto v = denseVectorLinear<ValueType>( dist, ValueType(i), ValueType(1) / NV );
        std::ostringstream out;
        out << "out_" << i << "_.mtx";
        v.writeToFile( out.str() );
        auto v1 = read<SparseVector<ValueType>>( out.str() );
        v1.redistribute( dist );
        sum += v1.sum();  
        result += 1.0;
    }
 
    ValueType allSum = interProblem->sum( sum );
    interProblem->sumArray( result.getLocalValues() );

    std::cout << *commWorld << ": local sum = " << sum << ", all = " << allSum << std::endl;
}
