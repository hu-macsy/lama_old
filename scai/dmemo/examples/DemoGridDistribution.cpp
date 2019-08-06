/**
 * @file DemoGridDistribution.cpp
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
 * @brief Demo program for grid distribution
 * @author Thomas Brandes
 * @date 12.09.2018
 */

#include <scai/dmemo.hpp>

#include <scai/dmemo/GridDistribution.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Settings.hpp>

using namespace scai;
using namespace dmemo;

int main( int argc, const char* argv[] )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    common::Settings::setRank( comm->getNodeRank() );

    common::Settings::parseArgs( argc, argv );

    const IndexType N1 = 200;
    const IndexType N2 = 100;
    const IndexType N3 = 200;

    common::Grid3D domain( N1, N2, N3 );

    auto dist = gridDistribution( domain );   

    std::cout << dist->getCommunicator() << " : " << *dist << std::endl;
}
