/**
 * @file DemoGenBlocklDistribution.cpp
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
 * @brief Demo program for geneneral block distribution
 * @author Thomas Brandes
 * @date 12.09.2018
 */

#include <scai/dmemo.hpp>

#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Settings.hpp>

using namespace scai;
using namespace dmemo;

int main()
{
    const IndexType N = 1000;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    float weight = 0.5;

    common::Settings::setRank( comm->getRank() );
    common::Settings::getEnvironment( weight, "SCAI_WEIGHT" );

    auto dist = genBlockDistributionByWeight( N, weight, comm );

    std::cout << dist->getCommunicator() << " : " << *dist << " by weight " << weight << std::endl;
}
