/**
 * @file DemoDistribution.cpp
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
 * @brief Demo program for distribution + communication
 * @author Thomas Brandes
 * @date 10.02.2016
 */

#include <scai/dmemo.hpp>

using namespace scai;
using namespace dmemo;

int main()
{
    SCAI_LOG_THREAD( "Main" )
    // get the default communicator (usually MPI if it has been enabled, or set by SCAI_COMMUNICATOR
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    IndexType size = 71;
    float weight = 1.0;
    DistributionPtr dist ( Distribution::getDistributionPtr( "CYCLIC", comm, size, weight ) );
    // Note: distribution pointers are always const pointers, so distributions can never be changed
    std::cout << *comm << ", dist = " << *dist << std::endl;
}
