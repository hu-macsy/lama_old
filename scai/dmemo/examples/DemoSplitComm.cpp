/**
 * @file DemoSplitComm.cpp
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
 * @brief Demo program for splitting communicator
 * @author Thomas Brandes
 * @date 08.11.2018
 */

#include <scai/dmemo.hpp>

using namespace scai;
using namespace dmemo;

int main()
{
    // get the default communicator (usually MPI if it has been enabled, or set by SCAI_COMMUNICATOR

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId color = comm->getRank() % 2 == 0 ? 0 : 1;
    PartitionId key   = comm->getSize() - comm->getRank();  // reverse order
 
    CommunicatorPtr comm1 = comm->split( color, key );

    std::cout << *comm << ": color = " << color << ", new communicator = " << *comm1 << std::endl;
}
