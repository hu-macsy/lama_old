/**
 * @file DemoSettings.cpp
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
 * @brief Demo program for reading settings of individual processors
 * @author Thomas Brandes
 * @date 08.11.2018
 */

#include <scai/common/Settings.hpp>

#include <scai/dmemo.hpp>

using namespace scai;
using namespace dmemo;

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs(argc, argv);

    // get the default communicator (usually MPI if it has been enabled, or set by SCAI_COMMUNICATOR

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::string settingsFileName;

    if ( common::Settings::getEnvironment( settingsFileName, "SCAI_SETTINGS" ) )
    {
        int n = common::Settings::readSettingsFile( settingsFileName.c_str(), comm->getNodeName(), comm->getNodeRank() );

        std::cout << *comm << ", name = " << comm->getNodeName() << ", id = " << comm->getNodeId() 
                  << ", " << comm->getNodeRank() << " of " << comm->getNodeSize() 
                  << " : have got " << n << " settings from file " << settingsFileName << std::endl;
    }
    else
    {
        std::cerr << "ATTENTION: Environment variable SCAI_SETTINGS not set" << std::endl;
        return -1;
    }
    
    common::Settings::printEnvironment( std::cout );

    int domain;

    if ( common::Settings::getEnvironment( domain, "SCAI_DOMAIN" ) ) 
    {
    }
    else
    {
        std::cerr << "ATTENTION: environment variable SCAI_DOMAIN not set, comm = " << *comm << std::endl;
        return -1;
    }
    
    auto commDomain = comm->split( domain );

    std::cout << *comm << ", domain = " << domain << ", commDomain = " << *commDomain << std::endl;

    return 0;
}
