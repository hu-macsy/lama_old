/**
 * @file sectionExample.cpp
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
 * @brief Example program to show how to work on grid sections
 * @author Thomas Brandes
 * @date 16.05.2017
 */

#include <scai/lama/GridVector.hpp>
#include <scai/lama/GridReadAccess.hpp>
#include <scai/lama/GridWriteAccess.hpp>
#include <scai/lama/GridSection.hpp>

#include <scai/lama/io/ImageIO.hpp>
#include <scai/lama/io/MatlabIO.hpp>

#include <scai/common/Settings.hpp>

#include <scai/lama.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

using namespace common;

typedef double real;
typedef ComplexDouble complex;

SCAI_LOG_DEF_LOGGER( logger, "main" )

int main( int argc, const char* argv[] )
{
    // relevant SCAI arguments:
    //   SCAI_CONTEXT = ...    set default context
    //   SCAI_DEVICE  = ...    set default device

    common::Settings::parseArgs( argc, argv );

    if ( argc < 1 )
    {
        std::cout << "Wrong call, please use : " << argv[0] << " <inputFileName> <outputFileName>" << std::endl;
    }

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    IndexType np = comm->getSize();

    DistributionPtr dist( new GridDistribution( Grid1D( 100 ), comm, Grid1D( np ) ) );

    GridVector<double> x( dist, 1.0 );

    x( Range( 25, 45, 3 ) ) = 2.0;

    std::cout << "x = " << x << std::endl;

    x.writeToFile( "x.txt" );
}
