/**
 * @file dmemo/test/dmemoTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief ToDo: Missing description in ./dmemo/test/dmemoTest.cpp
 * @author Thomas Brandes
 * @date 10.02.2016
 */

#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE dmemoTest

// indicate that default main of Boost is not used here

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <scai/common/Settings.hpp>
#include <scai/dmemo.hpp>

#include <iostream>

bool init_function()
{
    try
    {
        scai::dmemo::CommunicatorPtr comm = scai::dmemo::Communicator::getCommunicatorPtr();
        return true;
    }
    catch ( scai::common::Exception& ex )
    {
        return false;
    }
}

int main( int argc, char* argv[] )
{
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    return boost::unit_test::unit_test_main( &init_function, argc, argv );
}
