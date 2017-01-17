/**
 * @file test/cuda/commonCUDATest.cpp
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
 * @brief ToDo: Missing description in ./test/cuda/commonCUDATest.cpp
 * @author Thomas Brandes
 * @date 08.03.2016
 */

#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif

// indicate that default main of Boost is not used here

#define BOOST_TEST_NO_MAIN

#define BOOST_TEST_MODULE CommonCUDATest

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <scai/common/Settings.hpp>

#include <iostream>

/** The init function just returns true */

bool init_function()
{
    return true;
}

int main( int argc, char* argv[] )
{
    // parse command line argument, SCAI_DEVICE may be set
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );
    return boost::unit_test::unit_test_main( &init_function, argc, argv );
}
