/**
 * @file tasking/test/cuda/taskingCUDATest.cpp
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
 * @brief ToDo: Missing description in ./tasking/test/cuda/taskingCUDATest.cpp
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
#include <scai/testsupport/commonTestMain.hpp>

int main( int argc, char* argv[] )
{
    // parse command line argument, SCAI_DEVICE may be set
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );
    return scai::testsupport::commonTestMain( argc, argv );
}
