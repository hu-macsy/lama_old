/**
 * @file BLASKernelTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief ToDo: Missing description in ./BLASKernelTest.cpp
 * @author Eric Schricker
 * @date 20.11.2015
 */

#ifndef BOOST_TEST_DYN_LINK
    #define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE blaskernelTest

// indicate that default main of Boost is not used here

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <scai/hmemo.hpp>

#include <scai/common/Settings.hpp>
#include <scai/logging.hpp>

#include <iostream>

#include <scai/blaskernel/test/TestMacros.hpp>

BOOST_GLOBAL_FIXTURE( ContextFix );

/** Static variables of ContextFix are defined here */

scai::hmemo::ContextPtr ContextFix::testContext;

/** The init function returns true if it can get the specified context. */

bool init_function()
{
    try
    {
        scai::hmemo::ContextPtr testContext = scai::hmemo::Context::getContextPtr();
        return true;
    }
    catch ( scai::common::Exception& ex )
    {
        std::cerr << "Could not get context for test: " << ex.what() << std::endl;
        return false;
    }
} 

int main( int argc, char* argv[] )
{
    SCAI_LOG_THREAD( "main" )

    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    int rc = boost::unit_test::unit_test_main( &init_function, argc, argv );
 
    return rc;
}
