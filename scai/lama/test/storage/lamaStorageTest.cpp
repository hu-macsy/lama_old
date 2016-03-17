/**
 * @file lamaStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Main program for test of LAMA storage classes
 * @author: Thomas Brandes
 * @date 16.03.2016
 **/

#ifndef BOOST_TEST_DYN_LINK
    #define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE lamaStorageTest

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <scai/hmemo/test/ContextFix.hpp>
#include <scai/common/Settings.hpp>

#include <cstdio>

/* Use global fixture for context and define the static variable.
 * Avoids expensive calls of reserve/release routines of the Context for each test.
 */

BOOST_GLOBAL_FIXTURE( ContextFix )

scai::hmemo::ContextPtr ContextFix::testContext;

/** The init function returns true if it can get the specified context. */

bool init_function()
{
    try
    {
        scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();
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
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    return boost::unit_test::unit_test_main( &init_function, argc, argv );
}
