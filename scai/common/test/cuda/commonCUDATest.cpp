
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
