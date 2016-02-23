
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MemoryTest

// indicate that default main of Boost is not used here

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <scai/common/Settings.hpp>
#include <scai/hmemo.hpp>

#include <iostream>

bool init_function()
{
    // scai::hmemo::ContextPtr ctx = scai::hmemo::Context::getContextPtr();
    // std::cout << "Running test for context " << *ctx << std::endl;

    return true;
} 

int main( int argc, char* argv[] )
{
    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    return boost::unit_test::unit_test_main( &init_function, argc, argv );
}

