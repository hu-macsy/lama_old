
#ifndef BOOST_TEST_DYN_LINK
    #define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE hmemoTest

// indicate that default main of Boost is not used here

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <scai/hmemo.hpp>
#include <scai/hmemo/test/MockContext.hpp>

#include <scai/common/Settings.hpp>
#include <scai/logging.hpp>

#include <iostream>

// one time registration of MockContext required via guard variable

template scai::hmemo::Context::Register<MockContext>::RegisterGuard 
    scai::hmemo::Context::Register<MockContext>::registerGuard;

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
    SCAI_LOG_THREAD( "main" )

    scai::common::Settings::parseArgs( argc, const_cast<const char**>( argv ) );

    return boost::unit_test::unit_test_main( &init_function, argc, argv );
}
