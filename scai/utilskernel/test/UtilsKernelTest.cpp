
#ifndef BOOST_TEST_DYN_LINK
    #define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE utilskernelTest

// indicate that default main of Boost is not used here

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include <scai/hmemo.hpp>

#include <scai/common/Settings.hpp>
#include <scai/logging.hpp>

#include <iostream>

scai::hmemo::ContextPtr testContext;

/** The init function returns true if it can get the specified context. */

bool init_function()
{
    try
    {
        testContext = scai::hmemo::Context::getContextPtr();
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
 
    testContext.reset();   // frees the context

    return rc;
}
