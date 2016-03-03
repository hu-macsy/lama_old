
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
