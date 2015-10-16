#include <boost/test/unit_test.hpp>

#include <scai/kernel/KernelContextFunction.hpp>

using namespace scai::common;
using namespace scai::interface;

static void dummyRoutine()
{
}

BOOST_AUTO_TEST_CASE( SimpleTest )
{
    // This simple test registers a function in the interface and uses it later

    KernelInterface::set( dummyRoutine, "dummy", scai::common::context::Host );
  
    KernelContextFunction<void(*)()> f( "dummy" );

    f[ context::Host ]();  // just call it

    // throw exception if called for CUDA, not registered

    BOOST_CHECK_THROW( 
    { 
        f[ context::CUDA ](); 

    }, Exception );

    BOOST_CHECK_THROW( 
    { 
        KernelContextFunction<void(*)()> g( "dummy1" );  // wrong name

    }, Exception );

    BOOST_CHECK_THROW( 
    { 
        KernelContextFunction<int(*)()> g( "dummy" );   // wrong signature

    }, Exception );
}

