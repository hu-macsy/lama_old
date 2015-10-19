#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace scai::common;
using namespace scai::kregistry;

static void dummyRoutine()
{
}

BOOST_AUTO_TEST_CASE( SimpleTest )
{
    // This simple test registers a function in the kernel registry and uses it later

    KernelRegistry::set( dummyRoutine, "dummy", context::Host );
  
    KernelContextFunction<void(*)()> f( "dummy" );

    f[ context::Host ]();  // just call it

    // throw exception if called for CUDA, not registered

    BOOST_CHECK_THROW( 
    { 
        f[ context::CUDA ](); 

    }, KernelRegistryException );

    BOOST_CHECK_THROW( 
    { 
        KernelContextFunction<void(*)()> g( "dummy1" );  // wrong name

    }, KernelRegistryException );

    BOOST_CHECK_THROW( 
    { 
        KernelContextFunction<int(*)()> g( "dummy" );   // wrong signature

    }, KernelRegistryException );
}
