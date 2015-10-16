#include <boost/test/unit_test.hpp>

#include <scai/kernel/KernelContextFunction.hpp>

using namespace scai::common;
using namespace scai::interface;

static int dummyRoutine()
{
    return 15;
}

struct TraitDummyRoutine
{
    typedef int ( *FuncType ) ();
    static const char* getId() { return "MyDummy"; }
};

BOOST_AUTO_TEST_CASE( TraitTest )
{
    // This simple test registers a function in the interface and uses it later

    KernelInterface::set<TraitDummyRoutine>( dummyRoutine, context::Host );

    KernelTraitContextFunction<TraitDummyRoutine> f;

    int x = f[ context::Host ]();  // just call it

    BOOST_CHECK_EQUAL( 15, x );

    // throw exception if called for CUDA, not registered

    BOOST_CHECK_THROW( 
    { 
        x = f[ context::CUDA ](); 

    }, Exception );

    // misspelling of name or signature is no more possible here, so everything is fine
}

