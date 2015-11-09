#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace scai::common;
using namespace scai::kregistry;

static int dummyRoutine()
{
    return 15;
}

/** Trait structure for registration of int routine with name "MyDummy" */

struct TraitDummyRoutine
{
    typedef int ( *FuncType ) ();    // signature of the function
    static const char* getId() { return "MyDummy"; }
};

BOOST_AUTO_TEST_CASE( TraitTest )
{
    // Same as simple test but uses a Trait for registration
    // The trait avoids misspelling of the routine name and the signature

    KernelRegistry::set<TraitDummyRoutine>( dummyRoutine, context::CUDA, KernelRegistry::KERNEL_ADD );

    KernelTraitContextFunction<TraitDummyRoutine> f;

    int x = f[ context::CUDA ]();  // just call it

    BOOST_CHECK_EQUAL( 15, x );

    // throw exception if called for MIC, not registered

    BOOST_CHECK_THROW( 
    { 
        x = f[ context::MIC ](); 

    }, KernelRegistryException );

    // misspelling of name or signature is no more possible here, so no further test for failure
}

