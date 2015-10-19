#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace scai::kregistry;

static double add1( const double x )
{
    return x + 1.0;
}

static double minus1( const double x )
{
    return x - 1.0;
}

/** Trait to handle function double ( fn ) ( double ) in KernelRegistry. */

struct UnaryAddTrait
{
    typedef double ( *FuncType ) ( double );
    static const char* getId() { return "add"; }
};

struct UnaryMinusTrait
{
    typedef double ( *FuncType ) ( double );
    static const char* getId() { return "minus"; }
};

BOOST_AUTO_TEST_CASE( ContextTest )
{
    // register context::Host

    KernelRegistry::set<UnaryAddTrait>( add1, context::Host );
    KernelRegistry::set<UnaryMinusTrait>( minus1, context::Host );

    // register context::CUDA

    KernelRegistry::set<UnaryAddTrait>( add1, context::CUDA );
    KernelRegistry::set<UnaryMinusTrait>( minus1, context::CUDA );

    // register context::MIC, only add1

    KernelRegistry::set<UnaryAddTrait>( add1, context::MIC );

    // register context::UserContext, only minus1

    KernelRegistry::set<UnaryMinusTrait>( minus1, context::UserContext );

    KernelRegistry::printAll();

    KernelTraitContextFunction<UnaryAddTrait> add;
    KernelTraitContextFunction<UnaryMinusTrait> minus;

    // add can be alled at context::MIC

    BOOST_CHECK_EQUAL( context::MIC, add.validContext( context::MIC ) );

    // minus must be called at context::Host

    BOOST_CHECK_EQUAL( context::Host, minus.validContext( context::MIC ) );

    // add, minus can be called together @ CUDA

    BOOST_CHECK_EQUAL( context::CUDA, add.validContext( minus, context::CUDA ) );
    BOOST_CHECK_EQUAL( context::Host, add.validContext( minus, context::MIC ) );
    BOOST_CHECK_EQUAL( context::Host, add.validContext( minus, context::UserContext ) );
    BOOST_CHECK_EQUAL( context::Host, add.validContext( minus, context::Host ) );
}

