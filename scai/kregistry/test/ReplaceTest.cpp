#include <boost/test/unit_test.hpp>

#include <scai/kregistry/KernelContextFunction.hpp>

using namespace scai;
using namespace scai::kregistry;
using namespace scai::common;

template<typename ValueType>
static ValueType add1( const ValueType x )
{
    return x + 1.0;
}

template<typename ValueType>
static ValueType minus1( const ValueType x ) 
{
    return x - 1.0;
}

/** Trait to handle function ValueType ( fn ) ( ValueType ) in KernelRegistry. */

template<typename ValueType>
struct UnaryOpTrait
{
    typedef ValueType ( *FuncType ) ( ValueType );
    static const char* getId() { return "UnaryOp"; }  
};

BOOST_AUTO_TEST_CASE( ReplaceTest )
{
    // register unary operator for double

    KernelRegistry::set<UnaryOpTrait<double> >( add1<double>, context::Host );
    KernelRegistry::set<UnaryOpTrait<double> >( minus1<double>, context::Host );  // does not overwrite add1
  
    // register unary operator for float

    KernelRegistry::set<UnaryOpTrait<float> >( add1<float>, context::Host );
    KernelRegistry::set<UnaryOpTrait<float> >( minus1<float>, context::Host, true );  // overrides add1

    KernelTraitContextFunction<UnaryOpTrait<float> > opFloat;
    KernelTraitContextFunction<UnaryOpTrait<double> > opDouble;

    double xd = opDouble[ context::Host ]( 1.0 );
    float  xf = opFloat [ context::Host ]( 1.0f );

    BOOST_CHECK_EQUAL( 2.0, xd );
    BOOST_CHECK_EQUAL( 0.0, xf );
}

