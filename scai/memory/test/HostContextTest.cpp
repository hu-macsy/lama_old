
#include <boost/test/unit_test.hpp>
 
#include <scai/memory/Context.hpp>

BOOST_AUTO_TEST_CASE( host_context )
{
    SCAI_LOG_DEF_LOGGER( logger, "Test" )

    SCAI_LOG_THREAD( "main" )

    using namespace scai::memory;

    BOOST_CHECK( Context::canCreate( context::Host ) );

    ContextPtr host = Context::create( context::Host, -1 );

    BOOST_CHECK( host.get() );

    SCAI_LOG_INFO( logger, "host context = " << *host )
}

