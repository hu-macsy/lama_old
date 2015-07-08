
#include <boost/test/unit_test.hpp>
 
#include <memory/Context.hpp>

BOOST_AUTO_TEST_CASE( host_context )
{
    LAMA_LOG_DEF_LOGGER( logger, "Test" )

    LAMA_LOG_THREAD( "main" )

    using namespace memory;

    BOOST_CHECK( Context::canCreate( context::Host ) );

    ContextPtr host = Context::create( context::Host, -1 );

    BOOST_CHECK( host.get() );

    LAMA_LOG_INFO( logger, "host context = " << *host )
}

