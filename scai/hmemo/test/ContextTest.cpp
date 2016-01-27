#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>

#include "MockContext.hpp"

using namespace scai::common::context;

BOOST_AUTO_TEST_CASE( ContextTest )
{
    SCAI_LOG_DEF_LOGGER( logger, "Test" )

    SCAI_LOG_THREAD( "main" )

    using namespace scai::hmemo;

    ContextPtr userContext  = Context::getContextPtr( UserContext, 1 );
    ContextPtr userContext2 = Context::getContextPtr( UserContext, 2 );
    ContextPtr hostContext  = Context::getContextPtr( Host );

    SCAI_LOG_INFO( logger, "userContext = " << *userContext );

    HArray<double> X( 10, 5.0 );

    {
        WriteAccess<double> write( X, userContext );  
    }

    // read @ userContext2: valid data is transfered from userContext to here

    ReadAccess<double> read( X, userContext2 );

    const double* vals = read.get();

    for ( int i = 0; i < 10; ++i )
    {
        SCAI_ASSERT_EQUAL( vals[i], 5.0, "check" )
    }


    Y.clear();

    Y.purge();

    HArray<float> v ( 4, 1.0f );  

    {
        ReadAccess<float> read( v, userContext );
        WriteAccess<float> write( v, userContext );
    }

    try
    {
        ReadAccess<float> read( v, userContext );
        WriteAccess<float> write( v, userContext2 );
        COMMON_THROWEXCEPTION( "read and write access at same time not possible" )
    }
    catch ( scai::common::Exception& ex )
    {
        //std::cout << "Exception caught: " << ex.what() << std::endl;
    }
}

