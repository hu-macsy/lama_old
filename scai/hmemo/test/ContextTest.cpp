#include <boost/test/unit_test.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/LAMAArray.hpp>
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

    LAMAArray<double> X( 10, 5.0 );

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

    // Now make some checks

    std::cout << "X @ " << *userContext << ", valid = " << X.isValid( userContext )
              << ", capacity = " << X.capacity( userContext ) << std::endl;

    std::cout << "X @ " << *userContext2 << ", valid = " << X.isValid( userContext2 )
              << ", capacity = " << X.capacity( userContext2 ) << std::endl;

    std::cout << "X @ " << *hostContext << ", valid = " << X.isValid( hostContext )
              << ", capacity = " << X.capacity( hostContext ) << std::endl;

    LAMAArray<double> Y( X );

    // valid should be the same for Y, capacity should be 0 if not valid

    std::cout << "Y @ " << *userContext << ", valid = " << Y.isValid( userContext )
              << ", capacity = " << Y.capacity( userContext ) << std::endl;

    std::cout << "Y @ " << *userContext2 << ", valid = " << Y.isValid( userContext2 )
              << ", capacity = " << Y.capacity( userContext2 ) << std::endl;

    std::cout << "Y @ " << *hostContext << ", valid = " << Y.isValid( hostContext )
              << ", capacity = " << Y.capacity( hostContext ) << std::endl;

    Y.clear();

    std::cout << "Y cleared now" << std::endl;

    // valid should be the same for Y, capacity should be 0 if not valid

    std::cout << "Y @ " << *userContext << ", valid = " << Y.isValid( userContext )
              << ", capacity = " << Y.capacity( userContext ) << std::endl;

    std::cout << "Y @ " << *userContext2 << ", valid = " << Y.isValid( userContext2 )
              << ", capacity = " << Y.capacity( userContext2 ) << std::endl;

    std::cout << "Y @ " << *hostContext << ", valid = " << Y.isValid( hostContext )
              << ", capacity = " << Y.capacity( hostContext ) << std::endl;

    Y.purge();

    std::cout << "Y purged now" << std::endl;

    // valid should be the same for Y, capacity should be 0 if not valid

    std::cout << "Y @ " << *userContext << ", valid = " << Y.isValid( userContext )
              << ", capacity = " << Y.capacity( userContext ) << std::endl;

    std::cout << "Y @ " << *userContext2 << ", valid = " << Y.isValid( userContext2 )
              << ", capacity = " << Y.capacity( userContext2 ) << std::endl;

    std::cout << "Y @ " << *hostContext << ", valid = " << Y.isValid( hostContext )
              << ", capacity = " << Y.capacity( hostContext ) << std::endl;

    int values[] = { 1, 2, 3, 4 };

    LAMAArray<float> v ( 4, values );   // implicit type conversion allowed

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
        std::cout << "Exception caught: " << ex.what() << std::endl;
    }
}

