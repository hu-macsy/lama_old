#include <boost/test/unit_test.hpp>


#include "MockContext.hpp"

#include <scai/hmemo.hpp>
#include <scai/common/ContextType.hpp>

using namespace scai;
using namespace scai::hmemo;
using namespace scai::common::context;

BOOST_AUTO_TEST_CASE( ReleaseTest )
{
    ContextPtr userContext  = Context::getContextPtr( UserContext, 1 );
    ContextPtr hostContext  = Context::getContextPtr( Host );

    HArray<double> X( 10, 5.0 );

    // write access on UserContext
    {
        WriteAccess<double> write( X, userContext );  
    }
    // read access on Host, okay as write is released
    {
        ReadAccess<double> read( X, hostContext );  
    }

    common::function<void()> delay;

    // write access on UserContext, but delay the release
    {
        WriteAccess<double> write( X, userContext );  
        delay = write.releaseDelayed();

        // write access can no more be used

        try
        {
            double* ptr = write.get();  
            ptr[0] = 0;
        }
        catch ( scai::common::Exception& ex )
        {
            std::cout << "Exception caught: " << ex.what() << std::endl;
        }
    }

    // read access on Host, not possible as release on write is not done
    try
    {
        ReadAccess<double> read( X, hostContext );  
    }
    catch ( scai::common::Exception& ex )
    {
        std::cout << "Exception caught: " << ex.what() << std::endl;
    }

    delay();  // now it is okay

    {   
        ReadAccess<double> read( X, hostContext );
    }

    // read access on UserContext, but delay the release
    {
        ReadAccess<double> read( X, userContext );
        delay = read.releaseDelayed();
    }

    // write access on Host, not possible as release on read is not done
    try
    {
        WriteAccess<double> write( X, hostContext );
    }
    catch ( scai::common::Exception& ex )
    {
        std::cout << "Exception caught: " << ex.what() << std::endl;
    }

    delay();  // now it is okay

    {
        WriteAccess<double> write( X, hostContext );
    }
}
