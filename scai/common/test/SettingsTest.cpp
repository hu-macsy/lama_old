#include <boost/test/unit_test.hpp>

#include <scai/common/Settings.hpp>

BOOST_AUTO_TEST_CASE( SettingsTest )
{

using scai::common::Settings;

const char* args[] = { "--SCAI_DEVICE=0,1,2", "--SOLVER=cg" };

int nargs = 2;

Settings::parseArgs( nargs, args );

BOOST_CHECK_EQUAL( 1, nargs );

int device = -1;

Settings::setRank( 3 );   

bool set = Settings::getEnvironment( device, "SCAI_DEVICE" );

BOOST_CHECK( set );
BOOST_CHECK_EQUAL( 0, device );

Settings::putEnvironment( "SCAI_DEVICE", "dummy" );

set = Settings::getEnvironment( device, "SCAI_DEVICE" );
BOOST_CHECK( !set );

}

