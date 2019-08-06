/**
 * @file SettingsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Test routines for class Settings
 * @author Thomas Brandes
 * @date 05.02.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/exception/Exception.hpp>
#include <fstream>

using scai::common::Settings;
using scai::common::Exception;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SettingsTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( parseArgTest )
{
    const char* args[] = { "--SCAI_DEVICE=0,1,2", "--SOLVER=cg" };
    int nargs = 2;
    Settings::parseArgs( nargs, args );
    // exactly one argument is taken for parsing
    BOOST_CHECK_EQUAL( 1, nargs );
    int device = -1;
    Settings::setRank( 3 );  // choice for multiple arguments separated by ,
    bool set = Settings::getEnvironment( device, "SCAI_DEVICE" );
    BOOST_CHECK( set );
    BOOST_CHECK_EQUAL( 0, device );
    Settings::putEnvironment( "SCAI_DEVICE", "dummy" );
    set = Settings::getEnvironment( device, "SCAI_DEVICE" );
    BOOST_CHECK( !set );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getBoolTest )
{
    static char var[] = "Dummy";
    bool flag;
    bool set;
    Settings::putEnvironment( var, "Yes" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && flag );
    Settings::putEnvironment( var, "No" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && !flag );
    Settings::putEnvironment( var, "1" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && flag );
    Settings::putEnvironment( var, "false" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && !flag );
    Settings::putEnvironment( var, "True" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && flag );
    Settings::putEnvironment( var, "0" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && !flag );
    Settings::putEnvironment( var, "Ja" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && flag );
    Settings::putEnvironment( var, "shit" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( !set );
    Settings::putEnvironment( var, 0 );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && !flag );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( printSettingTest )
{

    static std::string var  = "SCAI_XYZ";
    static std::string val  = "XYZ_Done";

    Settings::putEnvironment( var.c_str(), val.c_str() );

    std::ostringstream out;

    Settings::printEnvironment( out );

    BOOST_CHECK( out.str().find( var ) != std::string::npos );
    BOOST_CHECK( out.str().find( val ) != std::string::npos );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getIntTest )
{
    Settings::putEnvironment( "SCAI_INT_VAL", "1" );

    {
        int i = 0;
        int one = 1;
        bool set = Settings::getEnvironment( i, "SCAI_INT_VAL" );
        BOOST_CHECK( set );
        BOOST_CHECK_EQUAL( one, i );
    }
    {
        unsigned int i = 0;
        unsigned int one = 1;
        bool set = Settings::getEnvironment( i, "SCAI_INT_VAL" );
        BOOST_CHECK( set );
        BOOST_CHECK_EQUAL( one, i );
    }
    {
        long i = 0;;
        bool set = Settings::getEnvironment( i, "SCAI_INT_VAL" );
        BOOST_CHECK( set );
        long one = 1;
        BOOST_CHECK_EQUAL( one, i );
    }
    {
        unsigned long i = 0;
        unsigned long one = 1;
        bool set = Settings::getEnvironment( i, "SCAI_INT_VAL" );
        BOOST_CHECK( set );
        BOOST_CHECK_EQUAL( one, i );
    }

    /* --------------------------------------------------------------------- */

}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( settingsFileTest )
{
     // generate an example file

     std::fstream file;
     file.open( "test_settings.txt", std::ios::out | std::ios::trunc );
     file << "dummy0 0 WEIGHT=0.4\n";
     file << "error1 5 WEIGHT DOMAIN=5\n";
     file << "error2 X DOMAIN=5\n";
     file << "dummy1 5 WEIGHT=0.3 DOMAIN=5\n";
     file << "dummy0 0-2 DOMAIN=4\n";
     file.close();

     Settings::readSettingsFile( "test_settings.txt", "dummy0", 0 );

     float weight = 0;
     float expWeight = 0.4f;
     bool set = Settings::getEnvironment( weight, "WEIGHT" );
     BOOST_CHECK( set );
     BOOST_CHECK_EQUAL( weight, expWeight );

     int domain = 0;
     int expDomain = 4;
     set = Settings::getEnvironment( domain, "DOMAIN" );
     BOOST_CHECK( set );
     BOOST_CHECK_EQUAL( domain, expDomain );

     BOOST_CHECK_THROW( 
     {
         Settings::readSettingsFile( "test_settings.txt", "error1", 5 );
     }, Exception );

     BOOST_CHECK_THROW( 
     {
         Settings::readSettingsFile( "test_settings.txt", "error2", 5 );
     }, Exception );

     BOOST_CHECK_THROW( 
     {
         Settings::readSettingsFile( "test_settings.txt", "none", 0 );
     }, Exception );

     Settings::readSettingsFile( "test_settings.txt", "dummy1", 5 );

     expWeight = 0.3f;
     set = Settings::getEnvironment( weight, "WEIGHT" );
     BOOST_CHECK( set );
     BOOST_CHECK_EQUAL( weight, expWeight );
}

BOOST_AUTO_TEST_SUITE_END()
