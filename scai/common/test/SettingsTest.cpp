/**
 * @file SettingsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Test routines for class Settings
 *
 * @author Thomas Brandes
 * @date 05.02.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Settings.hpp>

BOOST_AUTO_TEST_CASE( SettingsTest )
{
    using scai::common::Settings;

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

BOOST_AUTO_TEST_CASE( SettingsConvertTest )
{
    static char var[] = "Dummy";

    using scai::common::Settings;

    bool flag;
    bool set;

    Settings::putEnvironment( var, "Yes" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && flag );
 
    Settings::putEnvironment( var, "No" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && !flag );

    Settings::putEnvironment( var, "shit" );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( !set );

    Settings::putEnvironment( var, 0 );
    set = Settings::getEnvironment( flag, var );
    BOOST_CHECK( set && !flag );
}

