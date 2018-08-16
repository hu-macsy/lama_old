/**
 * @file ExceptionTest.cpp
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
 * @brief Test routines for dervied exception classes.
 * @author Thomas Brandes
 * @date 28.12.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/exception/IOException.hpp>
#include <scai/common/exception/AssertException.hpp>
#include <scai/common/exception/UnsupportedException.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unsupported.hpp>

#include <string>

using namespace scai::common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ExceptionTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( IOExceptionTest )
{
    IOException e( "File error" );
    std::string msg = e.what();
    BOOST_CHECK( msg.find( "IO" ) != std::string::npos );
    BOOST_CHECK( msg.find( "File error" ) != std::string::npos );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( AssertExceptionTest )
{
    AssertException e( "xyz" );
    std::string msg = e.what();
    BOOST_CHECK( msg.find( "Assert" ) != std::string::npos );
    BOOST_CHECK( msg.find( "xyz" ) != std::string::npos );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( UnsupportedExceptionTest )
{
    UnsupportedException e( "zyx" );
    std::string msg = e.what();
    BOOST_CHECK( msg.find( "zyx" ) != std::string::npos );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( UnsupportedTest )
{
    UnsupportedException::resetSetting();

    const char var[] = "SCAI_UNSUPPORTED";

    std::string savedValue;
    bool wasSet = Settings::getEnvironment( savedValue, var );

    bool replace = true;
    Settings::putEnvironment( var, "ERROR", replace );

    BOOST_CHECK_EQUAL( UnsupportedException::getUnsupportedSetting(),
                       UnsupportedException::UNSUPPORTED_ERROR );

    BOOST_CHECK_THROW(
    {
        SCAI_UNSUPPORTED( "Fail for this stuff here" )
    }, UnsupportedException );

    Settings::putEnvironment( var, "IGNORE", replace );

    UnsupportedException::resetSetting();  // otherwise environment variable is not read again

    BOOST_CHECK_EQUAL( UnsupportedException::getUnsupportedSetting(),
                       UnsupportedException::UNSUPPORTED_IGNORE );

    SCAI_UNSUPPORTED( "This stuff is now ignored" )

    Settings::putEnvironment( var, "WARN", replace );

    UnsupportedException::resetSetting();  // otherwise environment variable is not read again

    BOOST_CHECK_EQUAL( UnsupportedException::getUnsupportedSetting(),
                       UnsupportedException::UNSUPPORTED_WARN );

    if ( wasSet )
    {
        Settings::putEnvironment( var, savedValue.c_str(), replace );
    }
    else
    {
        Settings::putEnvironment( var, "", replace );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

