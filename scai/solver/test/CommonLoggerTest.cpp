/**
 * @file CommonLoggerTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Test routines for the class CommonLogger.
 * @author Matthias Makulla
 * @date 02.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/logger/FileLogger.hpp>

#include <scai/solver/logger/LoggerWriteBehaviour.hpp>

#include <scai/common/test/Configuration.hpp>

#include <memory>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

using std::shared_ptr;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CommonLoggerTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CommonLoggerTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( LoggerIdTest )
{
    CommonLogger consoleLogger( "<CommonLoggerTest>: ", LogLevel::convergenceHistory,
                                LoggerWriteBehaviour::toConsoleOnly );
    std::string s = consoleLogger.id();
    BOOST_CHECK_EQUAL( s, "<CommonLoggerTest>: " );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( setAndGetLogLevelTest )
{
    CommonLogger consoleLogger( "<CommonLoggerTest>: ", LogLevel::convergenceHistory,
                                LoggerWriteBehaviour::toConsoleOnly, shared_ptr<Timer>( new Timer() ) );
    BOOST_CHECK_EQUAL( consoleLogger.getLogLevel(), LogLevel::convergenceHistory );
    consoleLogger.setLogLevel( LogLevel::completeInformation );
    BOOST_CHECK_EQUAL( consoleLogger.getLogLevel(), LogLevel::completeInformation );
}

/* --------------------------------------------------------------------- */

void logMessageTest( std::string logFileName, LoggerWriteBehaviour lwb )
{
    SCAI_LOG_DEBUG( logger, "CommonLoggerTest with LoggerWriteBehaviour: " << lwb );
    CommonLogger consoleAndFileLogger( "<CommonLoggerTest>: ", LogLevel::noLogging, lwb, logFileName,
                                       shared_ptr<Timer>( new Timer() ) );
    FileLogger::getFileLogger().setLogFile( logFileName );
    consoleAndFileLogger.logMessage( LogLevel::completeInformation, std::string( "OmittedMessage\n" ) );
    FileLogger::getFileLogger().closeLogFile();
}

BOOST_AUTO_TEST_CASE( ConsoleAndFileLoggingTest )
{
    std::string testMessage( "ConsoleAndFileLoggerTestMessage\n" );
    const std::string path = scai::test::Configuration::getPath();
    SCAI_LOG_INFO( logger, "Configuration path = " << path );
    std::string logFileName( path + "/" + "LogFileCommonLogger.log" );
    logMessageTest( logFileName, LoggerWriteBehaviour::toFileAndConsole );
    logMessageTest( logFileName, LoggerWriteBehaviour::toFileOnly );
    logMessageTest( logFileName, LoggerWriteBehaviour::toConsoleOnly );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
