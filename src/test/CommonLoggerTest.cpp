/**
 * @file CommonLoggerTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Contains the implementation of the class CommonLoggerTest
 * @author Alexander Büchel, Matthias Makulla
 * @date 02.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <lama/solver/logger/Timer.hpp>
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/solver/logger/FileLogger.hpp>

#include <lama/solver/logger/LoggerWriteBehaviour.hpp>

#include <test/Configuration.hpp>

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CommonLoggerTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.CommonLoggerTest" )

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
                                LoggerWriteBehaviour::toConsoleOnly, std::auto_ptr<Timer>( new Timer() ) );

    BOOST_CHECK_EQUAL( consoleLogger.getLogLevel(), LogLevel::convergenceHistory );

    consoleLogger.setLogLevel( LogLevel::completeInformation );

    BOOST_CHECK_EQUAL( consoleLogger.getLogLevel(), LogLevel::completeInformation );
}

/* --------------------------------------------------------------------- */

void logMessageTest( std::string logFileName, LoggerWriteBehaviour::LoggerWriteBehaviour lwb )
{
    LAMA_LOG_DEBUG( logger, "CommonLoggerTest with LoggerWriteBehaviour: " << lwb );

    CommonLogger consoleAndFileLogger( "<CommonLoggerTest>: ", LogLevel::noLogging, lwb, logFileName,
                                       std::auto_ptr<Timer>( new Timer() ) );

    FileLogger::getFileLogger().setLogFile( logFileName );

    consoleAndFileLogger.logMessage( LogLevel::completeInformation, std::string( "OmittedMessage\n" ) );

    FileLogger::getFileLogger().closeLogFile();
}

BOOST_AUTO_TEST_CASE( ConsoleAndFileLoggingTest )
{
    std::string testMessage( "ConsoleAndFileLoggerTestMessage\n" );
    const std::string path = Configuration::getInstance().getPath();
    LAMA_LOG_INFO( logger, "Configuration path = " << path );
    std::string logFileName( path + "/" + "LogFileCommonLogger.log" );

    logMessageTest( logFileName, LoggerWriteBehaviour::toFileAndConsole );
    logMessageTest( logFileName, LoggerWriteBehaviour::toFileOnly );
    logMessageTest( logFileName, LoggerWriteBehaviour::toConsoleOnly );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
