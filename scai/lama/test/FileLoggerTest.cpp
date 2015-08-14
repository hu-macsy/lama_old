/**
 * @file FileLoggerTest.cpp
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
 * @brief Contains the implementation of the class FileLoggerTest
 * @author Alexander Büchel, Matthias Makulla
 * @date 02.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/solver/logger/FileLogger.hpp>

#include <fstream>

#include <test/TestMacros.hpp>
#include <test/Configuration.hpp>
#include <scai/common/unique_ptr.hpp>

using namespace lama;
using namespace memory;
using common::Exception;
using common::unique_ptr;
using common::scoped_array;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( FileLoggerTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.FileLoggerTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( LoggingTest )
{
    std::string testMessage( "FileLoggerTestMessage\n" );
    FileLogger& flogger = FileLogger::getFileLogger();
    // This should throw an exception
    LAMA_CHECK_THROW( flogger.setLogFile( "/15/16/17" ), Exception );
    const std::string path = Configuration::getInstance().getPath();
    LAMA_LOG_INFO( logger, "Configuration path = " << path );
    std::string logFileName( path + "/" + "FileLoggerTestFile.log" );
    LAMA_LOG_INFO( logger, "Log file name = " << logFileName );
    flogger.setLogFile( logFileName );
    // Setting same name twice should be okay
    flogger.setLogFile( logFileName );
    // Setting other name should throw an exception
    LAMA_CHECK_THROW( flogger.setLogFile( logFileName + "1" ), Exception );
    flogger.logMessage( testMessage );
    flogger.closeLogFile();
    scoped_array<char> fileInput( new char[testMessage.length()] );
    std::fstream fileStream;
    fileStream.open( logFileName.c_str(), std::fstream::in );
    fileStream.read( fileInput.get(), testMessage.length() );

    for ( IndexType i = 0; i < ( IndexType ) testMessage.length(); ++i )
    {
        BOOST_CHECK_EQUAL( testMessage[i], fileInput[i] );
    }

    fileStream.close();
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
