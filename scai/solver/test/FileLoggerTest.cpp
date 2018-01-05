/**
 * @file FileLoggerTest.cpp
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
 * @brief Test routines for the class FileLogger.
 * @author Matthias Makulla
 * @date 02.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/Communicator.hpp>

#include <scai/solver/logger/FileLogger.hpp>

#include <fstream>

#include <scai/solver/test/TestMacros.hpp>

#include <scai/logging.hpp>

#include <scai/common/test/Configuration.hpp>
#include <scai/common/SCAITypes.hpp>

#include <scai/testsupport/GlobalTempDir.hpp>
#include <scai/testsupport/uniquePathComm.hpp>

#include <memory>

using namespace scai;

using namespace solver;

using common::Exception;
using std::unique_ptr;

using scai::dmemo::Communicator;

using scai::testsupport::uniquePathPerNode;
using scai::testsupport::GlobalTempDir;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( FileLoggerTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.FileLoggerTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( LoggingTest )
{
    const auto logFileName = uniquePathPerNode(GlobalTempDir::getPath(),
                                               *Communicator::getDefaultCommunicatorPtr(),
                                               "FileLoggerTest.log");

    BOOST_TEST_MESSAGE("Log filename: " << logFileName);

    std::string testMessage( "FileLoggerTestMessage\n" );
    FileLogger& flogger = FileLogger::getFileLogger();
    // This should throw an exception
    SCAI_CHECK_THROW( flogger.setLogFile( "/15/16/17" ), Exception );
    flogger.setLogFile( logFileName );
    // Setting same name twice should be okay
    flogger.setLogFile( logFileName );
    // Setting other name should throw an exception
    SCAI_CHECK_THROW( flogger.setLogFile( logFileName + "1" ), Exception );
    flogger.logMessage( testMessage );
    flogger.closeLogFile();
    unique_ptr<char[]> fileInput( new char[testMessage.length()] );
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
