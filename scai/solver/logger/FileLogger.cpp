/**
 * @file FileLogger.cpp
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
 * @brief This file contains the implementation of the FileLogger.hpp file
 * @author Matthias Makulla
 * @date 06.04.2011
 */

// hpp
#include <scai/solver/logger/FileLogger.hpp>

// internal scai libraries
#include <scai/common/macros/throw.hpp>

// std
#include <sstream>

namespace scai
{

namespace solver
{

FileLogger::~FileLogger()
{
    if ( mFileStream.is_open() )
    {
        mFileStream.close();
    }
}

FileLogger& FileLogger::getFileLogger()
{
    static FileLogger logger;
    return logger;
}

void FileLogger::logMessage( const std::string& message )
{
    if ( mFileStream.is_open() )
    {
        mFileStream << message;
    }
    else
    {
        COMMON_THROWEXCEPTION( "FileLogger: no open file" )
    }
}

void FileLogger::setLogFile( const std::string& logFileName )
{
    if ( !mFileStream.is_open() )
    {
        mFileStream.open( logFileName.c_str(), std::fstream::out );

        if ( mFileStream.fail() )
        {
            COMMON_THROWEXCEPTION( "Could not open log file " << logFileName );
        }

        mFileName = logFileName;
    }
    else if ( logFileName != mFileName )
    {
        COMMON_THROWEXCEPTION( "Tried to set new log file " << logFileName << ", but current log file " << mFileName << " is still open" );
    }
}

void FileLogger::closeLogFile()
{
    if ( mFileStream.is_open() )
    {
        mFileStream.close();
        mFileName = "";
    }
}

FileLogger::FileLogger()
{
}

} /* end namespace solver */

} /* end namespace scai */
