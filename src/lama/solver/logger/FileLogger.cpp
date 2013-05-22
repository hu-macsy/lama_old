/**
 * @file FileLogger.cpp
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
 * @brief This file contains the implementation of the FileLogger.hpp file
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

// hpp
#include <lama/solver/logger/FileLogger.hpp>

// others
#include <lama/exception/Exception.hpp>

#include <sstream>

namespace lama
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
    mFileStream << message;
}

void FileLogger::setLogFile( const std::string& logFileName )
{
    if ( !mFileStream.is_open() )
    {
        mFileStream.open( logFileName.c_str(), std::fstream::out );
        if ( mFileStream.fail() )
        {
            LAMA_THROWEXCEPTION( "Could not open log file " << logFileName );
        }
        mFileName = logFileName;
    }
    else if ( logFileName != mFileName )
    {
        LAMA_THROWEXCEPTION( "Tried to set the log file of the logger to two different files." );
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

} // namespace lama
