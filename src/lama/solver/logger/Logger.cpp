/**
 * @file Logger.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Implementation of the functionality of the Logger
 * @author Matthias Makulla
 * @date 06.04.2011
 * $Id$
 */

// hpp
#include <lama/solver/logger/Logger.hpp>

// default timer
#include <lama/solver/logger/Timer.hpp>

// others
#include <lama/solver/logger/FileLogger.hpp>
#include <lama/solver/Solver.hpp>

#include <iostream>

namespace lama
{

LAMA_LOG_DEF_LOGGER( Logger::logger, "Logger" );

const std::string& Logger::id() const
{
    //static const std::string id = "Logger";
    return mId;
}

Logger::Logger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    bool ignoreRank )
    : mId( id ), mLogLevel( level ), mWriteBehaviour( writeBehaviour ), mIgnoreRank( ignoreRank )
{
    mTimer.reset( new Timer() );
}

Logger::Logger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    std::auto_ptr<Timer> timer,
    bool ignoreRank )
    : mTimer( timer ), mId( id ), mLogLevel( level ), mWriteBehaviour( writeBehaviour ), mIgnoreRank(
        ignoreRank )
{
}

Logger::Logger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    const std::string& logFileName,
    std::auto_ptr<Timer> timer,
    bool ignoreRank )
    : mTimer( timer ), mId( id ), mLogLevel( level ), mWriteBehaviour( writeBehaviour ), mIgnoreRank(
        ignoreRank )
{
    FileLogger::getFileLogger().setLogFile( logFileName );
}

Logger::~Logger()
{
}

LogLevel::LogLevel Logger::getLogLevel() const
{
    return mLogLevel;
}

void Logger::setLogLevel( LogLevel::LogLevel level )
{
    mLogLevel = level;
}

void Logger::logString( const std::string& message )
{
    std::string logMessage = this->createPrefix() + message;

    switch ( mWriteBehaviour )
    {
    case LoggerWriteBehaviour::toFileAndConsole:
    {
        FileLogger::getFileLogger().logMessage( logMessage );
        std::cout << logMessage;
    }
    break;
    case LoggerWriteBehaviour::toConsoleOnly:
    {
        std::cout << logMessage;
    }
    break;
    case LoggerWriteBehaviour::toFileOnly:
    {
        FileLogger::getFileLogger().logMessage( logMessage );
    }
    break;
    }
}

void Logger::logString( LogLevel::LogLevel level, const std::string& message )
{
    if ( level <= mLogLevel && !mIgnoreRank )
    {
        logString( message );
    }
}

void Logger::logMessage( LogLevel::LogLevel level, const std::string& message )
{
    logString( level, message );
}

void Logger::logNewLine( LogLevel::LogLevel level )
{
    logString( level, "\n" );
}

void Logger::logResidual(
    LogLevel::LogLevel level,
    const Solver& solver,
    const Norm& norm,
    const std::string iterationPrefix )
{
    if ( level <= mLogLevel )
    {
        std::stringstream residualStream;
        residualStream << iterationPrefix;
        residualStream << "Residual: ";
        residualStream << norm( solver.getResidual() );
        residualStream << "\n";
        logString( level, residualStream.str() );
    }
}

void Logger::logTime( const std::string& timerId, LogLevel::LogLevel level, const std::string& message )
{
    LAMA_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );

    double time = mTimer->getTime( timerId );

    if ( level <= mLogLevel )
    {
        std::stringstream timeStringStream;
        timeStringStream << message;
        timeStringStream << time;
        timeStringStream << "\n";
        logString( level, timeStringStream.str() );
    }
}

void Logger::startTimer( const std::string& timerId )
{
    LAMA_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );

    mTimer->start( timerId );
}

void Logger::stopTimer( const std::string& timerId )
{
    LAMA_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );

    mTimer->stop( timerId );
}

void Logger::stopAndResetTimer( const std::string& timerId )
{
    LAMA_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );

    mTimer->stopAndReset( timerId );
}

} // namespace lama
