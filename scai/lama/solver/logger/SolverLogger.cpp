/**
 * @file SolverLogger.cpp
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
 * @brief Implementation of the functionality of the SolverLogger
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/logger/SolverLogger.hpp>

// local library
#include <scai/lama/solver/logger/Timer.hpp>
#include <scai/lama/solver/logger/FileLogger.hpp>
#include <scai/lama/solver/Solver.hpp>

// std
#pragma offload_attribute (push, target(mic))
#include <iostream>
#pragma offload_attribute (pop)

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( SolverLogger::logger, "SolverLogger" );

const std::string& SolverLogger::id() const
{
    //static const std::string id = "SolverLogger";
    return mId;
}

SolverLogger::SolverLogger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    bool ignoreRank )
    : mId( id ), mLogLevel( level ), mWriteBehaviour( writeBehaviour ), mIgnoreRank( ignoreRank )
{
    mTimer.reset( new Timer() );

    SCAI_LOG_INFO( logger, "SolverLogger created, id = " << mId << ", level = " << mLogLevel << ", writeBehaviour = "
                           << mWriteBehaviour << ", ignore rank = " << mIgnoreRank )
}

SolverLogger::SolverLogger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    common::shared_ptr<Timer> timer,
    bool ignoreRank )
    : mTimer( timer ), mId( id ), mLogLevel( level ), mWriteBehaviour( writeBehaviour ), mIgnoreRank(
          ignoreRank )
{
    SCAI_LOG_INFO( logger, "SolverLogger with timer created, id = " << mId << ", level = " << mLogLevel << ", writeBehaviour = "
                           << mWriteBehaviour << ", ignore rank = " << mIgnoreRank )
}

SolverLogger::SolverLogger(
    const std::string& id,
    LogLevel::LogLevel level,
    LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
    const std::string& logFileName,
    common::shared_ptr<Timer> timer,
    bool ignoreRank )
    : mTimer( timer ), mId( id ), mLogLevel( level ), mWriteBehaviour( writeBehaviour ), mIgnoreRank(
          ignoreRank )
{
    FileLogger::getFileLogger().setLogFile( logFileName );

    SCAI_LOG_INFO( logger, "SolverLogger for file " << logFileName << " created, id = " << mId
                           << ", level = " << mLogLevel << ", writeBehaviour = " 
                           << mWriteBehaviour << ", ignore rank = " << mIgnoreRank )
}

SolverLogger::~SolverLogger()
{
}

LogLevel::LogLevel SolverLogger::getLogLevel() const
{
    return mLogLevel;
}

void SolverLogger::setLogLevel( LogLevel::LogLevel level )
{
    mLogLevel = level;
}

void SolverLogger::logString( const std::string& message )
{
    std::string logMessage = this->createPrefix() + message;

    switch( mWriteBehaviour )
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

        default:
 
            SCAI_LOG_ERROR( logger, "illegal write behaviour: " << mWriteBehaviour )
    }
}

void SolverLogger::logString( LogLevel::LogLevel level, const std::string& message )
{
    if( level <= mLogLevel && !mIgnoreRank )
    {
        logString( message );
    }
}

void SolverLogger::logMessage( LogLevel::LogLevel level, const std::string& message )
{
    SCAI_LOG_DEBUG( logger, "logMessage, level = " << level 
                             << " ( mLevel = " << mLogLevel << " ), msg = " << message )
    logString( level, message );
}

void SolverLogger::logNewLine( LogLevel::LogLevel level )
{
    logString( level, "\n" );
}

void SolverLogger::logResidual(
    LogLevel::LogLevel level,
    const Solver& solver,
    const Norm& norm,
    const std::string iterationPrefix )
{
    if( level <= mLogLevel )
    {
        std::stringstream residualStream;
        residualStream << iterationPrefix;
        residualStream << "Residual: ";
        residualStream << norm( solver.getResidual() );
        residualStream << "\n";
        logString( level, residualStream.str() );
    }
}

void SolverLogger::logTime( const std::string& timerId, LogLevel::LogLevel level, const std::string& message )
{
    SCAI_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );

    double time = mTimer->getTime( timerId );

    if( level <= mLogLevel )
    {
        std::stringstream timeStringStream;
        timeStringStream << message;
        timeStringStream << time;
        timeStringStream << "\n";
        logString( level, timeStringStream.str() );
    }
}

void SolverLogger::startTimer( const std::string& timerId )
{
    SCAI_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );

    mTimer->start( timerId );
}

void SolverLogger::stopTimer( const std::string& timerId )
{
    SCAI_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );

    mTimer->stop( timerId );
}

void SolverLogger::stopAndResetTimer( const std::string& timerId )
{
    SCAI_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );

    mTimer->stopAndReset( timerId );
}

} /* end namespace lama */

} /* end namespace scai */
