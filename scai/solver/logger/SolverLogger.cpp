/**
 * @file SolverLogger.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of the functionality of the SolverLogger
 * @author Matthias Makulla
 * @date 06.04.2011
 */

// hpp
#include <scai/solver/logger/SolverLogger.hpp>

// local library
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/FileLogger.hpp>
#include <scai/solver/Solver.hpp>

// std
#include <iostream>

namespace scai
{

namespace solver
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
    const lama::Norm& norm,
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

} /* end namespace solver */

} /* end namespace scai */
