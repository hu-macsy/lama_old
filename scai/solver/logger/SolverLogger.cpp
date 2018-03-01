/**
 * @file SolverLogger.cpp
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

#include <scai/common/macros/loop.hpp>

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
    LogLevel level,
    LoggerWriteBehaviour writeBehaviour,
    bool ignoreRank )
    : mId( id ), mLogLevel( level ), mWriteBehaviour( writeBehaviour ), mIgnoreRank( ignoreRank )
{
    mTimer.reset( new Timer() );
    SCAI_LOG_INFO( logger, "SolverLogger created, id = " << mId << ", level = " << mLogLevel << ", writeBehaviour = "
                   << mWriteBehaviour << ", ignore rank = " << mIgnoreRank )
}

SolverLogger::SolverLogger(
    const std::string& id,
    LogLevel level,
    LoggerWriteBehaviour writeBehaviour,
    std::shared_ptr<Timer> timer,
    bool ignoreRank )
    : mTimer( timer ), mId( id ), mLogLevel( level ), mWriteBehaviour( writeBehaviour ), mIgnoreRank(
        ignoreRank )
{
    SCAI_LOG_INFO( logger, "SolverLogger with timer created, id = " << mId << ", level = " << mLogLevel << ", writeBehaviour = "
                   << mWriteBehaviour << ", ignore rank = " << mIgnoreRank )
}

SolverLogger::SolverLogger(
    const std::string& id,
    LogLevel level,
    LoggerWriteBehaviour writeBehaviour,
    const std::string& logFileName,
    std::shared_ptr<Timer> timer,
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

LogLevel SolverLogger::getLogLevel() const
{
    return mLogLevel;
}

void SolverLogger::setLogLevel( LogLevel level )
{
    mLogLevel = level;
}

void SolverLogger::logString( const std::string& message )
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

        default:
            SCAI_LOG_ERROR( logger, "illegal write behaviour: " << mWriteBehaviour )
    }
}

void SolverLogger::logString( LogLevel level, const std::string& message )
{
    if ( level <= mLogLevel && !mIgnoreRank )
    {
        logString( message );
    }
}

void SolverLogger::logMessage( LogLevel level, const std::string& message )
{
    SCAI_LOG_DEBUG( logger, "logMessage, level = " << level
                    << " ( mLevel = " << mLogLevel << " ), msg = " << message )
    logString( level, message );
}

void SolverLogger::logNewLine( LogLevel level )
{
    logString( level, "\n" );
}

template<typename ValueType>
void SolverLogger::logResidual(
    LogLevel level,
    const Solver<ValueType>& solver,
    const lama::Norm<ValueType>& norm,
    const std::string iterationPrefix )
{
    if ( level <= mLogLevel )
    {
        const char* typeId = common::TypeTraits<ValueType>::id();
        std::stringstream residualStream;
        residualStream << iterationPrefix;
        residualStream << "Residual<" << typeId << "> : ";
        residualStream << norm( solver.getResidual() );
        residualStream << "\n";
        logString( level, residualStream.str() );
    }
}

void SolverLogger::logTime( const std::string& timerId, LogLevel level, const std::string& message )
{
    SCAI_ASSERT_DEBUG( mTimer.get(), "mTimer == NULL" );
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

/* ========================================================================= */
/*       Template method instantiation                                       */
/* ========================================================================= */

#define SOLVER_LOGGER_SPECIFIER( ValueType )             \
     template void SolverLogger::logResidual(            \
        LogLevel level,                                  \
        const Solver<ValueType>& solver,                 \
        const lama::Norm<ValueType>&,                    \
        const std::string );                             \

SCAI_COMMON_LOOP( SOLVER_LOGGER_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
