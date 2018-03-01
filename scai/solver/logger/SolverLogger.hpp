/**
 * @file SolverLogger.hpp
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
 * @brief contains the SolverLogger header.
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>

// local library
#include <scai/solver/logger/LogLevel.hpp>
#include <scai/solver/logger/LoggerWriteBehaviour.hpp>
#include <scai/solver/logger/Timer.hpp>

// internal scai libraries
#include <scai/lama/norm/Norm.hpp>

// std
#include <string>
#include <sstream>
#include <memory>

namespace scai
{

namespace solver
{

template<typename ValueType> class Solver;
class SolverLogger;

typedef std::shared_ptr<SolverLogger> LoggerPtr;

/**
 * @brief A logger abstraction.
 *
 * This class represents a logger abstraction. It defines common logger
 * operations. Derived classes may use the createPrefix() method to
 * customize messages.
 */
class COMMON_DLL_IMPORTEXPORT SolverLogger: private common::NonCopyable
{
public:

    /**
     * @brief Returns the id of this.
     *
     * @return the id of this.
     */
    const std::string& id() const;

    /**
     * @brief Creates a logger with the specified properties.
     *
     * This constructor creates a common logger with the specified
     * properties. These properties are the write behavior and the
     * log-level. No file name for a logfile is required.
     *
     * @param[in] id               TODO[doxy] Complete Description.
     * @param[in] level            The loglevel of the logger. Messages with a loglevel
     *                             greater than the level of the logger will be omitted.
     *                             Instead of a common logger with the "noLogging"
     *                             loglevel a NullLogger should be used.
     * @param[in] writeBehaviour   Specifies, if the logger shall write its
     *                             output to the console and a file or to the
     *                             console only
     * @param[in] ignoreRank       TODO[doxy] Complete Description.
     */
    SolverLogger(
        const std::string& id,
        LogLevel level,
        LoggerWriteBehaviour writeBehaviour,
        bool ignoreRank = false );

    /**
     * @brief Creates a logger with the specified properties.
     *
     * This constructor creates a common logger with the specified
     * properties. These properties are the write behavior and the
     * log-level. No file name for a logfile is required.
     *
     * @param[in] id               TODO[doxy] Complete Description.
     * @param[in] level            The loglevel of the logger. Messages with a loglevel
     *                             greater than the level of the logger will be omitted.
     *                             Instead of a common logger with the "noLogging"
     *                             loglevel a NullLogger should be used.
     * @param[in] writeBehaviour   Specifies, if the logger shall write its
     *                             output to the console and a file or to the
     *                             console only
     *
     * @param[in] timer            The timer which shall be used by this logger.
     * @param[in] ignoreRank       TODO[doxy] Complete Description.
     */
    SolverLogger(
        const std::string& id,
        LogLevel level,
        LoggerWriteBehaviour writeBehaviour,
        std::shared_ptr<Timer> timer,
        bool ignoreRank = false );

    /**
     * @brief Creates a logger with the specified properties.
     *
     * This constructor creates a common logger with the specified
     * properties. These properties are the write behavior and the
     * log-level. No file name for a logfile is required.
     *
     * @param[in] id               TODO[doxy] Complete Description.
     * @param[in] level            The loglevel of the logger. Messages with a loglevel
     *                             greater than the level of the logger will be omitted.
     *                             Instead of a common logger with the "noLogging"
     *                             loglevel a NullLogger should be used.
     * @param[in] writeBehaviour   Specifies, if the logger shall write its
     *                             output to the console and a file or to the
     *                             console only
     * @param[in] logFileName      The name of the logfile which shall be used by
     *                             the logger. WARNING: If the users wants so use
     *                             multiple loggers in different solvers this
     *                             constructor should only be used once (for
     *                             example with the top-level solver). The other
     *                             solvers should than use the other constructor
     *                             for the common logger.
     * @param[in] timer            The timer which shall be used by this logger.
     * @param[in] ignoreRank       TODO[doxy] Complete Description.
     */
    SolverLogger(
        const std::string& id,
        LogLevel level,
        LoggerWriteBehaviour writeBehaviour,
        const std::string& logFileName,
        std::shared_ptr<Timer> timer,
        bool ignoreRank = false );

    /**
     * Destructor.
     */
    virtual ~SolverLogger();

    /**
     * @brief Returns the log level of this logger.
     *
     * @return   The log level of this.
     */
    LogLevel getLogLevel() const;

    /**
     * @brief Sets the internal LogLevel
     *
     * @param level   The log level.
     */
    void setLogLevel( LogLevel level );

    /**
     * @brief Logs a user specified message
     *
     * This method logs a user specified message (string)
     *
     * @param level     The LogLevel of the message.
     * @param message   The message to log.
     */
    void logMessage( LogLevel level, const std::string& message );

    /**
     * @brief Logs an empty line.
     *
     * @param level   The LogLevel at which the empty line shall be logged
     */
    void logNewLine( LogLevel level );

    /**
     * @brief Logs the residual of the solver.
     *
     * @param[in] level             The LogLevel at which the residual shall be logged.
     * @param[in] solver            The solver which supplies the residual
     * @param[in] norm              The Norm used to calculate the residual
     * @param[in] iterationPrefix   A Prefix to put in front of the generated log message (Default: "" )
     */
    template<typename ValueType>
    void logResidual(
        LogLevel level,
        const Solver<ValueType>& solver,
        const lama::Norm<ValueType>& norm,
        const std::string iterationPrefix = "" );

    /**
     * @brief logs the elapsed time in combination with a user-specified message.
     *
     * Logs the elapsed time since the internal timer got started, does not
     * stop the timer. The caller may also specify an additional message.
     * The time and message will be logged in the format [message] [time].
     *
     * @param[in] timerId   The ID of the timer whoms time shall be logged.
     * @param[in] level     The loglevel of the message.
     * @param[in] message   The message to log.
     */
    void logTime( const std::string& timerId, LogLevel level, const std::string& message );

    /**
     * @brief Starts an internal timer of the logger.
     *
     * Starts the internal timer of the logger.
     *
     * @param[in] timerId The ID of the timer which shall be started.
     */
    void startTimer( const std::string& timerId );

    /**
     * @brief Stops an internal timer of the logger. Does not reset it.
     *
     * Stops an internal timer of the logger. Does not reset the timer.
     * The caller may call startTimer() again and the timer will continue
     * to measure the time, beginning with the amount of time it had
     * measured, before the timer got stopped.
     *
     * @param[in] timerId The ID of the timer which shall be stopped.
     */
    void stopTimer( const std::string& timerId );

    /**
     * @brief Stops and resets an internal timer.
     *
     * Stops and resets an internal timer.
     *
     * @param[in] timerId The ID of the timer which shall be stopped and reset.
     */
    void stopAndResetTimer( const std::string& timerId );

    /**
     * @brief logs a message followed by an custom type by using the << operator.
     *
     * Example: Can be used to log the number of iterations
     * made by an iterative solver.
     *
     * @param[in] level     The LogLevel
     * @param[in] message   The message which shall be logged
     * @param[in] arg       The argument which shall be logged after the message.
     *                      The << operator must be defined for it!
     */
    template<typename ValueType>
    void logType( LogLevel level, const std::string& message, ValueType arg );

protected:

    /**
     * @brief Performs a LogLevel check and logs the message depending on this
     *
     * @param[in] level     The LogLevel to log the string at.
     * @param[in] message   The string to be logged.
     */
    virtual void logString( LogLevel level, const std::string& message );

    /**
     * @brief Logs the message. DOES NOT perform a LogLevel check.
     *
     * @param[in] message   The message to be logged.
     */
    virtual void logString( const std::string& message );

    /**
     * @brief Template-Pattern Method - used by subclasses to create
     * a prefix which will be added to each message.
     */
    virtual std::string createPrefix() = 0;

    /**
     * @brief Timer used for timings.
     */
    std::shared_ptr<Timer> mTimer;

    std::string mId;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private    :
    LogLevel mLogLevel;
    LoggerWriteBehaviour mWriteBehaviour;
    bool mIgnoreRank;
};

template<typename ValueType>
void SolverLogger::logType( LogLevel level, const std::string& message, ValueType arg )
{
    if ( level <= mLogLevel )
    {
        std::stringstream intStream;
        intStream << message;
        intStream << arg;
        intStream << "\n";
        logString( level, intStream.str() );
    }
}

} /* end namespace solver */

} /* end namespace scai */
