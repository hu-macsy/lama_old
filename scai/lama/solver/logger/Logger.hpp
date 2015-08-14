/**
 * @file Logger.hpp
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
 * @brief contains the Logger header.
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */
#ifndef LAMA_SOLVER_LOGGER_HPP_
#define LAMA_SOLVER_LOGGER_HPP_

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>
#include <scai/common/shared_ptr.hpp>

// others
#include <scai/lama/solver/logger/LogLevel.hpp>
#include <scai/lama/solver/logger/LoggerWriteBehaviour.hpp>
#include <scai/lama/solver/logger/Timer.hpp>

#include <scai/lama/norm/Norm.hpp>

#include <string>
#include <sstream>
#include <memory>

namespace lama
{

class Solver;
class Logger;

typedef common::shared_ptr<Logger> LoggerPtr;

/**
 * @brief A logger abstraction.
 *
 * This class represents a logger abstraction. It defines common logger
 * operations. Derived classes may use the createPrefix() method to
 * customize messages.
 */
class COMMON_DLL_IMPORTEXPORT Logger: private common::NonCopyable
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
    Logger(
        const std::string& id,
        LogLevel::LogLevel level,
        LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
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
    Logger(
        const std::string& id,
        LogLevel::LogLevel level,
        LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
        common::shared_ptr<Timer> timer,
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
    Logger(
        const std::string& id,
        LogLevel::LogLevel level,
        LoggerWriteBehaviour::LoggerWriteBehaviour writeBehaviour,
        const std::string& logFileName,
        common::shared_ptr<Timer> timer,
        bool ignoreRank = false );

    /**
     * Destructor.
     */
    virtual ~Logger();

    /**
     * @brief Returns the log level of this logger.
     *
     * @return   The log level of this.
     */
    LogLevel::LogLevel getLogLevel() const;

    /**
     * @brief Sets the internal LogLevel
     *
     * @param level   The log level.
     */
    void setLogLevel( LogLevel::LogLevel level );

    /**
     * @brief Logs a user specified message
     *
     * This method logs a user specified message (string)
     *
     * @param level     The LogLevel of the message.
     * @param message   The message to log.
     */
    void logMessage( LogLevel::LogLevel level, const std::string& message );

    /**
     * @brief Logs an empty line.
     *
     * @param level   The LogLevel at which the empty line shall be logged
     */
    void logNewLine( LogLevel::LogLevel level );

    /**
     * @brief Logs the residual of the solver.
     *
     * @param[in] level             The LogLevel at which the residual shall be logged.
     * @param[in] solver            The solver which supplies the residual
     * @param[in] norm              The Norm used to calculate the residual
     * @param[in] iterationPrefix   A Prefix to put in front of the generated log message (Default: "" )
     */
    void logResidual(
        LogLevel::LogLevel level,
        const Solver& solver,
        const Norm& norm,
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
    void logTime( const std::string& timerId, LogLevel::LogLevel level, const std::string& message );

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
    void logType( LogLevel::LogLevel level, const std::string& message, ValueType arg );

protected:

    /**
     * @brief Performs a LogLevel check and logs the message depending on this
     *
     * @param[in] level     The LogLevel to log the string at.
     * @param[in] message   The string to be logged.
     */
    virtual void logString( LogLevel::LogLevel level, const std::string& message );

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
    common::shared_ptr<Timer> mTimer;

    std::string mId;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private    :
    LogLevel::LogLevel mLogLevel;
    LoggerWriteBehaviour::LoggerWriteBehaviour mWriteBehaviour;
    bool mIgnoreRank;
};

template<typename ValueType>
void Logger::logType( LogLevel::LogLevel level, const std::string& message, ValueType arg )
{
    if( level <= mLogLevel )
    {
        std::stringstream intStream;
        intStream << message;
        intStream << arg;
        intStream << "\n";

        logString( level, intStream.str() );
    }
}

} // namespace lama

#endif // LAMA_SOLVER_LOGGER_HPP_
