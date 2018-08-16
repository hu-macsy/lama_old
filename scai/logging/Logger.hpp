/**
 * @file Logger.hpp
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
 * @brief Abstract logger class for a hierarchically organized logging system.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

// for dll import
#include <scai/common/config.hpp>

// local library
#include <scai/logging/SourceLocation.hpp>
#include <scai/logging/Level.hpp>

// std
#include <string>
#include <vector>

namespace scai
{

/** Namespace for logging data structures */

namespace logging
{

/**************************************************************************
 *  class Logger                                                           *
 **************************************************************************/

/** Logger is a an abstract class for hierarchical organization of logging objects

 LoggerProvider is allowed to traverse the hierarchy of loggers and
 to create missing ones (friend class).

 */

class COMMON_DLL_IMPORTEXPORT Logger
{

    friend class LoggerProvider;

public:

    /** Virtual constructor needed due to virtual functions */

    virtual ~Logger();

    /** get the full name of the logger, e.g. "X.Y.Z" */

    std::string getFullName() const;

    /** Check if logging statements of level TRACE are enabled. */

    inline bool isTraceEnabled() const
    {
        return mLevel <= level::TRACE;
    }

    /** Check if logging statements of level DEBUG are enabled. */

    inline bool isDebugEnabled() const
    {
        return mLevel <= level::DEBUG;
    }

    /** Check if logging statements of level INFO are enabled. */

    inline bool isInfoEnabled() const
    {
        return mLevel <= level::INFO;
    }

    /** Check if logging statements of level WARN are enabled. */

    inline bool isWarnEnabled() const
    {
        return mLevel <= level::WARN;
    }

    /** Check if logging statements of level ERROR are enabled. */

    inline bool isErrorEnabled() const
    {
        return mLevel <= level::SERROR;
    }

    /** Check if logging statements of level FATAL are enabled. */

    inline bool isFatalEnabled() const
    {
        return mLevel <= level::FATAL;
    }

    /** Getter routine for the logging level of this object. */

    inline level::Level getLevel() const
    {
        return mLevel;
    }

    /** Getter routine for the effective logging level of
     this object. If level has not been set explicitly it
     will ask the ancestors for the level
     */

    level::Level getEffectiveLevel() const;

    /** Setter routine for the logging level of this object.

     \param level is the new logging level for this object
     \param force true specifies that this is an explicit set

     This routine will set implicitly the levels recursively
     for the descendants whose level has not been set explicitly.
     */

    void setLevel( const level::Level level, const bool force = true );

    /** Logging output for level TRACE. This routine should
     only be called if TraceEnabled() returns true.

     \param loc is the file location of the logging statement
     \param msg is the message to be printed

     Each derived class has to implement this routine. This
     abstract class does not handle output of logging at all.
     */

    virtual void trace( SourceLocation loc, const std::string& msg ) = 0;

    /** Logging output for level DEBUG. This routine should
     only be called if isDebugEnabled() returns true.

     \param loc is the file location of the logging statement
     \param msg is the message to be printed

     Each derived class has to implement this routine. This
     abstract class does not handle output of logging at all.
     */

    virtual void debug( SourceLocation loc, const std::string& msg ) = 0;

    /** Logging output for level INFO. This routine should
     only be called if isInfoEnabled() returns true.

     \param loc is the file location of the logging statement
     \param msg is the message to be printed

     Each derived class has to implement this routine. This
     abstract class does not handle output of logging at all.
     */

    virtual void info( SourceLocation loc, const std::string& msg ) = 0;

    /** Logging output for level WARN. This routine should
     only be called if isWarnEnabled() returns true.

     \param loc is the file location of the logging statement
     \param msg is the message to be printed

     Each derived class has to implement this routine. This
     abstract class does not handle output of logging at all.
     */

    virtual void warn( SourceLocation loc, const std::string& msg ) = 0;

    /** Logging output for level ERROR. This routine should
     only be called if isErrorEnabled() returns true.

     \param loc is the file location of the logging statement
     \param msg is the message to be printed

     Each derived class has to implement this routine. This
     abstract class does not handle output of logging at all.
     */

    virtual void error( SourceLocation loc, const std::string& msg ) = 0;

    /** Logging output for level FATAL. This routine should
     only be called if isFatalEnabled() returns true.

     \param loc is the file location of the logging statement
     \param msg is the message to be printed

     Each derived class has to implement this routine. This
     abstract class does not handle output of logging at all.
     */

    virtual void fatal( SourceLocation loc, const std::string& msg ) = 0;

    /** Predicate that returns true if this logger is the root logger. */

    bool isRootLogger() const;

protected:

    /** Constructor for a logger.
     *
     \param name is the name of this logger (not full name)
     \param parent is pointer to the parent logger
     */

    Logger( const std::string& name, class Logger* parent );

    std::string mName; //!< name of the logger

    bool mSetFlag; //!< This flag indicates that level has been set explicitly.

    level::Level mLevel; //!< Current level of this logger.

    class Logger* const mParent; //!< points back to the parent logger, NULL for root

    std::vector<class Logger*> mSons; //!< sub loggers of this logger

private:

    Logger(); //!< Disable default logger

    Logger( const Logger& other );   // Disable default copy-constructor

    const Logger& operator=( const Logger& other ); // Disable default assignment operator
};

} /* end namespace logging */

} /* end namespace scai */
