/**
 * @file GenLogger.hpp
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
 * @brief Generic logger class that writes logging messages to stdout.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

// base classes
#include <logging/Logger.hpp>

namespace logging
{

/** GenLogger is a simple and generic solution of the abstract Logger class */

class GenLogger: public Logger
{

public:

    /** GenLoggerCreator is friend as constructor is only private. */

    friend class GenLoggerCreator;

    virtual ~GenLogger()
    {
    }

    /** Implementation of Logger::trace */

    virtual void trace( SourceLocation loc, const std::string& msg );

    /** Implementation of Logger::debug */

    virtual void debug( SourceLocation loc, const std::string& msg );

    /** Implementation of Logger::info */

    virtual void info( SourceLocation loc, const std::string& msg );

    /** Implementation of Logger::warn */

    virtual void warn( SourceLocation loc, const std::string& msg );

    /** Implementation of Logger::error */

    virtual void error( SourceLocation loc, const std::string& msg );

    /** Implementation of Logger::fatal */

    virtual void fatal( SourceLocation loc, const std::string& msg );

    /** Provide the root logger for this class.*/

    static Logger& getRoot();

    /** Enables or disables flushing */

    static void setFlush( bool flush );

    static void setFormat( const std::string& format );

    /** Output routine as global variable so it might be reset for other purposes. */

    static void ( *myPrintf ) ( const char* format, ... );

private:

    static bool sFlush; //!< if true flush each loggging output

    static std::vector<std::string> formatTokens;  //!< tokens of ouputline

    /** Constructor for a generic logger.
     *
     *  The constructor is private, only GenLoggerCreator can create it.

     \param name is the name of the logger at this level
     \param parent is a pointer to the ancestor logger (NULL for root)

     \sa Logger::Logger
     */

    GenLogger( const std::string& name, class Logger* parent );

    /** Generic routine for logging output that can be used for all levels

     \param level is the string representation of the level
     \param loc is the file location of the logging statement
     \param msg is the message output of the logging statement
     */

    void log( const char* level, SourceLocation& loc, const std::string& msg );

    /** Configuration of all generic loggers. */

    static void configure();

    /** Read configuration of loggers from a file.

     \return number of evaluated entries in the file.
     */

    static int readConfig( const char* filename );

    static GenLogger* rootLogger; //!< pointer to root logger for all GenLoggers

    /** Traversing this logger and all output loggers, can be used for DEBUG of Logger */

    void traverse();
};

}
