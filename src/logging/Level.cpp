/**
 * @file Level.cpp
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
 * @brief Level.cpp
 * @author Thomas Brandes
 * @date 01.03.2011
 * @since 1.0.0
 */

#include <logging/Level.hpp>

#include <cstring>

namespace log4lama
{

/** Attention: these strings are of no help as initialization might be done after use. */

static const std::string TRACE_ID = "TRACE";
static const std::string DEBUG_ID = "DEBUG";
static const std::string INFO_ID = "INFO";
static const std::string WARN_ID = "WARN";
static const std::string SERROR_ID = "ERROR";
static const std::string FATAL_ID = "FATAL";
static const std::string OFF_ID = "OFF";
static const std::string UNKOWN_ID = "UNKNOWN";

const char* level2str( const Level level )
{
    switch ( level )
    {

    case TRACE:
        return "TRACE";
    case DEBUG:
        return "DEBUG";
    case INFO:
        return "INFO";
    case WARN:
        return "WARN";
    case SERROR:
        return "ERROR";
    case FATAL:
        return "FATAL";
    case OFF:
        return "OFF";
    default:
        return "UNKNOWN";
    }
}

std::ostream& operator<<( std::ostream& os, const Level& level )
{
    os << level2str( level );
    return os;
}

Level str2level( const std::string& value )
{
    const char* cvalue = value.c_str();

    Level level;

    if ( strcmp( cvalue, "TRACE" ) == 0 )
    {
        level = TRACE;
    }
    else if ( strcmp( cvalue, "DEBUG" ) == 0 )
    {
        level = DEBUG;
    }
    else if ( strcmp( cvalue, "INFO" ) == 0 )
    {
        level = INFO;
    }
    else if ( strcmp( cvalue, "WARN" ) == 0 )
    {
        level = WARN;
    }
    else if ( strcmp( cvalue, "ERROR" ) == 0 )
    {
        level = SERROR;
    }
    else if ( strcmp( cvalue, "FATAL" ) == 0 )
    {
        level = FATAL;
    }
    else if ( strcmp( cvalue, "OFF" ) == 0 )
    {
        level = OFF;
    }
    else
    {
        level = MAXLEVEL; // used as error indicator
    }

    return level;
}

} //namespace log4lama
