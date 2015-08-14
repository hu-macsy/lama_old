/**
 * @file Level.cpp
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
 * @brief Implementation of conversion routines between level enum and strings.
 * @author Thomas Brandes
 * @date 01.03.2011
 */

#include <scai/logging/Level.hpp>

#include <cstring>

namespace logging
{

// Attention: Do not use any static (string) variables here 
// as init might be done after first loggings

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

} //namespace logging
