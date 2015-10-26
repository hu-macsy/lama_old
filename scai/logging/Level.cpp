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

// local library
#include <scai/logging/Level.hpp>

// std
#include <cstring>

namespace scai
{

namespace logging
{

// Attention: Do not use any static (string) variables here 
// as init might be done after first loggings

const char* level2str( const level::Level level )
{
    switch ( level )
    {
        case level::TRACE:
            return "TRACE";

        case level::DEBUG:
            return "DEBUG";

        case level::INFO:
            return "INFO";

        case level::WARN:
            return "WARN";

        case level::SERROR:
            return "ERROR";

        case level::FATAL:
            return "FATAL";

        case level::OFF:
            return "OFF";

        default:
            return "UNKNOWN";
    }
}

level::Level str2level( const std::string& value )
{
    const char* cvalue = value.c_str();
    level::Level level;

    if ( strcmp( cvalue, "TRACE" ) == 0 )
    {
        level = level::TRACE;
    }
    else if ( strcmp( cvalue, "DEBUG" ) == 0 )
    {
        level = level::DEBUG;
    }
    else if ( strcmp( cvalue, "INFO" ) == 0 )
    {
        level = level::INFO;
    }
    else if ( strcmp( cvalue, "WARN" ) == 0 )
    {
        level = level::WARN;
    }
    else if ( strcmp( cvalue, "ERROR" ) == 0 )
    {
        level = level::SERROR;
    }
    else if ( strcmp( cvalue, "FATAL" ) == 0 )
    {
        level = level::FATAL;
    }
    else if ( strcmp( cvalue, "OFF" ) == 0 )
    {
        level = level::OFF;
    }
    else
    {
        level = level::MAXLEVEL; // used as error indicator
    }

    return level;
}

namespace level
{

std::ostream& operator<<( std::ostream& os, const Level& level )
{
    os << scai::logging::level2str( level );
    return os;
}

} /* end namespace level */

} /* end namespace logging */

} /* end namespace scai */
