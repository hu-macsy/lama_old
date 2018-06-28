/**
 * @file Level.cpp
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
