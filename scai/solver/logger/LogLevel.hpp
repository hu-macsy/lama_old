/**
 * @file LogLevel.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Contains the logLevels for different loggers
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

#include <iostream>
#include <cstring>

namespace scai
{

namespace solver
{

/**
 * @brief Contains the logLevels for different loggers.
 */
namespace LogLevel
{

/**
 * @brief The different logLevels for loggers.
 */
typedef enum
{

    /**
     * @brief At this log level no information will be logged, used
     *        NullLogger instead.
     */
    noLogging = 0,

    /**
     * @brief Logs the convergence history of the solver - iterations
     *        and residuals.
     */
    convergenceHistory = 1,

    /**
     * @brief More information about the solver will be logged.
     */
    solverInformation = 2,

    /**
     * @brief Advanced solver information like residual requests
     *        and stopping criteria checks will belogged.
     */
    advancedInformation = 3,

    /**
     * @brief Logs every log message of the solver.
     */
    completeInformation = 4,

    /**
     * @brief unspecified
     */
    UNKNOWN = 5

} LogLevel;


} /* end namespace LogLevel */

inline const char* logLevel2str( const LogLevel::LogLevel level )
{
    switch ( level )
    {
        case LogLevel::noLogging:
            return "noLogging";

        case LogLevel::convergenceHistory:
            return "convergenceHistory";

        case LogLevel::solverInformation:
            return "solverInformation";

        case LogLevel::advancedInformation:
            return "advancedInformation";

        case LogLevel::completeInformation:
            return "completeInformation";

        default:
            return "Unknown";
    }
}

namespace LogLevel
{
/*
 * Output of ScalarType in stream by writing strings instead of numbers
 *
 */

inline std::ostream& operator<<( std::ostream& stream, const LogLevel& object )
{
    stream << logLevel2str( object );
    return stream;
}

}

inline LogLevel::LogLevel str2LogLevel( const char* str )
{
    for ( int level = LogLevel::noLogging; level < LogLevel::UNKNOWN; ++level )
    {
        if ( strcmp( logLevel2str( LogLevel::LogLevel( level ) ), str ) == 0 )
        {
            return LogLevel::LogLevel( level );
        }
    }

    return LogLevel::UNKNOWN;
}

} /* end namespace solver */

} /* end namespace scai */
