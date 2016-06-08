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
 * @endlicense
 *
 * @brief Contains the logLevels for different loggers
 * @author Matthias Makulla
 * @date 06.04.2011
 */

#pragma once

#include <iostream>

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
enum LogLevel
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
    completeInformation = 4
};

/*
 * Output of ScalarType in stream by writing strings instead of numbers
 *
 */

inline std::ostream& operator<<( std::ostream& stream, const LogLevel& object )
{
    switch ( object )
    {
        case noLogging:
            stream << "noLogging";
            break;

        case convergenceHistory:
            stream << "convergenceHistory";
            break;

        case solverInformation:
            stream << "solverInformation";
            break;

        case advancedInformation:
            stream << "advancedInformation";
            break;

        case completeInformation:
            stream << "completeInformation";
            break;

        default:
            stream << "unknown_LogLevel";
    }

    return stream;
}

} /* end namespace LogLevel */

} /* end namespace solver */

} /* end namespace scai */
