/**
 * @file LogLevel.hpp
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
 * @brief Contains the logLevels for different loggers
 * @author Matthias Makulla
 * @date 06.04.2011
 * @since 1.0.0
 */

#pragma once

#pragma offload_attribute (push, target(mic))
#include <iostream>
#pragma offload_attribute (pop)

namespace scai
{

namespace lama
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

} /* end namespace lama */

} /* end namespace scai */
