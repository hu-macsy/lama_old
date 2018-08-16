/**
 * @file LoggerWriteBehaviour.hpp
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
 * @brief Contains an enumeration which identifies whether a logger
 *        logs its messages to the console only or to a file and the console
 * @author Matthias Makulla
 * @date 06.04.2011
 */
#pragma once

namespace scai
{

namespace solver
{

/**
 * @brief Contains an enumeration which identifies whether a logger
 *        logs its messages to the console only or to a file and the console
 */
enum class LoggerWriteBehaviour
{
    /**
     * @brief Log messages will be written to standard out only.
     */
    toConsoleOnly,

    /**
     * @brief Log messages will be written to the log file only.
     */
    toFileOnly,

    /**
     * @brief Log messages will be written to the console and the logfile.
     */
    toFileAndConsole
};

/*
 * Output of LoggerWriteBehaviour in stream by writing strings instead of numbers
 *
 */
inline std::ostream& operator<<( std::ostream& stream, const LoggerWriteBehaviour& object )
{
    switch ( object )
    {
        case LoggerWriteBehaviour::toConsoleOnly:

            stream << "toConsoleOnly";
            break;

        case LoggerWriteBehaviour::toFileOnly:

            stream << "toFileOnly";
            break;
 
        case LoggerWriteBehaviour::toFileAndConsole:

            stream << "toFileOnly";
            break;

        default:
            stream << "UNKNOWN";
    }

    return stream;
}

} /* end namespace solver */

} /* end namespace scai */
