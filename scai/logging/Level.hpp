/**
 * @file Level.hpp
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
 * @brief Definition of different logging levels.
 * @author Thomas Brandes
 * @date 10.06.2015
 */

#pragma once

// internal scai libraries
#include <scai/common/config.hpp>

// std
#include <ostream>
#include <string>

namespace scai
{

namespace logging
{

/** Own namespace for enum type Level and its values. */

namespace level
{

/** @brief Enumeration type for different logging levels.
 *
 *  The values are ordered. Setting a lower level also enables all
 *  higher levels.
 */
typedef enum
{
    TRACE,   //!< even more detailed than DEBUG
    DEBUG,   //!< designates fine-grained informational events
    INFO,    //!< informational messages highlighting progress
    WARN,    //!< for potentially harmful situations
    SERROR,  //!< for errors that might still allow the application to continue
    FATAL,   //!< severe errors that will presumably lead to aborts
    OFF,     //!< turn logging off
    MAXLEVEL //!< end value for enumeration, also used as unknown level
} Level;

/** Output of level in a stream.
 */
COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& os, const Level& level );

} /* end namespace level */

level::Level str2level( const std::string& value );

/** Translate level to a string.
 *
 *  \return the logging level as a string.
 */

const char* level2str( const level::Level level );

} /* end namespace logging */

} /* end namespace scai */
