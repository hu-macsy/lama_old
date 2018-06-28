/**
 * @file LogLevels.cpp
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
 * @brief Simple example that shows using the logging library.
 *        /
 * @author Thomas Brandes
 * @date 21.08.2015
 */

#include <scai/logging.hpp>

SCAI_LOG_DEF_LOGGER( myLogger, "Demo" )

int main( int, char** )
{
    // macro to give the current thread a name that appears in further logs
    SCAI_LOG_THREAD( "main" )
    SCAI_LOG_INFO( myLogger, "a message about progress in the program" )
    SCAI_LOG_DEBUG( myLogger, "a message useful to find bugs in the program" )
    SCAI_LOG_TRACE( myLogger, "a message with very detailled info, usually not compiled" )
    SCAI_LOG_WARN( myLogger, "a message with a warning, but execution is still possible" )
    SCAI_LOG_ERROR( myLogger, "a message for an error, error handling will be invoked" )
    SCAI_LOG_FATAL( myLogger, "a message for a fatal error, execution will stop" )
}
