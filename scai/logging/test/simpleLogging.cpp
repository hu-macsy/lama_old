/**
 * @file scai/logging/test/simpleLogging.cpp
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
 * @brief simple executable that uses each log level of SCAI logging
 * @author Jan Ecker
 * @date 03.09.2015
 */

#include <scai/logging.hpp>

SCAI_LOG_DEF_LOGGER( logger, "Test" )

int main()
{
    SCAI_LOG_TRACE( logger, "trace message" )
    SCAI_LOG_DEBUG( logger, "debug message" )
    SCAI_LOG_INFO( logger, "info message" )
    SCAI_LOG_WARN( logger, "warn message" )
    SCAI_LOG_ERROR( logger, "error message" )
    SCAI_LOG_FATAL( logger, "fatal message" )
}
