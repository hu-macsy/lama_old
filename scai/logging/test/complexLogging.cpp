/**
 * @file scai/logging/test/complexLogging.cpp
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
 * @brief simple executable that uses the hierarchical logging of SCAI logging
 * @author Jan Ecker
 * @date 03.09.2015
 */

#include <scai/logging.hpp>

SCAI_LOG_DEF_LOGGER( logger_c1,   "Class1" )
SCAI_LOG_DEF_LOGGER( logger_c1m1, "Class1.method1" )
SCAI_LOG_DEF_LOGGER( logger_c1m2, "Class1.method2" )
SCAI_LOG_DEF_LOGGER( logger_c1m2r1, "Class1.method2.region1" )
SCAI_LOG_DEF_LOGGER( logger_c2,   "Class2" )
SCAI_LOG_DEF_LOGGER( logger_c2m1, "Class2.method1" )
SCAI_LOG_DEF_LOGGER( logger_c2m2, "Class2.method2" )

int main()
{
    SCAI_LOG_DEBUG( logger_c1,   "message class1" )
    SCAI_LOG_TRACE( logger_c1m1,  "message class1 method1" )
    SCAI_LOG_INFO( logger_c1m2,  "message class1 method2" )
    SCAI_LOG_INFO( logger_c1m2r1, "message class1 method2 region1" )
    SCAI_LOG_WARN( logger_c2,    "message class2" )
    SCAI_LOG_INFO( logger_c2m1,  "message class2 method1" )
    SCAI_LOG_TRACE( logger_c2m2, "message class2 method2" )
}
