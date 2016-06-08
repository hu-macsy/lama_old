/**
 * @file GPIUtils.hpp
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
 * @brief Utility macros for GPI calls
 * @author Lauretta Schubert
 * @date 25.02.2014
 */

#pragma once

#include <scai/dmemo/gpi/GPIException.hpp>

#include <sstream>
#include <cstdio>
#include <iostream>
#include <GASPI.h>

/** This directive might be enabled for debugging  */

#define SCAI_GASPI_CALL( call )                                           \
    {                                                                     \
        gaspi_return_t status = call;                                     \
                                                                          \
        if ( status != GASPI_SUCCESS )                                    \
        {                                                                 \
            std::ostringstream errorStr;                                  \
            errorStr << "GPI error in line " << __LINE__ ;                \
            errorStr << " of file " << __FILE__ << ": ";                  \
            errorStr << #call << "\n";                                    \
            scai::common::Exception::addCallStack( errorStr );            \
            std::fprintf( stderr, "%s\n", errorStr.str().c_str() );       \
            gaspi_printf( "%s\n", errorStr.str().c_str() );               \
            throw GPIException( errorStr.str(), status );                 \
        }                                                                 \
    }

