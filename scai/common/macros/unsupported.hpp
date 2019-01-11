/**
 * @file unsupported.hpp
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
 * @brief Macro to give user info about an unsupported feature.
 * @author Thomas Brandes
 * @date 10.08.2015
 */

#pragma once

#include <scai/common/exception/UnsupportedException.hpp>

#include <iostream>
#include <sstream>

/** This macro should be used to give hints about unsupported features.
 *
 *  It should only be used in cases where less performant solutions are still
 *  available.
 *
 *  By setting the environment variable COMMON_UNSUPPORTED to WARN only warnings
 *  will be given. For IGNORE no message is given at all. Otherwise an exception
 *  is thrown.
 */

#define SCAI_UNSUPPORTED( msg )                                                \
    {                                                                              \
        if ( scai::common::UnsupportedException::getUnsupportedSetting() !=        \
                scai::common::UnsupportedException::UNSUPPORTED_IGNORE )           \
        {                                                                          \
            std::ostringstream errorStr;                                           \
            errorStr << "Unsupported at line ";                                    \
            errorStr << __LINE__ << " of file " << __FILE__ << "\n";               \
            errorStr << "    Message: " << msg << std::endl;                       \
            errorStr << "Use environment variable SCAI_UNSUPPORTED";               \
            errorStr << " (WARN or IGNORE) to get rid of this message";            \
            errorStr << std::endl;                                                 \
            if ( scai::common::UnsupportedException::getUnsupportedSetting() ==    \
                    scai::common::UnsupportedException::UNSUPPORTED_ERROR )        \
            {                                                                      \
                throw scai::common::UnsupportedException( errorStr.str() );        \
            }                                                                      \
            else                                                                   \
            {                                                                      \
                std::cout << errorStr.str() << std::endl;                          \
            }                                                                      \
        }                                                                          \
    }

