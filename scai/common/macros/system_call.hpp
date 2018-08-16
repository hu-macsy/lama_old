/**
 * @file system_call.hpp
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
 * @brief Definition of macros for system calls
 * @author Thomas Brandes
 * @date 17.05.2106
 */
#pragma once

#include <cstring>
#include <sstream>
#include <iostream>

#include <scai/common/exception/Exception.hpp>

/**
 * @brief Macro for system call, throws error message in case of error
 */

#define SCAI_SYSTEM_CALL( call, msg )                                                   \
    {                                                                                   \
        int rc = call;                                                                  \
        if ( 0 != rc )                                                                  \
        {                                                                               \
            std::ostringstream errorStr;                                                \
            errorStr << "error in line " << __LINE__;                                   \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Msg  : " << msg << std::endl;                                \
            errorStr << "  Call : " #call << std::endl;                                 \
            errorStr << "  Error: " << strerror( rc );                                  \
            errorStr << ", rc = " << rc << std::endl;                                   \
            scai::common::Exception::addCallStack( errorStr );                          \
            throw scai::common::Exception( errorStr.str() );                            \
        }                                                                               \
    }

/**
 * @brief Macro for system call, prints error message in case of error
 *
 * This macro should be used in destructors where no more exceptions should be throwns
 */

#define SCAI_SYSTEM_CALL_NOTHROW( call, msg )                                           \
    {                                                                                   \
        int rc = call;                                                                  \
        if ( 0 != rc )                                                                  \
        {                                                                               \
            std::ostringstream errorStr;                                                \
            errorStr << "error in line " << __LINE__;                                   \
            errorStr << " of file " << __FILE__ << std::endl;                           \
            errorStr << "  Call : " #call << std::endl;                                 \
            errorStr << "  Error: " << strerror( rc );                                  \
            errorStr << "  Msg  : " << msg << std::endl;                                \
            errorStr << ", rc = " << rc << std::endl;                                   \
            scai::common::Exception::addCallStack( errorStr );                          \
            std::cerr << errorStr.str() << std::endl;                                   \
        }                                                                               \
    }

