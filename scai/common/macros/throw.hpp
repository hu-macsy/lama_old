/**
 * @file throw.hpp
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
 * @brief Definition of macros that throw exceptions
 * @author Thomas Brandes
 * @date 11.11.2015
 */
#pragma once

#include <scai/common/exception/Exception.hpp>

/**
 * @brief The macro SCAI_THROWEXCEPTION throws an exception that contains
 *        source code file and line as well as call stack in its message.
 *
 * @param[in] ExceptionClass must be same or derived class from scai::common::Exception
 * @param[in] msg   message to indicate reason for the exception
 * @throws    ExceptionClass (derived from scai::common:Exception, derived from std::exception)
 */

#define SCAI_THROWEXCEPTION( ExceptionClass, msg )                             \
    {                                                                              \
        std::ostringstream errorStr;                                               \
        errorStr<<"Exception in line "<<__LINE__<<" of file "<<__FILE__<<"\n";     \
        errorStr<<"    Message: "<<msg<<"\n";                                      \
        scai::common::Exception::addCallStack( errorStr );                         \
        throw ExceptionClass( errorStr.str() );                                    \
    }

/** COMMON_THROWEXCEPTION just throws a simple exception */

#define COMMON_THROWEXCEPTION( msg )   \
    SCAI_THROWEXCEPTION( scai::common::Exception, msg )
