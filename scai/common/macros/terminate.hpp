/**
 * @file terminate.hpp
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
 * @brief Definition of macro that terminates application
 * @author Thomas Brandes
 * @date 04.01.2017
 */
#pragma once

#include <scai/common/exception/Exception.hpp>
#include <iostream>

#define SCAI_TERMINATE( msg )                                                            \
    {                                                                                    \
        std::ostringstream errorStr;                                                     \
        errorStr<<"Terminate in line " << __LINE__ << " of file " << __FILE__ << "\n";   \
        errorStr << "    Message: " << msg << "\n";                                      \
        scai::common::Exception::addCallStack( errorStr );                               \
        std::cerr << errorStr.str()  << std::endl;                                       \
        std::terminate();                                                                \
   }
