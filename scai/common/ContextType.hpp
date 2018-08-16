/**
 * @file ContextType.hpp
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
 * @brief Definition of enum type for scalar value types typically used in numerical applications.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <ostream>

namespace scai
{

namespace common
{

/** Enumeration type for the supported contexts. The type is used to select
 *  the appropriate code that will be used for the computations in the context.
 *
 *  The same context type does not imply that two different contexts can use
 *  the same data. Two CUDA contexts might allocate their own data where data
 *  must be transfered explicitly.
 */
enum class ContextType
{
    Host,          //!< context for cpu + main memory
    CUDA,          //!< CUDA GPU device
    OpenCL,        //!< OpenCL GPU device, currently not supported
    UserContext,   //!< can be used for a new derived Context class
    MaxContext     //!< used for dimension of ContextType arrays
};

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const ContextType& type );

COMMON_DLL_IMPORTEXPORT const char* contextType2str( const ContextType type );

COMMON_DLL_IMPORTEXPORT ContextType str2ContextType( const char* str );

} /* end namespace common */

} /* end namespace scai */
