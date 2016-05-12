/**
 * @file ContextType.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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

/** Own struct for the enumeration type ContextType and its values. 
 *
 *  Note: This enumeration type was originally defined in the class Context so
 *        it did not require an own struct. But this is no longer possible as
 *        Context is now a factory where ContextType is the input type to create
 *        a new Context.
 */

struct context
{
    /** Enumeration type for the supported contexts. The type is used to select
     *  the appropriate code that will be used for the computations in the context.
     *
     *  The same context type does not imply that two different contexts can use
     *  the same data. Two CUDA contexts might allocate their own data where data
     *  must be transfered explicitly.
     */
    typedef enum 
    {
        Host,          //!< context for cpu + main memory
        CUDA,          //!< CUDA GPU device
        OpenCL,        //!< OpenCL GPU device, currently not supported
        MIC,           //!< Intel MIC
        UserContext,   //!< can be used for a new derived Context class
        MaxContext     //!< used for dimension of ContextType arrays
    } ContextType;

    /** Enumeration type for access kind, may be read or write */

    typedef enum
    {
        Read, //!<  read access to the array, can be multiple
        Write, //!<  write access to the array, only one at a time
        MaxAccessKind //!<  internal use for dimension of arrays
    } AccessKind;

}; /* end struct context */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const context::ContextType& type );

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const context::AccessKind& kind );

} /* end namespace common */

} /* end namespace scai */
