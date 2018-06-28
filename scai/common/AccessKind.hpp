/**
 * @file AccessKind.hpp
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
 * @brief Definition of enum type for read and write access
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

/** Enumeration type for access kind, may be read or write */

enum class AccessKind
{
    Read, //!<  read access to the array, can be multiple
    Write, //!<  write access to the array, only one at a time
    MaxAccessKind //!<  internal use for dimension of arrays
};

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const AccessKind& kind );

} /* end namespace common */

} /* end namespace scai */
