/**
 * @file VectorKind.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Enumeration class for different kinds of vector, e.g. SPARSE, DENSE
 * @author Thomas Brandes, Jiri Kraus
 * @date 31.10.2017
 */
#pragma once

#include <scai/common/config.hpp>

#include <ostream>

namespace scai
{

namespace lama
{

/**
 * @brief Enum class for the different kind of vector
 */
enum class VectorKind
{
    DENSE,      //!< vector format for a dense vector
    SPARSE,     //!< vector format for a sparse vector
    JOINED,     //!< vector format for a joined vector
    UNDEFINED   //!< for convenience, always the last entry, stands also for number of entries
}; 

COMMON_DLL_IMPORTEXPORT const char* vectorKind2Str( const VectorKind vectorKind );

COMMON_DLL_IMPORTEXPORT VectorKind str2VectorKind( const char* str );

/** @brief Output operator<< for VectorKind prints meaningful names instead of int values */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const VectorKind& kind );

} /* end namespace lama */

} /* end namespace scai */
