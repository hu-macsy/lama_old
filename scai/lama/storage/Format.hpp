/**
 * @file Format.hpp
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
 * @brief Definition of enum class for supported matrix storage formats.
 * @author Thomas Brandes
 * @date 31.10.2017
 */

#pragma once

#include <scai/common/config.hpp>

#include <ostream>

namespace scai
{

namespace lama
{

/** Enumeration type for different matrix storage formats.
 *
 *  Note: operator<< is implemented for this type and should be adapted in case of modifications.
 */
enum class Format
{
    DENSE,    //!< Dense, all elements are stored
    CSR,      //!< Compressed Sparse Row
    ELL,      //!< ELLPack
    DIA,      //!< Diagonal
    JDS,      //!< Jagged Diagonal Storage
    COO,      //!< Coordinate list
    STENCIL,  //!< stencil pattern
    UNDEFINED //!< Default value
};

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const Format& storageFormat );

COMMON_DLL_IMPORTEXPORT const char* format2Str( const Format storageFormat );

COMMON_DLL_IMPORTEXPORT Format str2Format( const char* str );

} /* end namespace lama */

} /* end namespace scai */
