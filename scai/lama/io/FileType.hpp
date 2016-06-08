/**
 * @file FileType.hpp
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
 * @brief FileType.hpp
 * @author Kai Buschulte
 * @date 12.05.2010
 */

#pragma once

// internal scai libraries

#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/ScalarTypeHelper.hpp>

namespace scai
{

namespace lama
{

/** Define an own namespace for enumeration types. */

namespace File
{
/**
 * @brief Defines the supported file types
 */
enum FileType
{
    /**
     * @brief the SAMG format
     */
    SAMG_FORMAT,

    /**
     * @brief the Matrix Market Format
     *        (see http://math.nist.gov/matrixMarket for details).
     */
    MATRIX_MARKET,

    /**
     * @brief unspecified, used internally
     */
    DEFAULT
};

/*
 * Output of ScalarType in stream by writing strings instead of numbers
 */

static inline std::ostream& operator<<( std::ostream& stream, const FileType& object )
{
    switch ( object )
    {
        case SAMG_FORMAT:
            stream << "SAMG";
            break;

        case MATRIX_MARKET:
            stream << "MATRIX_MARKET";
            break;

        default:
            stream << "<unknown_file_type>";
    }

    return stream;
}

}  // namespace File

} /* end namespace lama */

} /* end namespace scai */
