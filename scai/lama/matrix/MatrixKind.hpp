/**
 * @file MatrixKind.hpp
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
 * @brief Abstract base class for all matrices supported by LAMA.
 * @author Jiri Kraus, Thomas Brandes
 * @date 22.02.2011
 */
#pragma once

#include <scai/common/config.hpp>

#include <ostream>

namespace scai
{

namespace lama
{

/**
 * @brief Enum class for the different kind of matrix
 */
enum class MatrixKind
{
    DENSE,      //!< matrix format for a dense vector
    SPARSE,     //!< matrix format for a sparse matrix
    JOINED,     //!< matrix format for a joined matrix
    UNDEFINED   //!< for convenience, always the last entry, stands also for number of entries
};

/** This function prints a MatrixKind on an output stream.
 *
 *  \param stream   is the reference to the output stream
 *  \param kind      is the enum value that is printed
 */
inline std::ostream& operator<<( std::ostream& stream, const MatrixKind& kind )
{
    switch ( kind )
    {
        case scai::lama::MatrixKind::DENSE:
        {
            stream << "DENSE";
            break;
        }

        case scai::lama::MatrixKind::SPARSE:
        {
            stream << "SPARSE";
            break;
        }

        default:
        {
            stream << "<unknown matrix kind>";
            break;
        }
    }

    return stream;
}

} /* end namespace lama */

} /* end namespace scai */
