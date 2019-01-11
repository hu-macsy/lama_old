/**
 * @file SyncKind.hpp
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
 * @brief SyncKind describes if the communication and computation should be done synchronously or asynchronously.
 */
enum class SyncKind
{
    SYNCHRONOUS,  // exchange communication - sync - local + halo computation
    ASYNC_COMM,  // asynchronous execution of communication
    ASYNC_LOCAL  // asynchronous execution of local computation
};

/** This function prints a SyncKind on an output stream.
 *
 *  \param stream  is the reference to the output stream
 *  \param kind    is the enum value that is printed
 */
inline std::ostream& operator<<( std::ostream& stream, const SyncKind& kind )
{
    switch ( kind )
    {
        case SyncKind::SYNCHRONOUS:
        {
            stream << "SYNCHRONOUS";
            break;
        }

        case SyncKind::ASYNC_COMM:
        {
            stream << "ASYN_COMM";
            break;
        }

        case SyncKind::ASYNC_LOCAL:
        {
            stream << "ASYNC_LOCAL";
            break;
        }

        default:
        {
            stream << "<unknown sync kind>";
            break;
        }
    }

    return stream;
}


} /* end namespace lama */

} /* end namespace scai */
