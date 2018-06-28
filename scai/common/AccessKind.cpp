/**
 * @file AccessKind.cpp
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
 * @brief Implementation of operations on enum type AccessKind.
 * @author Jiri Kraus
 * @date 07.11.2011
 */

#include <scai/common/AccessKind.hpp>

#include <cstring>

namespace scai
{

namespace common
{

/* -----------------------------------------------------------------------------*/

std::ostream& operator<<( std::ostream& stream, const AccessKind& kind )
{
    switch ( kind )
    {
        case AccessKind::Write :
            stream << "Write";
            break;

        case AccessKind::Read :
            stream << "Read";
            break;

        default:
            stream << "AccessKind_" << static_cast<int>( kind );
    }

    return stream;
}

/* -----------------------------------------------------------------------------*/

} /* end namespace common */

} /* end namespace scai */
