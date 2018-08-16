/**
 * @file VectorKind.cpp
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
 * @brief Implementations of functions for enum class VectorKind
 * @author Thomas Brandes, Jiri Kraus
 * @date 31.10.2017
 */

#include <scai/lama/VectorKind.hpp>

#include <cstring>

namespace scai
{

namespace lama
{

/* ---------------------------------------------------------------------------------- */

const char* vectorKind2Str( const VectorKind vectorKind )
{
    switch ( vectorKind )
    {
        case VectorKind::DENSE:
            return "DENSE";
            break;

        case VectorKind::SPARSE:
            return "SPARSE";
            break;

        case VectorKind::JOINED:
            return "JOINED";
            break;

        case VectorKind::UNDEFINED:
            return "Undefined";
            break;
    }

    return "<illegal_vector_kind>";
}

VectorKind str2VectorKind( const char* str )

{
    for ( int kind = 0; kind < static_cast<int>( VectorKind::UNDEFINED ); ++kind )
    {
        if ( strcmp( vectorKind2Str( VectorKind( kind ) ), str ) == 0 )
        {
            return VectorKind( kind );
        }
    }

    return VectorKind::UNDEFINED;
}

/* ---------------------------------------------------------------------------------------*/
/*    VectorKind opertor<<                                                                */
/* ---------------------------------------------------------------------------------------*/

std::ostream& operator<<( std::ostream& stream, const VectorKind& kind )
{
    stream << vectorKind2Str( kind );
    return stream;
}

/* ---------------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
