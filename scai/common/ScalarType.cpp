/**
 * @file ScalarType.cpp
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
 * @brief Implementation of operations on ScalarType.
 * @author Thomas Brandes
 * @date 07.11.2011
 */

// hpp
#include <scai/common/ScalarType.hpp>

#include <cstring>

namespace scai
{

namespace common
{

const char* scalarType2str( const scalar::ScalarType stype )
{
    switch ( stype )
    {
        case scalar::INT:
            return "int";

        case scalar::LONG:
            return "long";

        case scalar::FLOAT:
            return "float";

        case scalar::DOUBLE:
            return "double";

        case scalar::INDEX_TYPE:
            return "IndexType";

        case scalar::LONG_DOUBLE:
            return "LongDouble";

        case scalar::COMPLEX:
            return "ComplexFloat";

        case scalar::DOUBLE_COMPLEX:
            return "ComplexDouble";

        case scalar::LONG_DOUBLE_COMPLEX:
            return "ComplexLongDouble";

        case scalar::INTERNAL:
            return "_Internal";

        case scalar::PATTERN:
            return "_Pattern";

        default:
            return "Unknown";
    }
}

scalar::ScalarType str2ScalarType( const char* str )
{
    for ( int type = scalar::INDEX_TYPE; type < scalar::UNKNOWN; ++type )
    {
        if ( strcmp( scalarType2str( scalar::ScalarType( type ) ), str ) == 0 )
        {
            return scalar::ScalarType( type );
        }
    }

    return scalar::UNKNOWN;
}

bool isComplex( const scalar::ScalarType t )
{
    bool is = false;

    switch ( t )
    {
        case scalar::DOUBLE_COMPLEX:
        case scalar::COMPLEX:
        case scalar::LONG_DOUBLE_COMPLEX :
            is = true;
            break;

        default:
            is = false;
    }

    return is;
}

std::ostream& operator<<( std::ostream& stream, const scalar::ScalarType& object )
{
    stream << scalarType2str( object );
    return stream;
}

} /* end namespace common */

} /* end namespace scai */
