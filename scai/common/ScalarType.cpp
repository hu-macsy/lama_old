/**
 * @file ScalarType.cpp
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
 * @brief Implementation of operations on ScalarType.
 * @author Thomas Brandes
 * @date 07.11.2011
 */

// hpp
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/ScalarTypeHelper.hpp>
#include <scai/common/macros/assert.hpp>

#include <cstring>

namespace scai
{

namespace common
{

const char* scalarType2str( const ScalarType stype )
{
    switch ( stype )
    {
        case ScalarType::INT:
            return "int";

        case ScalarType::LONG:
            return "long";

        case ScalarType::FLOAT:
            return "float";

        case ScalarType::DOUBLE:
            return "double";

        case ScalarType::LONG_DOUBLE:
            return "LongDouble";

        case ScalarType::COMPLEX:
            return "ComplexFloat";

        case ScalarType::DOUBLE_COMPLEX:
            return "ComplexDouble";

        case ScalarType::LONG_DOUBLE_COMPLEX:
            return "ComplexLongDouble";

        case ScalarType::CHAR:
            return "char";

        case ScalarType::UNSIGNED_INT:
            return "uint";

        case ScalarType::UNSIGNED_LONG:
            return "ulong";

        case ScalarType::INDEX_TYPE:
            return "IndexType";

        case ScalarType::INTERNAL:
            return "_Internal";

        case ScalarType::PATTERN:
            return "_Pattern";

        default:
            return "Unknown";
    }
}

/*************************************************************************/

ScalarType str2ScalarType( const char* str )
{
    for ( int stype = 0; stype < static_cast<int>( ScalarType::UNKNOWN ); ++stype )
    {
        if ( strcmp( scalarType2str( ScalarType( stype ) ), str ) == 0 )
        {
            return ScalarType( stype );
        }
    }

    return ScalarType::UNKNOWN;
}

/*************************************************************************/

std::ostream& operator<<( std::ostream& stream, const ScalarType& object )
{
    stream << scalarType2str( object );
    return stream;
}

/*************************************************************************/

bool isComplex( const ScalarType stype )
{
    return common::mepr::ScalarTypeHelper<SCAI_ALL_TYPES_LIST>::isComplex( stype );
}

/*************************************************************************/

bool isNumeric( const ScalarType stype )
{
    return common::mepr::ScalarTypeHelper<SCAI_ALL_TYPES_LIST>::isNumeric( stype );
}

/*************************************************************************/


int precision( const ScalarType stype )
{
    // call with INTERNAL not allowed, replace it at call site with precision( TypeTraits<ValueType>::sid )

    SCAI_ASSERT( stype != ScalarType::INTERNAL, "precision of INTERNAL unknown" )

    // loop over all supported types and query its precision

    return common::mepr::ScalarTypeHelper<SCAI_ALL_TYPES_LIST>::precision( stype );
}

/*************************************************************************/

size_t typeSize( const ScalarType stype )
{
    // call with INTERNAL not allowed, replace it at cal size with typeSize( TypeTraits<ValueType>::sid )

    SCAI_ASSERT( stype != ScalarType::INTERNAL, "typeSize of INTERNAL unknown" )

    if ( stype == ScalarType::PATTERN )
    {
        return 0;      // allowed, stands for a dummy type
    }

    if ( stype == ScalarType::INDEX_TYPE )
    {
        return sizeof( IndexType );   // must be handled separately here
    }

    // for all types use TypeTraits as they contain already this information

    return common::mepr::ScalarTypeHelper<SCAI_ALL_TYPES_LIST>::sizeOf( stype );
}

/*************************************************************************/


} /* end namespace common */

} /* end namespace scai */
