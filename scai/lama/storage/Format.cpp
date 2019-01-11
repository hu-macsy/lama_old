/**
 * @file Format.cpp
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
 * @brief Implementation of functions for enum class Format
 * @author Thomas Brandes
 * @date 31.10.2017
 */

// local library
#include <scai/lama/storage/Format.hpp>


#include <cstring>

namespace scai
{

namespace lama
{

const char* format2Str( const Format storageFormat )
{
    switch ( storageFormat )
    {
        case Format::CSR:
            return "CSR";
            break;

        case Format::ELL:
            return "ELL";
            break;

        case Format::DIA:
            return "DIA";
            break;

        case Format::JDS:
            return "JDS";
            break;

        case Format::COO:
            return "COO";
            break;

        case Format::DENSE:
            return "Dense";
            break;

        case Format::STENCIL:
            return "Stencil";
            break;

        case Format::UNDEFINED:
            return "UNDEFINED";
            break;
    }

    return "UNDEFINED";
}

Format str2Format( const char* str )
{
    for ( int format = 0; format < static_cast<int>( Format::UNDEFINED ); ++format )
    {
        if ( strcmp( format2Str( Format( format ) ), str ) == 0 )
        {
            return Format( format );
        }
    }

    return Format::UNDEFINED;
}

std::ostream& operator<<( std::ostream& stream, const Format& storageFormat )
{
    stream << scai::lama::format2Str( storageFormat );
    return stream;
}


} /* end namespace lama */

} /* end namespace scai */
