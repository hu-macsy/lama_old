/**
 * @file ContextType.cpp
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
 * @brief Implementation of operations on ScalarType.
 * @author Jiri Kraus
 * @date 07.11.2011
 */

// hpp
#include <scai/common/ContextType.hpp>

#include <cstring>

namespace scai
{

namespace common
{

const char* contextType2str( const ContextType type )
{
    switch ( type )
    {
        case ContextType::Host :
            return "Host";

        case ContextType::CUDA :
            return "CUDA";

        case ContextType::OpenCL :
            return "OpenCL";

        case ContextType::UserContext :
            return "UserContext";

        default:
            return "UnknownContext";
    }
}

/* -----------------------------------------------------------------------------*/

ContextType str2ContextType( const char* str )
{
    for ( int type = 0; type < static_cast<int>( ContextType::MaxContext ); ++type )
    {
        ContextType ctx = ContextType( type );

        if ( strcmp( contextType2str( ctx ), str ) == 0 )
        {
            return ctx;
        }
    }

    return ContextType::MaxContext;
}

/* -----------------------------------------------------------------------------*/

std::ostream& operator<<( std::ostream& stream, const ContextType& type )
{
    stream << contextType2str( type );
    return stream;
}

} /* end namespace context */

} /* end namespace scai */
