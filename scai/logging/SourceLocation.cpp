/**
 * @file SourceLocation.cpp
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
 * @brief Implementation of methods for class SourceLocation.
 * @author Thomas Brandes
 * @date 01.03.2011
 */

// hpp
#include <scai/logging/SourceLocation.hpp>

// std
#include <string.h>

namespace scai
{

namespace logging
{

SourceLocation::SourceLocation( const char* const filename, const char* const funcname, const int line )
    : mFileName( filename ), mFuncName( funcname ), mLine( line )
{
    const char* shortname = strrchr( filename, '/' );

    if ( shortname )
    {
        mFileName = shortname + 1;
    }
}

std::ostream& operator<<( std::ostream& os, const SourceLocation& loc )
{
    os << loc.mFileName << "::" << loc.mLine << ", funct=" << loc.mFuncName;
    return os;
}

} /* end namespace logging */

} /* end namespace scai */
