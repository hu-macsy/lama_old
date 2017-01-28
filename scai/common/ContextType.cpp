/**
 * @file ContextType.cpp
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
 * @author Jiri Kraus
 * @date 07.11.2011
 */

// hpp
#include <scai/common/ContextType.hpp>

namespace scai
{

namespace common
{

const char* contextType2str( const context::ContextType& type )
{
    switch ( type )
    {
        case context::Host :
            return "Host";

        case context::CUDA :
            return "CUDA";

        case context::MIC :
            return "MIC";

        case context::OpenCL :
            return "OpenCL";

        case context::UserContext :
            return "UserContext";

        default:
            return "UnknownContext";
    }
}

std::ostream& operator<<( std::ostream& stream, const context::ContextType& type )
{
    stream << contextType2str( type );
    return stream;
}

/* -----------------------------------------------------------------------------*/

std::ostream& operator<<( std::ostream& stream, const context::AccessKind& kind )
{
    switch ( kind )
    {
        case context::Write :
            stream << "Write";
            break;

        case context::Read :
            stream << "Read";
            break;

        default:
            stream << "AccessKind_" << ( int ) kind;
    }

    return stream;
}

} /* end namespace context */

} /* end namespace scai */
