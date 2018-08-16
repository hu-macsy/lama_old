/**
 * @file ContextAccess.cpp
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
 * @brief Implemenations for class ContextAccess.
 * @author Thomas Brandes
 * @date 14.07.2011
 */

// hpp
#include <scai/hmemo/ContextAccess.hpp>

namespace scai
{

namespace hmemo
{

SCAI_LOG_DEF_LOGGER( ContextAccess::logger, "ContextAccess" )

ContextAccess::ContextAccess( ContextPtr context, const char* file, int line )
    : mContext( *context ), mReleased( false ), mFile( file ), mLine( line )
{
    SCAI_LOG_INFO( logger, *this << " enabled" )
    mContext.enable( mFile, mLine );
}

void ContextAccess::release()
{
    if ( mReleased )
    {
        return;
    }

    SCAI_LOG_INFO( logger, *this << " released" )
    mContext.disable( mFile, mLine );
    mReleased = false;
}

ContextAccess::~ContextAccess()
{
    release();
}

void ContextAccess::writeAt( std::ostream& stream ) const
{
    stream << "Access of " << mContext;
    stream << " at " << mFile << "( line = " << mLine << " )";
}

} /* end namespace hmemo */

} /* end namespace scai */
