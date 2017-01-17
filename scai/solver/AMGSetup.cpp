/**
 * @file AMGSetup.cpp
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
 * @brief AMGSetup.cpp
 * @author Jiri Kraus
 * @date 28.10.2011
 */

// hpp
#include <scai/solver/AMGSetup.hpp>

namespace scai
{

namespace solver
{

AMGSetup::AMGSetup()
    : mHostOnlyLevel( std::numeric_limits<IndexType>::max() ), mHostOnlyVars( 0 ), mReplicatedLevel(
        std::numeric_limits<IndexType>::max() )
{
}

AMGSetup::~AMGSetup()
{
}

void AMGSetup::setHostOnlyLevel( IndexType hostOnlyLevel )
{
    mHostOnlyLevel = hostOnlyLevel;
}

void AMGSetup::setHostOnlyVars( IndexType hostOnlyVars )
{
    mHostOnlyLevel = hostOnlyVars;
}

void AMGSetup::setReplicatedLevel( IndexType replicatedLevel )
{
    mReplicatedLevel = replicatedLevel;
}

void AMGSetup::writeAt( std::ostream& stream ) const
{
    stream << "AMGSetup( ... )";
}

} /* end namespace solver */

} /* end namespace scai */
