/**
 * @file Distributed.cpp
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
 * @brief Distributed.cpp
 * @author Jiri Kraus
 * @date 22.02.2011
 */

// hpp
#include <scai/dmemo/Distributed.hpp>

namespace scai
{

namespace dmemo
{

Distributed::Distributed( DistributionPtr distribution )

    : mDistribution( distribution )
{
    if ( !distribution )
    {
        COMMON_THROWEXCEPTION( "Distributed object must not have NULL distribution" )
    }
}

Distributed::Distributed( const Distributed& other )

    : mDistribution( other.mDistribution )
{
    // copy shared pointer is okay, mDistribution can never be NULL
}

Distributed::~Distributed()
{
}

void Distributed::swap( Distributed& other )
{
    mDistribution.swap( other.mDistribution );
}

void Distributed::setDistributionPtr( DistributionPtr distributionPtr )
{
    mDistribution = distributionPtr;
}

} /* end namespace dmemo */

} /* end namespace scai */
