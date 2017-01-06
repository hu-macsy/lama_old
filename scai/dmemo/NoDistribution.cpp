/**
 * @file NoDistribution.cpp
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
 * @brief Implementation of methods for class NoDistribution.
 * @author Thomas Brandes
 * @date 14.03.2011
 */

// hpp
#include <scai/dmemo/NoDistribution.hpp>

// std
#include <fstream>

#define MASTER 0

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( NoDistribution::logger, "Distribution.NoDistribution" )

/* ---------------------------------------------------------------------- */

NoDistribution::NoDistribution( const IndexType globalSize )
    : Distribution( globalSize )
{
}

/* ---------------------------------------------------------------------- */

NoDistribution::~NoDistribution()
{
}
 
/* ---------------------------------------------------------------------- */

bool NoDistribution::isLocal( const IndexType /* index */ ) const
{
    return true;
}

/* ---------------------------------------------------------------------- */

IndexType NoDistribution::getLocalSize() const
{
    return mGlobalSize;
}

/* ---------------------------------------------------------------------- */

IndexType NoDistribution::local2global( const IndexType localIndex ) const
{
    return localIndex;
}

/* ---------------------------------------------------------------------- */

IndexType NoDistribution::global2local( const IndexType globalIndex ) const
{
    return globalIndex;
}

/* ---------------------------------------------------------------------- */

IndexType NoDistribution::getBlockDistributionSize() const
{
    return mGlobalSize;
}

/* ---------------------------------------------------------------------- */

bool NoDistribution::isEqual( const Distribution& other ) const
{
    bool isSame = false;

    bool proven = proveEquality( isSame, other );

    if ( proven )
    {
        return isSame;
    }

    return false;
}

/* ---------------------------------------------------------------------- */

void NoDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "NoDistribution( size = " << mGlobalSize << " )";
}

/* ---------------------------------------------------------------------- */

void NoDistribution::computeOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes ) const
{
    PartitionId root = 0;
    owners.init( root, indexes.size() ); 
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

Distribution* NoDistribution::create( const DistributionArguments arg )
{
    // Note: weight argument is not used here
    //       same is true for matrix, commonunicationPtr
    return new NoDistribution( arg.globalSize );
}

} /* end namespace dmemo */

} /* end namespace scai */
