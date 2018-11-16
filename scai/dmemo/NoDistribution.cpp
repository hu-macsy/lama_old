/**
 * @file NoDistribution.cpp
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

NoDistribution::NoDistribution( const IndexType globalSize, CommunicatorPtr comm ) : 

    Distribution( globalSize, comm )
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

bool NoDistribution::hasAnyAddressing() const
{
    return true;
}

void NoDistribution::enableAnyAddressing() const
{
}

IndexType NoDistribution::getAnyLocalSize( const PartitionId ) const
{
    return getGlobalSize();
}

PartitionId NoDistribution::getAnyOwner( const IndexType ) const
{
    return 0;
}

IndexType NoDistribution::getAnyLocalIndex( const IndexType globalIndex, const PartitionId ) const
{
    return globalIndex;
}

IndexType NoDistribution::getAnyGlobalIndex( const IndexType localIndex, const PartitionId ) const
{
    return localIndex;
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
    stream << "NoDistribution( size = " << mGlobalSize << ", comm = " << getCommunicator() << " )";
}

/* ---------------------------------------------------------------------- */

void NoDistribution::computeOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes ) const
{
    PartitionId root = 0;
    owners.setSameValue( indexes.size(), root );
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

const char* NoDistribution::getId()
{
    static const char id[] = "NO";
    return id;
}

} /* end namespace dmemo */

} /* end namespace scai */
