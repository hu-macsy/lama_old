/**
 * @file SingleDistribution.cpp
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
 * @brief Implementation of methods for single distribution class.
 * @author Thomas Brandes
 * @date 30.01.2017
 */

// hpp
#include <scai/dmemo/SingleDistribution.hpp>

// common

#include <scai/common/macros/assert.hpp>

// std
#include <fstream>

#define MASTER 0

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( SingleDistribution::logger, "Distribution.SingleDistribution" )

SingleDistribution::~SingleDistribution()
{
    SCAI_LOG_DEBUG( logger, "~SingleDistribution" )
}

SingleDistribution::SingleDistribution( const IndexType globalSize, const CommunicatorPtr communicator, const PartitionId owner ) :

    Distribution( globalSize, communicator ),
    mOwner( owner )

{
    SCAI_ASSERT_VALID_INDEX_ERROR( owner, communicator->getSize(), "owner is not a valid rank of any processor" )

    // owner must have the same value on all processors, otherwise this would be rather strange

    SCAI_LOG_INFO( logger,
                   "SingleDistribution of " << getGlobalSize() << " elements" << ", owner " << owner )
}

bool SingleDistribution::isLocal( const IndexType ) const
{
    return mCommunicator->getRank() == mOwner;
}

/* ---------------------------------------------------------------------- */

PartitionId SingleDistribution::findOwner( const IndexType ) const
{
    return mOwner;
}

/* ---------------------------------------------------------------------- */

IndexType SingleDistribution::getLocalSize() const
{
    if ( mOwner == mCommunicator->getRank() )
    {
        return mGlobalSize;
    }
    else
    {
        return 0;
    }
}

/* ---------------------------------------------------------------------- */

IndexType SingleDistribution::getBlockDistributionSize() const
{
    return getLocalSize();
}

/* ---------------------------------------------------------------------- */

IndexType SingleDistribution::getMaxLocalSize() const
{
    return mGlobalSize;
}

/* ---------------------------------------------------------------------- */

IndexType SingleDistribution::local2global( const IndexType localIndex ) const
{
    return localIndex;
}

/* ---------------------------------------------------------------------- */

IndexType SingleDistribution::global2local( const IndexType globalIndex ) const
{
    IndexType localIndex = invalidIndex;

    if ( mOwner == mCommunicator->getRank() )
    {
        localIndex = globalIndex;
    }

    return localIndex;
}

/* ---------------------------------------------------------------------- */

void SingleDistribution::computeOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes ) const
{
    ContextPtr ctx = Context::getHostPtr();    // currently only available @ Host

    const IndexType n = indexes.size();

    WriteOnlyAccess<PartitionId> wOwners( owners, ctx, n );

    // ToDo: call a kernel and allow arbitrary context

    for ( IndexType i = 0; i < n; i++ )
    {
        wOwners[i] = mOwner;
    }
}

/* ---------------------------------------------------------------------- */

void SingleDistribution::getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const
{
    const IndexType nLocal  = getLocalSize();

    SCAI_LOG_INFO( logger, getCommunicator() << ": getOwnedIndexes, have " << nLocal << " of " << mGlobalSize )

    WriteOnlyAccess<IndexType> wGlobalIndexes( myGlobalIndexes, nLocal );

    for ( IndexType i = 0; i < nLocal; ++i )
    {
        wGlobalIndexes[ i ] = i;
    }
}

/* ---------------------------------------------------------------------- */

bool SingleDistribution::hasAnyAddressing() const
{
    return true;
}

void SingleDistribution::enableAnyAddressing() const
{
    // nothing to do, there are simple formulas to compute it
}

IndexType SingleDistribution::getAnyLocalSize( const PartitionId rank ) const
{
    return rank == mOwner ? getGlobalSize() : 0;
}

PartitionId SingleDistribution::getAnyOwner( const IndexType globalIndex ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( globalIndex, getGlobalSize(), "illegal globalIndex" )
    return mOwner;
}

IndexType SingleDistribution::getAnyLocalIndex( const IndexType globalIndex, const PartitionId ) const
{
    return globalIndex;
}

IndexType SingleDistribution::getAnyGlobalIndex( const IndexType localIndex, const PartitionId rank ) const
{
    SCAI_ASSERT_EQ_DEBUG( rank, mOwner, "only elements @ partition " << mOwner << ", is single owner" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( localIndex, getGlobalSize(), "illegal local index" )
    return localIndex;
}

/* ---------------------------------------------------------------------- */

bool SingleDistribution::isEqual( const Distribution& other ) const
{
    bool isSame = false;

    bool proven = proveEquality( isSame, other );

    if ( proven )
    {
        return isSame;
    }

    if ( other.getKind() == getKind() )
    {
        const SingleDistribution& sOther = reinterpret_cast<const SingleDistribution&>( other );
        isSame = sOther.mOwner == mOwner;
    }

    // we know already that global size and communicator are equal

    return isSame;
}

/* ---------------------------------------------------------------------- */

void SingleDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "SingleDistribution( size = " << mGlobalSize << ", comm = " << *mCommunicator << ", owner = " << mOwner << " )";
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string SingleDistribution::createValue()
{
    return getId();
}

Distribution* SingleDistribution::create( const DistributionArguments arg )
{
    SCAI_LOG_INFO( logger, "create" )

    // by default we create one where processor 0 is owner

    return new SingleDistribution( arg.globalSize, arg.communicator, 0 );
}

} /* end namespace dmemo */

} /* end namespace scai */
