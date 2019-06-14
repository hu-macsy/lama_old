/**
 * @file BlockDistribution.cpp
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
 * @brief Implementation of methods for block distribution class.
 * @author Thomas Brandes
 * @date 18.03.2011
 */

// hpp
#include <scai/dmemo/BlockDistribution.hpp>

// std
#include <fstream>

#define MASTER 0

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( BlockDistribution::logger, "Distribution.BlockDistribution" )

void BlockDistribution::getLocalRange(
    IndexType& lb,
    IndexType& ub,
    const IndexType n,
    const PartitionId rank,
    const PartitionId size )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( rank, size, "illegal rank specified" )
    IndexType blockSize = ( n + size - 1 ) / size;
    lb = rank * blockSize;
    ub = ( rank + 1 ) * blockSize;
    ub = std::min( ub, n );
}

BlockDistribution::~BlockDistribution()
{
    SCAI_LOG_INFO( logger, "~BlockDistribution" )
}

BlockDistribution::BlockDistribution( const IndexType globalSize, const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator )
{
    PartitionId size = mCommunicator->getSize();
    PartitionId rank = mCommunicator->getRank();
    SCAI_LOG_DEBUG( logger, "BlockDistribution of " << getGlobalSize() << " elements" )
    mBlockSize = ( globalSize + size - 1 ) / size;
    getLocalRange( mLB, mUB, globalSize, rank, size );
    SCAI_LOG_INFO( logger,
                   "BlockDistribution of " << getGlobalSize() << " elements" << ", me has " << mLB << " : " << mUB )
}

bool BlockDistribution::isLocal( const IndexType globalIndex ) const
{
    return globalIndex >= mLB && globalIndex < mUB;
}

/* ---------------------------------------------------------------------- */

PartitionId BlockDistribution::getOwner( const IndexType globalIndex ) const
{
    return globalIndex / mBlockSize;
}

/* ---------------------------------------------------------------------- */

IndexType BlockDistribution::getLocalSize() const
{
    IndexType localSize = 0;

    if ( mLB < mUB )
    {
        localSize = mUB - mLB;
    }

    return localSize;
}

/* ---------------------------------------------------------------------- */

IndexType BlockDistribution::getMaxLocalSize() const
{
    return mBlockSize;
}

/* ---------------------------------------------------------------------- */

IndexType BlockDistribution::getBlockDistributionSize() const
{
    return getLocalSize();
}

/* ---------------------------------------------------------------------- */

IndexType BlockDistribution::local2Global( const IndexType localIndex ) const
{
    return mLB + localIndex;
}

/* ---------------------------------------------------------------------- */

IndexType BlockDistribution::global2Local( const IndexType globalIndex ) const
{
    IndexType localIndex = invalidIndex;

    if ( globalIndex >= mLB && globalIndex < mUB )
    {
        localIndex = globalIndex - mLB;
    }

    return localIndex;
}

/* ---------------------------------------------------------------------- */

void BlockDistribution::computeOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes ) const
{
    ContextPtr ctx = Context::getHostPtr();    // currently only available @ Host

    const IndexType n = indexes.size();

    ReadAccess<IndexType> rIndexes( indexes, ctx );
    WriteOnlyAccess<PartitionId> wOwners( owners, ctx, n );

    // ToDo: call a kernel and allow arbitrary context

    for ( IndexType i = 0; i < n; i++ )
    {
        wOwners[i] = rIndexes[i] / mBlockSize;   // same as getOwner
    }
}

/* ---------------------------------------------------------------------- */

void BlockDistribution::getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const
{
    const IndexType nLocal  = getLocalSize();

    SCAI_LOG_INFO( logger, getCommunicator() << ": getOwnedIndexes, have " << nLocal << " of " << mGlobalSize )

    IndexType i = mLB;

    for ( IndexType& globalIndex : hostWriteOnlyAccess( myGlobalIndexes, nLocal ) )
    {
        globalIndex = i++;
    }
}

/* ---------------------------------------------------------------------- */

bool BlockDistribution::isEqual( const Distribution& other ) const
{
    bool isSame = false;

    bool proven = proveEquality( isSame, other );

    if ( proven )
    {
        return isSame;
    }

    if ( strcmp( other.getKind(), getKind() ) == 0 )
    {
        // both BlockDistribution, we know already that global size and communicator are equal

        isSame = true;
    }
    else if ( strcmp( other.getKind(), "GEN_BLOCK" ) == 0 )
    {
        // general block distribution can check if it is just a block distribution

        isSame = other.isEqual( *this );
    }

    return isSame;
}

/* ---------------------------------------------------------------------- */

void BlockDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "BlockDistribution( comm = " << *mCommunicator << ", block = " << mBlockSize
           << ", size = " << mLB << ":" << mUB << " of " << mGlobalSize <<  " )";
}

/* ---------------------------------------------------------------------- */

bool BlockDistribution::hasAnyAddressing() const
{
    return true;
}

void BlockDistribution::enableAnyAddressing() const
{
    // done, here we have closed formulas
}

IndexType BlockDistribution::getAnyLocalSize( const PartitionId rank ) const
{
    PartitionId nP = mCommunicator->getSize();

    SCAI_ASSERT_VALID_INDEX_DEBUG( rank, nP, "illegal rank specified" )

    IndexType lb;
    IndexType ub;

    getLocalRange( lb, ub, getGlobalSize(), rank, nP );

    // Attention:  lb = 14 > ub = 12 might be possible

    if ( lb < ub )
    {
        return ub - lb;
    }
    else
    {
        return 0;
    }
}

PartitionId BlockDistribution::getAnyOwner( const IndexType globalIndex ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( globalIndex, mGlobalSize, "global index out of range" )
    return globalIndex / mBlockSize;
}

IndexType BlockDistribution::getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( globalIndex, mGlobalSize, "global index out of range" )

    if ( owner == invalidIndex )
    {
        return globalIndex % mBlockSize;
    }
    else
    {
        // avoid expensive ownership

        return globalIndex - owner * mBlockSize;
    }
}

IndexType BlockDistribution::getAnyGlobalIndex( const IndexType localIndex, const PartitionId owner ) const
{
    IndexType localSize = getAnyLocalSize( owner );

    SCAI_ASSERT_VALID_INDEX_DEBUG( localIndex, localSize, "local index out of range for owner = " << owner )

    return localIndex + owner * mBlockSize;
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string BlockDistribution::createValue()
{
    return getId();
}

DistributionPtr BlockDistribution::create( const DistributionArguments arg )
{
    SCAI_LOG_INFO( logger, "create" )
    // Note: weight argument is not used here
    return blockDistribution( arg.globalSize, arg.communicator );
}

} /* end namespace dmemo */

} /* end namespace scai */
