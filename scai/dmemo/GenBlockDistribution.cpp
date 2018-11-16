/**
 * @file GenBlockDistribution.cpp
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
 * @brief GenBlockDistribution.cpp
 * @author Thomas Brandes
 * @date 18.03.2011
 */

// hpp
#include <scai/dmemo/GenBlockDistribution.hpp>

// local library
#include <scai/dmemo/Distributed.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

// std
#include <fstream>
#include <memory>

#define MASTER 0

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GenBlockDistribution::logger, "Distribution.GenBlockDistribution" )

GenBlockDistribution::~GenBlockDistribution()
{
    SCAI_LOG_INFO( logger, "~GenBlockDistribution" )
}

/* ---------------------------------------------------------------------- */

void GenBlockDistribution::setOffsets(
    const PartitionId rank,
    const PartitionId numPartitions,
    const IndexType localSizes[] )
{
    mOffsets.reset( new IndexType[numPartitions + 1] );
    IndexType sumSizes = 0;

    mOffsets[0] = 0;   // just filling element for more convenient use of offset array

    for ( PartitionId p = 0; p < numPartitions; p++ )
    {
        sumSizes += localSizes[p];
        mOffsets[p + 1] = sumSizes;
        SCAI_LOG_TRACE( logger,
                        "Partition " << p << ": local size =  " << localSizes[p] << ", offset = " << mOffsets[p + 1] )
    }

    SCAI_ASSERT_EQUAL( sumSizes, getGlobalSize(), "sum over local sizes must be global size" )

    // mUB is first element not in the local range

    mUB = mOffsets[rank + 1];
    mLB = mOffsets[rank];
}

/* ---------------------------------------------------------------------- */

void GenBlockDistribution::setOffsets( const PartitionId rank, const PartitionId numPartitions, const IndexType mySize )
{
    std::unique_ptr<IndexType[]> localSizes( new IndexType[numPartitions] );
    // rank 0 is root
    mCommunicator->gather( localSizes.get(), 1, 0, &mySize );
    mCommunicator->bcast( localSizes.get(), numPartitions, 0 );
    SCAI_ASSERT_EQ_DEBUG( localSizes[rank], mySize, "wrongly gathered values" )
    setOffsets( rank, numPartitions, localSizes.get() );
}

/* ---------------------------------------------------------------------- */

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const std::vector<IndexType>& localSizes,
    const CommunicatorPtr communicator ) :

    Distribution( globalSize, communicator )
{
    PartitionId size = mCommunicator->getSize();
    PartitionId rank = mCommunicator->getRank();
    SCAI_LOG_INFO( logger, "GenBlockDistribution of " << getGlobalSize() << " elements" )
    SCAI_ASSERT_EQ_ERROR( size, static_cast<PartitionId>( localSizes.size() ), "size mismatch" )
    setOffsets( rank, size, &localSizes[0] );
    SCAI_LOG_INFO( logger, *this << ": constructed by local sizes" )
}

/* ---------------------------------------------------------------------- */

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const IndexType firstGlobalIdx,
    const IndexType lastGlobalIdx,
    bool,
    const CommunicatorPtr communicator ) :

    Distribution( globalSize, communicator )

{
    SCAI_ASSERT_LE_ERROR( firstGlobalIdx, lastGlobalIdx, "illegal local range" )

    PartitionId size = mCommunicator->getSize();
    PartitionId rank = mCommunicator->getRank();
    setOffsets( rank, size, lastGlobalIdx - firstGlobalIdx );
    SCAI_ASSERT_EQUAL( mLB, firstGlobalIdx, "serious mismatch in index range" )
    SCAI_ASSERT_EQUAL( mUB, lastGlobalIdx, "serious mismatch in index range" )
    SCAI_LOG_INFO( logger, *this << ": constructed by local range " << firstGlobalIdx << ":" << lastGlobalIdx )
}

/* ---------------------------------------------------------------------- */

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const IndexType localSize,
    const CommunicatorPtr communicator ) :

    Distribution( globalSize, communicator )

{
    int size = mCommunicator->getSize();
    int rank = mCommunicator->getRank();
    setOffsets( rank, size, localSize );
}

/* ---------------------------------------------------------------------- */

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const float weight,
    const CommunicatorPtr communicator ) :

    Distribution( globalSize, communicator )

{
    PartitionId size = mCommunicator->getSize();
    PartitionId rank = mCommunicator->getRank();
    SCAI_LOG_DEBUG( logger, "GenBlockDistribution of " << getGlobalSize() << " elements" << ", my weight = " << weight )
    std::vector<float> allWeights( size );
    communicator->allgather( &allWeights[0], 1, &weight );
    float totalWeight = 0;

    for ( PartitionId p = 0; p < size; p++ )
    {
        if ( allWeights[p] < /*=*/0 )
        {
            COMMON_THROWEXCEPTION( "Weight of partition " << p << " = " << allWeights[p] << " illegal, must be positive" );
        }

        totalWeight += allWeights[p];
    }

    SCAI_LOG_INFO( logger,
                   "GenBlockDistribution of " << getGlobalSize() << " elements" << ", total weight = " << totalWeight )

    mOffsets.reset( new IndexType[size + 1] );
    float sumWeight = 0.0f;

    mOffsets[0] = 0;

    for ( PartitionId p = 0; p < size; p++ )
    {
        sumWeight += allWeights[p];
        mOffsets[p + 1] = static_cast<IndexType>( sumWeight / totalWeight * globalSize + 0.5 );
    }

    mLB = mOffsets[rank];
    mUB = mOffsets[rank + 1];

    SCAI_LOG_INFO( logger, *this << " constructed by weight factors" )
}

bool GenBlockDistribution::isLocal( const IndexType globalIndex ) const
{
    return globalIndex >= mLB && globalIndex < mUB;
}

void GenBlockDistribution::getLocalRange( IndexType& lb, IndexType& ub ) const
{
    // keep in mind that lb == ub implies zero range, ub < lb can never happen

    lb = mLB;
    ub = mUB;
}

/* ---------------------------------------------------------------------- */

PartitionId GenBlockDistribution::getOwner( const IndexType globalIndex ) const
{
    // owner of an index can be computed by each processor without communication

    PartitionId first = 1;
    PartitionId last  = mCommunicator->getSize();

    if ( ! common::Utils::validIndex( globalIndex, mOffsets[last] ) )
    {
        return invalidPartition; // out of range
    }

    // binary search in the array mOffsets

    while ( first < last )
    {
        PartitionId mid = ( first + last + 1 ) / 2;

        if ( globalIndex < mOffsets[mid - 1] )
        {
            last = mid - 1;
        }
        else if ( globalIndex >= mOffsets[mid] )
        {
            first = mid + 1;
        }
        else
        {
            first = mid;
            last = mid;
        }
    }

    return first - 1;
}

/* ---------------------------------------------------------------------- */

IndexType GenBlockDistribution::getLocalSize() const
{
    IndexType localSize = 0;

    if ( mLB < mUB )
    {
        localSize = mUB - mLB;
    }

    return localSize;
}

/* ---------------------------------------------------------------------- */

IndexType GenBlockDistribution::local2global( const IndexType localIndex ) const
{
    return mLB + localIndex;
}

/* ---------------------------------------------------------------------- */

IndexType GenBlockDistribution::global2local( const IndexType globalIndex ) const
{
    IndexType localIndex = invalidIndex;   // default value if globalIndex is not local

    if ( globalIndex >= mLB && globalIndex < mUB )
    {
        localIndex = globalIndex - mLB;
    }

    return localIndex;
}

/* ---------------------------------------------------------------------- */

void GenBlockDistribution::computeOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes ) const
{
    ContextPtr ctx = Context::getHostPtr();    // currently only available @ Host

    const IndexType n = indexes.size();

    ReadAccess<IndexType> rIndexes( indexes, ctx );
    WriteOnlyAccess<PartitionId> wOwners( owners, ctx, n );

    // ToDo: call a kernel and allow arbitrary context

    for ( IndexType i = 0; i < n; i++ )
    {
        wOwners[i] = getOwner( rIndexes[i] );
    }
}

/* ---------------------------------------------------------------------- */

void GenBlockDistribution::getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const
{
    const IndexType nLocal = getLocalSize();

    SCAI_LOG_INFO( logger, getCommunicator() << ": getOwnedIndexes, have " << nLocal << " of " << mGlobalSize )

    const IndexType one = 1;   // avoids cast in argument list of setSequence

    utilskernel::HArrayUtils::setSequence( myGlobalIndexes, mLB, one, nLocal );
}

/* ---------------------------------------------------------------------- */

IndexType GenBlockDistribution::getBlockDistributionSize() const
{
    return getLocalSize();
}

/* ---------------------------------------------------------------------- */

bool GenBlockDistribution::hasAnyAddressing() const
{
    return true;
}

void GenBlockDistribution::enableAnyAddressing() const
{
    // done, here we have closed formulas
}

IndexType GenBlockDistribution::getAnyLocalSize( const PartitionId rank ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( rank, mCommunicator->getSize(), "illegal rank" )

    return mOffsets[rank + 1 ] - mOffsets[rank];
}

PartitionId GenBlockDistribution::getAnyOwner( const IndexType globalIndex ) const
{
    return getOwner( globalIndex );
}

IndexType GenBlockDistribution::getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( globalIndex, mGlobalSize, "illegal index for distribution" )

    // here the owner is very helpful

    return globalIndex - mOffsets[ owner ];
}

IndexType GenBlockDistribution::getAnyGlobalIndex( const IndexType localIndex, const PartitionId rank ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( localIndex, getAnyLocalSize( rank ), "Illegal local index for rank = " << rank )

    return localIndex + mOffsets[ rank ];
}

/* ---------------------------------------------------------------------- */

bool GenBlockDistribution::isEqual( const Distribution& other ) const
{
    bool isSame = false;

    bool proven = proveEquality( isSame, other );

    if ( proven )
    {
        return isSame;
    }

    if ( other.getKind() != getKind() )
    {
        return false;
    }

    const GenBlockDistribution& otherBlock = reinterpret_cast<const GenBlockDistribution&>( other );

    bool equal = true;

    for ( PartitionId p = 0; p < mCommunicator->getSize() + 1; ++p )
    {
        if ( mOffsets[p] != otherBlock.mOffsets[p] )
        {
            equal = false;
            break;
        }
    }

    return equal;
}

/* ---------------------------------------------------------------------- */

void GenBlockDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "GenBlockDistribution( gsize = " << mGlobalSize << ", local = " << mLB << ":" << mUB << ")";
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string GenBlockDistribution::createValue()
{
    return getId();
}

Distribution* GenBlockDistribution::create( const DistributionArguments arg )
{
    if ( arg.matrix != NULL )
    {
        SCAI_LOG_WARN( logger, "matrix argument ignored to create GEN_BLOCK distribution" )
    }

    return new GenBlockDistribution( arg.globalSize, arg.weight, arg.communicator );
}


} /* end namespace dmemo */

} /* end namespace scai */
