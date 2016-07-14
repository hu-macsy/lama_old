/**
 * @file GenBlockDistribution.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief GenBlockDistribution.cpp
 * @author Thomas Brandes
 * @date 18.03.2011
 */

// hpp
#include <scai/dmemo/GenBlockDistribution.hpp>

// local library
#include <scai/dmemo/Distributed.hpp>

// internal scai libraries
#include <scai/common/unique_ptr.hpp>

// std
#include <fstream>

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

void GenBlockDistribution::setOffsets(
    const IndexType rank,
    const IndexType numPartitions,
    const IndexType localSizes[] )
{
    mOffsets.reset( new IndexType[numPartitions] );
    IndexType sumSizes = 0;

    for ( PartitionId p = 0; p < numPartitions; p++ )
    {
        sumSizes += localSizes[p];
        mOffsets[p] = sumSizes;
        SCAI_LOG_TRACE( logger,
                        "Partition " << p << ": local size =  " << localSizes[p] << ", offset = " << mOffsets[p] )
    }

    SCAI_ASSERT_EQUAL( sumSizes, getGlobalSize(), "sum over local sizes must be global size" )
    mUB = mOffsets[rank] - 1;
    mLB = mOffsets[rank] - localSizes[rank];
}

void GenBlockDistribution::setOffsets( const IndexType rank, const IndexType numPartitions, const IndexType mySize )
{
    common::scoped_array<IndexType> localSizes( new IndexType[numPartitions] );
    // rank 0 is root
    mCommunicator->gather( localSizes.get(), 1, 0, &mySize );
    mCommunicator->bcast( localSizes.get(), numPartitions, 0 );
    SCAI_ASSERT_EQ_DEBUG( localSizes[rank], mySize, "wrongly gathered values" )
    setOffsets( rank, numPartitions, localSizes.get() );
}

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const std::vector<IndexType>& localSizes,
    const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator ), mOffsets( new IndexType[mCommunicator->getSize()] )
{
    PartitionId size = mCommunicator->getSize();
    PartitionId rank = mCommunicator->getRank();
    SCAI_LOG_INFO( logger, "GenBlockDistribution of " << getGlobalSize() << " elements" )
    SCAI_ASSERT_EQ_ERROR( size, static_cast<PartitionId>( localSizes.size() ), "size mismatch" )
    setOffsets( rank, size, &localSizes[0] );
    SCAI_LOG_INFO( logger, *this << ": constructed by local sizes" )
}

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const IndexType firstGlobalIdx,
    const IndexType lastGlobalIdx,
    const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator ), mOffsets( new IndexType[mCommunicator->getSize()] )
{
    PartitionId size = mCommunicator->getSize();
    PartitionId rank = mCommunicator->getRank();
    setOffsets( rank, size, lastGlobalIdx - firstGlobalIdx + 1 );
    SCAI_ASSERT_EQUAL( mLB, firstGlobalIdx, "serious mismatch in index range" )
    SCAI_ASSERT_EQUAL( mUB, lastGlobalIdx, "serious mismatch in index range" )
    SCAI_LOG_INFO( logger, *this << ": constructed by local range " << firstGlobalIdx << ":" << lastGlobalIdx )
}

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const IndexType localSize,
    const CommunicatorPtr communicator )
    : Distribution( globalSize, communicator ), mOffsets( new IndexType[mCommunicator->getSize()] )
{
    int size = mCommunicator->getSize();
    int rank = mCommunicator->getRank();
    setOffsets( rank, size, localSize );
}

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const float weight,
    const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator ), mOffsets( new IndexType[mCommunicator->getSize()] )
{
    int size = mCommunicator->getSize();
    int rank = mCommunicator->getRank();
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
    mOffsets.reset( new IndexType[size] );
    float sumWeight = 0.0f;

    for ( PartitionId p = 0; p < size; p++ )
    {
        sumWeight += allWeights[p];
        mOffsets[p] = static_cast<IndexType>( sumWeight / totalWeight * globalSize + 0.5 );
    }

    mLB = 0;
    mUB = mOffsets[rank] - 1;

    if ( rank > 0 )
    {
        mLB = mOffsets[rank - 1];
    }

    SCAI_LOG_INFO( logger, *this << " constructed by weight factors" )
}

bool GenBlockDistribution::isLocal( const IndexType globalIndex ) const
{
    return globalIndex >= mLB && globalIndex <= mUB;
}

void GenBlockDistribution::getLocalRange( IndexType& lb, IndexType& ub ) const
{
    lb = mLB;
    ub = mUB;
}

PartitionId GenBlockDistribution::getOwner( const IndexType globalIndex ) const
{
    int first = 0;
    int last = mCommunicator->getSize() - 1;

    if ( globalIndex < 0 )
    {
        return -1; // out of range
    }

    if ( globalIndex >= mOffsets[last] )
    {
        return -1; // out of range
    }

    // binary search in the array mOffsets

    while ( first < last )
    {
        int mid = ( first + last + 1 ) / 2;

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

    return first;
}

IndexType GenBlockDistribution::getLocalSize() const
{
    IndexType localSize = 0;

    if ( mLB <= mUB )
    {
        localSize = mUB - mLB + 1;
    }

    return localSize;
}

IndexType GenBlockDistribution::local2global( const IndexType localIndex ) const
{
    return mLB + localIndex;
}

/* ---------------------------------------------------------------------- */

IndexType GenBlockDistribution::global2local( const IndexType globalIndex ) const
{
    IndexType localIndex = nIndex; // default value if globalIndex is not local

    if ( globalIndex >= mLB && globalIndex <= mUB )
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

void GenBlockDistribution::computeOwners1(
    const std::vector<IndexType>& requiredIndexes,
    std::vector<PartitionId>& owners ) const
{
    owners.clear();
    owners.reserve( requiredIndexes.size() );
    SCAI_LOG_INFO( logger, "compute " << requiredIndexes.size() << " owners for " << *this )

    for ( unsigned int i = 0; i < requiredIndexes.size(); ++i )
    {
        IndexType requiredIndex = requiredIndexes[i];
        PartitionId owner = getOwner( requiredIndex );
        owners.push_back( owner );
        SCAI_LOG_TRACE( logger, "owner of required index = " << requiredIndex << " is " << owner )
    }
}

bool GenBlockDistribution::isEqual( const Distribution& other ) const
{
    if ( this == &other )
    {
        return true; // pointer equality, is always okay
    }

    if ( *mCommunicator != other.getCommunicator() )
    {
        return false;
    }

    const GenBlockDistribution* otherBlock = dynamic_cast<const GenBlockDistribution*>( &other );

    if ( !otherBlock )
    {
        return false;
    }

    bool equal = true;

    for ( PartitionId p = 0; p < mCommunicator->getSize(); ++p )
    {
        if ( mOffsets[p] != otherBlock->mOffsets[p] )
        {
            equal = false;
            break;
        }
    }

    return equal;
}

void GenBlockDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "GenBlockDistribution( gsize = " << mGlobalSize << ", local = " << mLB << ":" << mUB << ")";
}

void GenBlockDistribution::printDistributionVector( std::string name ) const
{
    IndexType myRank = mCommunicator->getRank();
    IndexType parts = mCommunicator->getSize();
    IndexType myLocalSize = getLocalSize();
    std::vector<IndexType> localSizes( parts );
    mCommunicator->gather( &localSizes[0], 1, MASTER, &myLocalSize );

    if ( myRank == MASTER ) // process 0 is MASTER process
    {
        std::ofstream file;
        file.open( ( name + ".part" ).c_str() );

        // print row - partition mapping
        for ( IndexType i = 0; i < parts; ++i )
        {
            for ( IndexType j = 0; j < localSizes[i]; j++ )
            {
                file << i << std::endl;
            }
        }

        file.close();
    }
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

const char GenBlockDistribution::theCreateValue[] = "GEN_BLOCK";

std::string GenBlockDistribution::createValue()
{
    return theCreateValue;
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
