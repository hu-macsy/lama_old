/**
 * @file GenBlockDistribution.cpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief GenBlockDistribution.cpp
 * @author Thomas Brandes
 * @date 18.03.2011
 * $Id$
 */

// hpp
#include <lama/distribution/GenBlockDistribution.hpp>

#include <fstream>

namespace lama
{

LAMA_LOG_DEF_LOGGER( GenBlockDistribution::logger, "Distribution.GenBlockDistribution" )

GenBlockDistribution::~GenBlockDistribution()
{
    LAMA_LOG_INFO( logger, "~GenBlockDistribution" )
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
        LAMA_LOG_TRACE( logger,
                        "Partition " << p << ": local size =  " << localSizes[p] << ", offset = " << mOffsets[p] )
    }

    LAMA_ASSERT_EQUAL_ERROR( sumSizes, getGlobalSize() )

    mUB = mOffsets[rank] - 1;
    mLB = mOffsets[rank] - localSizes[rank];
}

void GenBlockDistribution::setOffsets( const IndexType rank, const IndexType numPartitions, const IndexType mySize )
{
    boost::scoped_array<IndexType> localSizes( new IndexType[numPartitions] );

    const PartitionId root = 0;
    mCommunicator->gather( localSizes.get(), 1, root, &mySize );
    mCommunicator->bcast( localSizes.get(), numPartitions, root );

    LAMA_ASSERT_EQUAL_DEBUG( localSizes[ rank ], mySize )

    setOffsets( rank, numPartitions, localSizes.get() );
}

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const std::vector<IndexType>& localSizes,
    const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator )
{
    PartitionId size = mCommunicator->getSize();
    PartitionId rank = mCommunicator->getRank();

    LAMA_LOG_INFO( logger, "GenBlockDistribution of " << getGlobalSize() << " elements" )

    LAMA_ASSERT_EQUAL_ERROR( size, static_cast<PartitionId>( localSizes.size() ) )

    setOffsets( rank, size, &localSizes[0] );

    LAMA_LOG_INFO( logger, *this << ": constructed by local sizes" )
}

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const IndexType firstGlobalIdx,
    const IndexType lastGlobalIdx,
    const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator )
{
    PartitionId size = mCommunicator->getSize();
    PartitionId rank = mCommunicator->getRank();

    setOffsets( rank, size, lastGlobalIdx - firstGlobalIdx + 1 );

    LAMA_ASSERT_EQUAL_ERROR( mLB, firstGlobalIdx )
    LAMA_ASSERT_EQUAL_ERROR( mUB, lastGlobalIdx )

    LAMA_LOG_INFO( logger, *this << ": constructed by local range " << firstGlobalIdx << ":" << lastGlobalIdx )
}

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const IndexType localSize,
    const CommunicatorPtr communicator )
    : Distribution( globalSize, communicator )
{
    int size = mCommunicator->getSize();
    int rank = mCommunicator->getRank();

    setOffsets( rank, size, localSize );
}

GenBlockDistribution::GenBlockDistribution(
    const IndexType globalSize,
    const float weight,
    const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator )
{
    int size = mCommunicator->getSize();
    int rank = mCommunicator->getRank();

    LAMA_LOG_DEBUG( logger,
                    "GenBlockDistribution of " << getGlobalSize() << " elements" << ", my weight = " << weight )

    std::vector<float> allWeights;
    communicator->gather( allWeights, weight );
    float totalWeight = 0;

    for ( PartitionId p = 0; p < size; p++ )
    {
        if ( allWeights[p] < /*=*/0 )
        {
            LAMA_THROWEXCEPTION(
                "Weight of partition " << p << " = " << allWeights[p] << " illegal, must be positive" );
        }

        totalWeight += allWeights[p];
    }

    LAMA_LOG_INFO( logger,
                   "GenBlockDistribution of " << getGlobalSize() << " elements" << ", total weight = " << totalWeight )
    mOffsets.reset( new IndexType[size] );
    float sumWeight = 0.0;

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

    LAMA_LOG_INFO( logger, *this << " constructed by weight factors" )
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

IndexType GenBlockDistribution::global2local( const IndexType globalIndex ) const
{
    IndexType localIndex = nIndex; // default value if globalIndex is not local

    if ( globalIndex >= mLB && globalIndex <= mUB )
    {
        localIndex = globalIndex - mLB;
    }

    return localIndex;
}

void GenBlockDistribution::computeOwners(
    const std::vector<IndexType>& requiredIndexes,
    std::vector<PartitionId>& owners ) const
{
    owners.clear();
    owners.reserve( requiredIndexes.size() );
    LAMA_LOG_INFO( logger, "compute " << requiredIndexes.size() << " owners for " << *this )

    for ( unsigned int i = 0; i < requiredIndexes.size(); ++i )
    {
        IndexType requiredIndex = requiredIndexes[i];
        PartitionId owner = getOwner( requiredIndex );
        owners.push_back( owner );
        LAMA_LOG_TRACE( logger, "owner of required index = " << requiredIndex << " is " << owner )
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
    mCommunicator->gather( &localSizes[0], 1, 0/*MASTER*/, &myLocalSize );

    if ( myRank == 0 ) // process 0 is MASTER process
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

}
