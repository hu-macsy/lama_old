/**
 * @file BlockDistribution.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Implementation of methods for block distribution class.
 * @author Thomas Brandes
 * @date 18.03.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/distribution/BlockDistribution.hpp>

// local library
#include <scai/lama/matrix/Matrix.hpp>

// std
#include <fstream>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( BlockDistribution::logger, "Distribution.BlockDistribution" )

void BlockDistribution::getRange(
    IndexType& lb,
    IndexType& ub,
    const IndexType n,
    const PartitionId rank,
    const PartitionId size )
{
    SCAI_ASSERT_DEBUG( rank < size, "illegal rank = " << rank << ", size = " << size )
    IndexType blockSize = ( n + size - 1 ) / size;
    lb = rank * blockSize;
    ub = ( rank + 1 ) * blockSize - 1;
    ub = std::min( ub, n - 1 );
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
    getRange( mLB, mUB, globalSize, rank, size );
    SCAI_LOG_INFO( logger,
                   "BlockDistribution of " << getGlobalSize() << " elements" << ", me has " << mLB << " : " << mUB )
}

bool BlockDistribution::isLocal( const IndexType globalIndex ) const
{
    return globalIndex >= mLB && globalIndex <= mUB;
}

PartitionId BlockDistribution::getOwner( const IndexType globalIndex ) const
{
    return globalIndex / mBlockSize;
}

IndexType BlockDistribution::getLocalSize() const
{
    IndexType localSize = 0;

    if( mLB <= mUB )
    {
        localSize = mUB - mLB + 1;
    }

    return localSize;
}

IndexType BlockDistribution::local2global( const IndexType localIndex ) const
{
    return mLB + localIndex;
}

IndexType BlockDistribution::global2local( const IndexType globalIndex ) const
{
    IndexType localIndex = nIndex;

    if( globalIndex >= mLB && globalIndex <= mUB )
    {
        localIndex = globalIndex - mLB;
    }

    return localIndex;
}

void BlockDistribution::computeOwners(
    const std::vector<IndexType>& requiredIndexes,
    std::vector<PartitionId>& owners ) const
{
    owners.clear();
    owners.reserve( requiredIndexes.size() );
    SCAI_LOG_INFO( logger, "compute " << requiredIndexes.size() << " owners for " << *this )

    for( size_t i = 0; i < requiredIndexes.size(); i++ )
    {
        PartitionId owner = getOwner( requiredIndexes[i] );
        owners.push_back( owner );
    }
}

bool BlockDistribution::isEqual( const Distribution& other ) const
{
    if( this == &other )
    {
        return true;
    }

    if( dynamic_cast<const BlockDistribution*>( &other ) )
    {
        return mGlobalSize == other.getGlobalSize();
    }

    return false;
}

void BlockDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "BlockDistribution(gsize=" << mGlobalSize << ",bsize=" << mBlockSize << ")";
}

void BlockDistribution::printDistributionVector( std::string name ) const
{
    PartitionId myRank = mCommunicator->getRank();
    PartitionId parts = mCommunicator->getSize();
    IndexType myLocalSize = getLocalSize();
    std::vector<IndexType> localSizes( parts );
    mCommunicator->gather( &localSizes[0], 1, 0/*MASTER*/, &myLocalSize );

    if( myRank == 0 ) // process 0 is MASTER process
    {
        std::ofstream file;
        file.open( ( name + ".part" ).c_str() );

        // print row - partition mapping
        for( IndexType i = 0; i < parts; ++i )
        {
            for( IndexType j = 0; j < localSizes[i]; j++ )
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

BlockDistribution* BlockDistribution::create(
    const CommunicatorPtr communicator,
    const IndexType globalSize,
    const float )
{
    // weight remains unused
    return new BlockDistribution( globalSize, communicator );
}

BlockDistribution* BlockDistribution::create( const CommunicatorPtr communicator, const Matrix& matrix, const float )
{
    // we only take the size of the matrix

    return new BlockDistribution( matrix.getNumRows(), communicator );
}

/* ---------------------------------------------------------------------------------*/

bool BlockDistribution::initialized = Distribution::registerCreator<BlockDistribution>( "BLOCK" );

} /* end namespace lama */

} /* end namespace scai */
