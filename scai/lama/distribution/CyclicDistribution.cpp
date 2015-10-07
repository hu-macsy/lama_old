/**
 * @file CyclicDistribution.cpp
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
 * @brief CyclicDistribution.cpp
 * @author schubert
 * @date 20.05.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/distribution/CyclicDistribution.hpp>

// internal scai library
#include <scai/common/Constants.hpp>

// std
#include <fstream>

#define MASTER Constants<IndexType>::zero

namespace scai
{

using common::Constants;

namespace lama
{

SCAI_LOG_DEF_LOGGER( CyclicDistribution::logger, "Distribution.CyclicDistribution" )

CyclicDistribution::~CyclicDistribution()
{
    SCAI_LOG_INFO( logger, "~CyclicDistribution" )
}

CyclicDistribution::CyclicDistribution(
    const IndexType globalSize,
    const IndexType chunkSize,
    const CommunicatorPtr communicator )

    : Distribution( globalSize, communicator ), mChunkSize( chunkSize )
{
    SCAI_LOG_INFO( logger, "CyclicDistribution of " << mGlobalSize << " elements " << " and chunk size " << mChunkSize )
}

PartitionId CyclicDistribution::getOwner( const IndexType globalIndex ) const
{
    IndexType size = mCommunicator->getSize();
    IndexType globalChunkIndex = globalIndex / mChunkSize;

    return globalChunkIndex % size;
}

bool CyclicDistribution::isLocal( const IndexType globalIndex ) const
{
    IndexType rank = mCommunicator->getRank();

    if( getOwner( globalIndex ) == rank )
    {
        SCAI_LOG_TRACE( logger, "global index " << globalIndex << " is local on partition " << rank )
        return true;
    }
    else
    {
        SCAI_LOG_TRACE( logger, "global index " << globalIndex << " is not local on partition " << rank )
        return false;
    }
}

IndexType CyclicDistribution::getLocalSize() const
{
    const PartitionId rank = mCommunicator->getRank();

    const IndexType elements = getPartitionSize( rank );

    SCAI_LOG_TRACE( logger, *mCommunicator << ": local size = " << elements << " elements" )

    return elements;
}

void CyclicDistribution::getChunkInfo( IndexType& localChunks, IndexType& extra, const PartitionId rank ) const
{
    const PartitionId size = mCommunicator->getSize();

    const IndexType chunks = mGlobalSize / mChunkSize;

    // mGlobalSize % mChunkSize elements remain, are counted later

    localChunks = chunks / size;

    IndexType remainChunks = chunks % size;
    extra = 0;

    if( rank < remainChunks )
    {
        localChunks++;
    }
    else if( rank == remainChunks )
    {
        extra = mGlobalSize % mChunkSize;
    }

    SCAI_LOG_TRACE( logger,
                    "Partition " << rank << " of " << size << ": " << localChunks << " of " << chunks << ", extra = " << extra )
}

IndexType CyclicDistribution::getNumChunks( const PartitionId rank ) const
{
    IndexType localChunks = 0;
    IndexType extra = 0;

    getChunkInfo( localChunks, extra, rank );

    if( extra )
    {
        ++localChunks; // count also the non-full chunk
    }

    return localChunks;
}

IndexType CyclicDistribution::getNumLocalChunks() const
{
    return getNumChunks( mCommunicator->getRank() );
}

IndexType CyclicDistribution::getNumTotalChunks() const
{
    IndexType numChunks = mGlobalSize / mChunkSize;

    if( mGlobalSize % mChunkSize != 0 )
    {
        ++numChunks; // count also the non-full chunk
    }

    return numChunks;
}

IndexType CyclicDistribution::getPartitionSize( const PartitionId partition ) const
{
    IndexType localChunks = 0;
    IndexType extra = 0;

    getChunkInfo( localChunks, extra, partition );

    IndexType elements = localChunks * mChunkSize + extra;

    return elements;
}

IndexType CyclicDistribution::local2global( const IndexType localIndex ) const
{
    IndexType size = mCommunicator->getSize();
    IndexType rank = mCommunicator->getRank();

    IndexType localChunk = localIndex / mChunkSize;
    IndexType localOffset = localIndex % mChunkSize;
    IndexType globalChunk = localChunk * size + rank;
    IndexType globalIndex = globalChunk * mChunkSize + localOffset;
    SCAI_LOG_TRACE( logger,
                    "local Index " << localIndex << " is with chunkSize " << mChunkSize << " and " << size << " partitions on partition " << rank << ": " << globalIndex )
    return globalIndex;
}

IndexType CyclicDistribution::allGlobal2local( const IndexType globalIndex ) const
{
    IndexType size = mCommunicator->getSize();

    IndexType globalChunkIndex = globalIndex / mChunkSize;
    IndexType localChunkIndex = globalChunkIndex / size;
    IndexType localIndex = localChunkIndex * mChunkSize + globalIndex % mChunkSize;
    SCAI_LOG_TRACE( logger,
                    "global Index " << globalIndex << " is with chunkSize " << mChunkSize << " and " << size << " partitions: " << localIndex )
    return localIndex;
}

IndexType CyclicDistribution::global2local( const IndexType globalIndex ) const
{
    if( isLocal( globalIndex ) )
    {
        return allGlobal2local( globalIndex );
    }
    else
    {
        return nIndex;
    }
}

void CyclicDistribution::computeOwners(
    const std::vector<IndexType>& requiredIndexes,
    std::vector<PartitionId>& owners ) const
{
    IndexType size = mCommunicator->getSize();
    owners.clear();
    owners.reserve( requiredIndexes.size() );

    SCAI_LOG_INFO( logger, "compute " << requiredIndexes.size() << " owners for " << *this )

    for( size_t i = 0; i < requiredIndexes.size(); i++ )
    {
        IndexType globalChunkIndex = requiredIndexes[i] / mChunkSize;
        IndexType owner = globalChunkIndex % size;
        SCAI_LOG_TRACE( logger, "owner of global index " << requiredIndexes[i] << " is " << owner )
        owners.push_back( owner );
    }
}

namespace
{
bool checkChunkSize( const CyclicDistribution& d, IndexType chunkSize )
{
    return chunkSize == d.chunkSize();
}
}

bool CyclicDistribution::isEqual( const Distribution& other ) const
{
    if( this == &other )
    {
        return true;
    }

    const CyclicDistribution* cycOther = dynamic_cast<const CyclicDistribution*>( &other );

    if( cycOther )
    {
        return ( mGlobalSize == other.getGlobalSize() && checkChunkSize( *cycOther, mChunkSize ) );
    }

    return false;
}

void CyclicDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "CyclicDistribution(gsize=" << mGlobalSize << ",chunkSize=" << mChunkSize << ",numLocalChunks="
           << getNumLocalChunks() << ")";
}

void CyclicDistribution::printDistributionVector( std::string name ) const
{
    IndexType myRank = mCommunicator->getRank();
    IndexType parts = mCommunicator->getSize();

    IndexType totalNumChunks = getNumTotalChunks();

    if( myRank == MASTER ) // process 0 is MASTER process
    {
        std::ofstream file;
        file.open( ( name + ".part" ).c_str() );
        // print row - partition mapping
        IndexType actualProcess = 0;

        for( IndexType i = 0; i < totalNumChunks; ++i )
        {
            for( IndexType j = 0; j < mChunkSize; j++ )
            {
                file << actualProcess << std::endl;
            }

            actualProcess = ( actualProcess + 1 ) % parts;
        }

        file.close();
    }
}

} /* end namespace lama */

} /* end namespace scai */
