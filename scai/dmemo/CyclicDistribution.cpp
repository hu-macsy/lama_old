/**
 * @file CyclicDistribution.cpp
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
 * @brief CyclicDistribution.cpp
 * @author schubert
 * @date 20.05.2011
 */

// hpp
#include <scai/dmemo/CyclicDistribution.hpp>

// std
#include <fstream>

#define MASTER 0

namespace scai
{

using namespace hmemo;

namespace dmemo
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
    const PartitionId rank = mCommunicator->getRank();

    if ( getOwner( globalIndex ) == rank )
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

IndexType CyclicDistribution::getMaxLocalSize() const
{
    // processor 0 has always the most elemements
    return getPartitionSize( 0 );
}

void CyclicDistribution::getChunkInfo( IndexType& localChunks, IndexType& extra, const PartitionId rank ) const
{
    const PartitionId size = mCommunicator->getSize();
    const IndexType chunks = mGlobalSize / mChunkSize;
    // mGlobalSize % mChunkSize elements remain, are counted later
    localChunks = chunks / size;
    IndexType remainChunks = chunks % size;
    extra = 0;

    if ( static_cast<IndexType>( rank ) < remainChunks )
    {
        localChunks++;
    }
    else if ( static_cast<IndexType>( rank ) == remainChunks )
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

    if ( extra )
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

    if ( mGlobalSize % mChunkSize != 0 )
    {
        ++numChunks; // count also the non-full chunk
    }

    return numChunks;
}

/* ---------------------------------------------------------------------- */

IndexType CyclicDistribution::getPartitionSize( const PartitionId partition ) const
{
    IndexType localChunks = 0;
    IndexType extra = 0;
    getChunkInfo( localChunks, extra, partition );
    IndexType elements = localChunks * mChunkSize + extra;
    return elements;
}

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

IndexType CyclicDistribution::global2local( const IndexType globalIndex ) const
{
    if ( isLocal( globalIndex ) )
    {
        return allGlobal2local( globalIndex );
    }
    else
    {
        return invalidIndex;
    }
}

/* ---------------------------------------------------------------------- */

void CyclicDistribution::computeOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes ) const
{
    ContextPtr ctx = Context::getHostPtr();    // currently only available @ Host

    const IndexType n = indexes.size();
    const IndexType size = mCommunicator->getSize();

    SCAI_LOG_INFO( logger, *this << ": compute owners, n = " << n << ", size = " << size )

    ReadAccess<IndexType> rIndexes( indexes, ctx );
    WriteOnlyAccess<PartitionId> wOwners( owners, ctx, n );

    // ToDo: call a kernel and allow arbitrary context

    #pragma omp parallel for

    for ( IndexType i = 0; i < n; i++ )
    {
        wOwners[i] = ( rIndexes[i] / mChunkSize ) % size;   // see getOwner( i )
    }
}

/* ---------------------------------------------------------------------- */

void CyclicDistribution::getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const
{
    const IndexType nLocal  = getLocalSize();

    const Communicator& comm = getCommunicator();

    const PartitionId rank = comm.getRank();
    const PartitionId size = comm.getSize();

    SCAI_LOG_INFO( logger, comm << ": getOwnedIndexes, have " << nLocal << " of " << mGlobalSize )

    WriteOnlyAccess<IndexType> wGlobalIndexes( myGlobalIndexes, nLocal );

    IndexType pos   = 0;

    // we can directly enumerate the indexes owned by this processor

    for ( IndexType first = rank * mChunkSize; first < mGlobalSize; first += size * mChunkSize )
    {
        for ( IndexType j = 0; j < mChunkSize; j++ )
        {
            IndexType myIndex = first + j;

            if ( myIndex < mGlobalSize )
            {
                wGlobalIndexes[ pos++ ] = myIndex;
            }
        }
    }

    SCAI_ASSERT_EQ_ERROR( pos, nLocal, "count mismatch" )
}

/* ---------------------------------------------------------------------- */

IndexType CyclicDistribution::getBlockDistributionSize() const
{
    const PartitionId numPartitions = getCommunicator().getSize();

    if ( numPartitions == 1 )
    {
        return mGlobalSize;
    }

    // the following bool expression is evaluated by all processor with the same result

    if ( mGlobalSize <= mChunkSize * numPartitions )
    {
        return getLocalSize();
    }
    else
    {
        return invalidIndex;
    }
}

/* ---------------------------------------------------------------------- */

bool CyclicDistribution::hasAnyAddressing() const
{
    return true;
}

void CyclicDistribution::enableAnyAddressing() const
{
    // nothing to do as we can compute owners, local indexes by closed formulas
}

IndexType CyclicDistribution::getAnyLocalSize( const PartitionId rank ) const
{
    return getPartitionSize( rank );
}

PartitionId CyclicDistribution::getAnyOwner( const IndexType globalIndex ) const
{
    return getOwner( globalIndex );
}

IndexType CyclicDistribution::getAnyLocalIndex( const IndexType globalIndex, const PartitionId ) const
{
    IndexType size = mCommunicator->getSize();
    IndexType globalChunkIndex = globalIndex / mChunkSize;
    IndexType localChunkIndex = globalChunkIndex / size;
    IndexType localIndex = localChunkIndex * mChunkSize + globalIndex % mChunkSize;
    return localIndex;
}

IndexType CyclicDistribution::getAnyGlobalIndex( const IndexType localIndex, const PartitionId owner ) const
{
    IndexType size = mCommunicator->getSize();

    IndexType localChunk  = localIndex / mChunkSize;
    IndexType chunkPos    = localIndex - mChunkSize * localChunk; // pos in chunk
    IndexType globalChunk = localChunk * size + owner;
    return globalChunk + chunkPos;
}

/* ---------------------------------------------------------------------- */

bool CyclicDistribution::isEqual( const Distribution& other ) const
{
    bool isSame = false;

    bool proven = proveEquality( isSame, other );

    if ( proven )
    {
        return isSame;
    }

    if ( other.getKind() == getKind() )
    {
        const CyclicDistribution& cycOther = reinterpret_cast<const CyclicDistribution&>( other );

        isSame = chunkSize() == cycOther.chunkSize();
    }

    return isSame;
}

/* ---------------------------------------------------------------------- */

void CyclicDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "CyclicDistribution(gsize=" << mGlobalSize << ",chunkSize=" << mChunkSize << ",numLocalChunks="
           << getNumLocalChunks() << ")";
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string CyclicDistribution::createValue()
{
    return getId();
}

Distribution* CyclicDistribution::create( const DistributionArguments arg )
{
    if ( arg.matrix != NULL )
    {
        SCAI_LOG_WARN( logger, "matrix argument ignored to create CYCLIC distribution" )
    }

    return new CyclicDistribution( arg.globalSize, 1, arg.communicator );
}

} /* end namespace dmemo */

} /* end namespace scai */
