/**
 * @file CyclicDistribution.cpp
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
    IndexType rank = mCommunicator->getRank();

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

void CyclicDistribution::getChunkInfo( IndexType& localChunks, IndexType& extra, const PartitionId rank ) const
{
    const PartitionId size = mCommunicator->getSize();
    const IndexType chunks = mGlobalSize / mChunkSize;
    // mGlobalSize % mChunkSize elements remain, are counted later
    localChunks = chunks / size;
    IndexType remainChunks = chunks % size;
    extra = 0;

    if ( rank < remainChunks )
    {
        localChunks++;
    }
    else if ( rank == remainChunks )
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
    if ( isLocal( globalIndex ) )
    {
        return allGlobal2local( globalIndex );
    }
    else
    {
        return nIndex;
    }
}

/* ---------------------------------------------------------------------- */

void CyclicDistribution::computeOwners1(
    const std::vector<IndexType>& requiredIndexes,
    std::vector<PartitionId>& owners ) const
{
    IndexType size = mCommunicator->getSize();
    owners.clear();
    owners.reserve( requiredIndexes.size() );
    SCAI_LOG_INFO( logger, "compute " << requiredIndexes.size() << " owners for " << *this )

    for ( size_t i = 0; i < requiredIndexes.size(); i++ )
    {
        IndexType globalChunkIndex = requiredIndexes[i] / mChunkSize;
        IndexType owner = globalChunkIndex % size;
        SCAI_LOG_TRACE( logger, "owner of global index " << requiredIndexes[i] << " is " << owner )
        owners.push_back( owner );
    }
}

/* ---------------------------------------------------------------------- */

void CyclicDistribution::computeOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes ) const
{   
    ContextPtr ctx = Context::getHostPtr();    // currently only available @ Host
    
    const IndexType n = indexes.size();
    const IndexType size = mCommunicator->getSize();
    
    ReadAccess<IndexType> rIndexes( indexes, ctx );
    WriteOnlyAccess<PartitionId> wOwners( owners, ctx, n );
    
    // ToDo: call a kernel and allow arbitrary context

    for ( IndexType i = 0; i < n; i++ )
    {   
        wOwners[i] = ( rIndexes[i] / mChunkSize ) % size;
    }
}

/* ---------------------------------------------------------------------- */

namespace
{
bool checkChunkSize( const CyclicDistribution& d, IndexType chunkSize )
{
    return chunkSize == d.chunkSize();
}
}

bool CyclicDistribution::isEqual( const Distribution& other ) const
{
    if ( this == &other )
    {
        return true;
    }

    const CyclicDistribution* cycOther = dynamic_cast<const CyclicDistribution*>( &other );

    if ( cycOther )
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

    if ( myRank == MASTER ) // process 0 is MASTER process
    {
        std::ofstream file;
        file.open( ( name + ".part" ).c_str() );
        // print row - partition mapping
        IndexType actualProcess = 0;

        for ( IndexType i = 0; i < totalNumChunks; ++i )
        {
            for ( IndexType j = 0; j < mChunkSize; j++ )
            {
                file << actualProcess << std::endl;
            }

            actualProcess = ( actualProcess + 1 ) % parts;
        }

        file.close();
    }
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string CyclicDistribution::createValue()
{
    return "CYCLIC";
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
