/**
 * @file BlockDistribution.cpp
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

    if ( mLB <= mUB )
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

    if ( globalIndex >= mLB && globalIndex <= mUB )
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

bool BlockDistribution::isEqual( const Distribution& other ) const
{
    if ( this == &other )
    {
        return true;
    }

    if ( dynamic_cast<const BlockDistribution*>( &other ) )
    {
        return mGlobalSize == other.getGlobalSize();
    }

    return false;
}

void BlockDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "BlockDistribution( comm = " << *mCommunicator << ", block = " << mBlockSize
           << ", size = " << mLB << ":" << mUB << " of " << mGlobalSize <<  " )";
}

void BlockDistribution::printDistributionVector( std::string name ) const
{
    PartitionId myRank = mCommunicator->getRank();
    PartitionId parts = mCommunicator->getSize();
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

const char BlockDistribution::theCreateValue[] = "BLOCK";

std::string BlockDistribution::createValue()
{
    return theCreateValue;
}

Distribution* BlockDistribution::create( const DistributionArguments arg )
{
    SCAI_LOG_INFO( logger, "create" )
    // Note: weight argument is not used here
    return new BlockDistribution( arg.globalSize, arg.communicator );
}

} /* end namespace dmemo */

} /* end namespace scai */
