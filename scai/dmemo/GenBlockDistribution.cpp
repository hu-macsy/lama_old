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
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

// std
#include <fstream>
#include <memory>

/** 
 *  common definition for root processor 0
 */
#define MASTER 0

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GenBlockDistribution::logger, "Distribution.GenBlockDistribution" )

/* ---------------------------------------------------------------------- */

GenBlockDistribution::~GenBlockDistribution()
{
    SCAI_LOG_INFO( logger, "~GenBlockDistribution" )
}

/* ---------------------------------------------------------------------- */

GenBlockDistribution::GenBlockDistribution( std::unique_ptr<IndexType[]>&& offsets, CommunicatorPtr comm ) :
  
    Distribution( offsets[comm->getSize()], comm ),
    mOffsets( std::move( offsets ) )
{
    PartitionId rank = comm->getRank();
    mLB = mOffsets[ rank ];
    mUB = mOffsets[ rank + 1 ];
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GenBlockDistribution> genBlockDistributionBySize( 
    const IndexType localSize, 
    CommunicatorPtr comm )
{
    SCAI_ASSERT_DEBUG( comm, "Null pointer for commmunicator" )

    PartitionId size = comm->getSize();

    std::unique_ptr<IndexType[]> offsets( new IndexType[ size + 1 ] );

    offsets[0] = 0;  // local sizes will be added after this entry

    IndexType* allLocalSizes = offsets.get() + 1;

    comm->gather( allLocalSizes, 1, 0, &localSize );
    comm->bcast( allLocalSizes, size, 0 );

    SCAI_ASSERT_EQ_DEBUG( allLocalSizes[comm->getRank()], localSize, "wrongly gathered values" )

    // sizes to offsets, e.g.  0  5  3  4  1  ->  0  5  8  12 13

    for ( PartitionId p = 0; p < size; ++p )
    {
        offsets[ p + 1 ] += offsets[p];
    }

    return std::make_shared<GenBlockDistribution>( std::move( offsets ), comm );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GenBlockDistribution> genBlockDistributionByOffset(
    const IndexType N,
    const IndexType offset,
    CommunicatorPtr comm )
{
    SCAI_ASSERT_DEBUG( comm, "Null pointer for commmunicator" )

    PartitionId size = comm->getSize();

    std::unique_ptr<IndexType[]> offsets( new IndexType[ size + 1 ] );

    comm->gather( offsets.get(), 1, MASTER, &offset );
    comm->bcast( offsets.get(), size, MASTER );

    offsets[size] = N;

    IndexType curOffset = N; 

    for ( PartitionId p = size; p-- > 0; )
    {
        if ( offsets[p] == invalidIndex )
        {
            // correct the offsets for invalidIndex (no elements)
            offsets[p] = curOffset;
        }
        else
        {
            // check that offset for p is valid 
            SCAI_ASSERT_LE_ERROR( offsets[p], curOffset, "invalid offset for processor p = " << p )
            curOffset = offsets[p];
        }
    }

    offsets[0] = 0;  

    return std::make_shared<GenBlockDistribution>( std::move( offsets ), comm );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GenBlockDistribution> genBlockDistributionBySize( 
    const IndexType globalSize,
    const IndexType localSize, 
    CommunicatorPtr comm )
{
    auto dist = genBlockDistributionBySize( localSize, comm );
    SCAI_ASSERT_EQUAL( dist->getGlobalSize(), globalSize, "local sizes do not sum up to expected global size" )
    return dist;
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GenBlockDistribution> genBlockDistributionBySizes( 
    const std::vector<IndexType>& localSizes,
    CommunicatorPtr comm )
{
    SCAI_ASSERT_DEBUG( comm, "Null pointer for commmunicator" )

    PartitionId size = comm->getSize();

    SCAI_ASSERT_EQ_DEBUG( static_cast<PartitionId>( localSizes.size() ), size, "illegal vector for local sizes" )

    std::unique_ptr<IndexType[]> offsets( new IndexType[ size + 1 ] );

    offsets[0] = 0;  // local sizes will be added after this entry

    IndexType* allLocalSizes = offsets.get() + 1;

    for ( PartitionId p = 0; p < size; ++p )
    {
        allLocalSizes[p] = localSizes[p];
    }

    // sizes to offsets, e.g.  0  5  3  4  1  ->  0  5  8  12 13

    for ( PartitionId p = 0; p < size; ++p )
    {
        offsets[ p + 1 ] += offsets[p];
    }

    // we make some simple test to verify that all processors have used the same array

    SCAI_ASSERT_EQ_ERROR( comm->max( offsets[size] ), offsets[size], "different global sizes" )

    return std::make_shared<GenBlockDistribution>( std::move( offsets ), comm );
}

/* ---------------------------------------------------------------------- */

bool GenBlockDistribution::isLocal( const IndexType globalIndex ) const
{
    return globalIndex >= mLB && globalIndex < mUB;
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

IndexType GenBlockDistribution::local2Global( const IndexType localIndex ) const
{
    return mLB + localIndex;
}

/* ---------------------------------------------------------------------- */

IndexType GenBlockDistribution::global2Local( const IndexType globalIndex ) const
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

bool GenBlockDistribution::isSameGenBlockDistribution( const GenBlockDistribution& other ) const
{
    // make at least sure that both distribution have same processor set

    if ( getCommunicator() != other.getCommunicator() )
    {
        return false;
    }

    PartitionId np = getCommunicator().getSize();

    bool equal = true;

    for ( PartitionId p = 0; p < np + 1; ++p )
    {
        if ( mOffsets[p] != other.mOffsets[p] )
        {
            equal = false;
            break;
        }
    }

    return equal;
}

/* ---------------------------------------------------------------------- */

bool GenBlockDistribution::isBlockDistribution() const
{
    const PartitionId np = getCommunicator().getSize();
    const IndexType   N  = getGlobalSize();

    const IndexType blockSize = ( N + np - 1 ) / np;

    bool equal = true;

    for ( PartitionId p = 0; p < np + 1; ++p )
    {
        IndexType otherOffset = std::min( p * blockSize, N );

        if ( mOffsets[p] != otherOffset )
        {
            equal = false;
            break;
        }
    }

    return equal;
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

    if ( strcmp( other.getKind(), BlockDistribution::getId() ) == 0 )
    {
        // other is block distribution, we know already same size and same number of processors

        return isBlockDistribution();
    }

    if ( strcmp( other.getKind(), GenBlockDistribution::getId() ) == 0 )
    {
        return isSameGenBlockDistribution( static_cast<const GenBlockDistribution&>( other ) );
    }

    return isSame;
}

/* ---------------------------------------------------------------------- */

void GenBlockDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "GenBlockDistribution( comm = " << *mCommunicator
           << ", size = " << mLB << ":" << mUB << " of " << mGlobalSize << " )";
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string GenBlockDistribution::createValue()
{
    return getId();
}

DistributionPtr GenBlockDistribution::create( const DistributionArguments arg )
{
    if ( arg.matrix != NULL )
    {
        SCAI_LOG_WARN( logger, "matrix argument ignored to create GEN_BLOCK distribution" )
    }
    
    return genBlockDistributionByWeight( arg.globalSize, arg.weight, arg.communicator );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GenBlockDistribution> genBlockDistributionByWeight(
    const IndexType globalSize,
    const float weight,
    const CommunicatorPtr comm )
{
    PartitionId size = comm->getSize();

    // each processor gets all values of the weights

    std::vector<float> allWeights( size );

    comm->allgather( &allWeights[0], 1, &weight );

    float totalWeight = 0.0f;

    for ( PartitionId p = 0; p < size; p++ )
    {
        SCAI_ASSERT_GE_ERROR( allWeights[p], 0.0f, 
                              "weight of processor " << p << " illegal, must not be negative" );

        totalWeight += allWeights[p];
    }

    std::unique_ptr<IndexType[]> offsets( new IndexType[size + 1] );

    float sumWeight = 0.0f;

    offsets[0] = 0;

    for ( PartitionId p = 0; p < size; p++ )
    {
        sumWeight += allWeights[p];
        offsets[p + 1] = static_cast<IndexType>( sumWeight / totalWeight * globalSize + 0.5 );
    }

    return std::make_shared<GenBlockDistribution>( std::move( offsets ), comm );
}

} /* end namespace dmemo */

} /* end namespace scai */
