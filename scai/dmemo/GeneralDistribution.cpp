/**
 * @file GeneralDistribution.cpp
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
 * @brief GeneralDistribution.cpp
 * @author brandes
 * @date 25.02.2011
 */

// hpp
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CommunicationPlan.hpp>
#include <scai/dmemo/GlobalAddressingPlan.hpp>

// internal scai libraries
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/common/macros/assert.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/tracing.hpp>

// std
#include <algorithm>
#include <functional>

#define MASTER 0

namespace scai
{

using namespace hmemo;
using namespace utilskernel;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GeneralDistribution::logger, "Distribution.General" )

const char GeneralDistribution::theCreateValue[] = "GENERAL";

/* ---------------------------------------------------------------------- */

GeneralDistribution::GeneralDistribution(
    const IndexType globalSize,
    HArray<IndexType> myIndexes,
    bool checkFlag,
    const CommunicatorPtr comm ) :

    Distribution( globalSize, comm ),
    mLocal2Global( std::move( myIndexes ) )

{
    SCAI_LOG_INFO( logger, "GeneralDistribution( size = " << globalSize
                   << ", myIndexes = " << myIndexes.size() << ", comm = " << *comm )

    SCAI_ASSERT_DEBUG( HArrayUtils::validIndexes( myIndexes, globalSize ), "myIndexes contains illegal values" )

    // sort the array ascending in-place, should be very fast if already sorted

    HArrayUtils::sort( NULL, &mLocal2Global, mLocal2Global, true );

    fillIndexMap();

    if ( checkFlag )
    {
        SCAI_ASSERT_EQ_ERROR( mGlobalSize, comm->sum( mLocal2Global.size() ), 
                              "illegal general distribution, total number of owned indexes not global size" )

        // this method also verifies that each global index appears only once

        enableBlockDistributedOwners();
    }
}

/* ---------------------------------------------------------------------- */

GeneralDistribution::GeneralDistribution( GeneralDistribution&& other ) :

    Distribution( other.getGlobalSize(), other.getCommunicatorPtr() ),
    mLocal2Global( std::move( other.mLocal2Global ) ),
    mGlobal2Local( std::move( other.mGlobal2Local ) )
{
    mBlockDistributedOwners.reset();
    mAnyAddressing.reset();
}    

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GeneralDistribution> generalDistribution( 
    const IndexType globalSize,
    HArray<IndexType> myGlobalIndexes,
    const CommunicatorPtr comm )
{
    return std::make_shared<GeneralDistribution>( globalSize, std::move( myGlobalIndexes ), true, comm );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GeneralDistribution> generalDistribution( 
    HArray<IndexType> myGlobalIndexes,
    const CommunicatorPtr comm )
{
    const bool checkFlag = true;
    const IndexType globalSize = comm->sum( myGlobalIndexes.size() );
    return std::make_shared<GeneralDistribution>( globalSize, std::move( myGlobalIndexes ), checkFlag, comm );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GeneralDistribution> generalDistributionUnchecked( 
    const IndexType globalSize,
    HArray<IndexType> myGlobalIndexes,
    const CommunicatorPtr comm )
{
    return std::make_shared<GeneralDistribution>( globalSize, std::move( myGlobalIndexes ), false, comm );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GeneralDistribution> generalDistributionUnchecked( std::shared_ptr<const Distribution> dist )
{
    const GeneralDistribution* genDist = dynamic_cast<const GeneralDistribution*>( dist.get() );

    if ( genDist )
    {
        if ( dist.use_count() == 1 )
        {
            // we can reuse member variables as there was only one reference
            GeneralDistribution& mutGenDist = const_cast<GeneralDistribution&>( *genDist );
            return std::make_shared<GeneralDistribution>( std::move( mutGenDist ) );
        }
        else
        {
            const HArray<IndexType>& ownedIndexes = genDist->getMyIndexes();
            return generalDistributionUnchecked( genDist->getGlobalSize(), ownedIndexes, genDist->getCommunicatorPtr() );
        }
    }
    else
    {
        return generalDistributionUnchecked( dist->getGlobalSize(), dist->ownedGlobalIndexes(), dist->getCommunicatorPtr() );
    }
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GeneralDistribution> generalDistributionBySingleOwners(
    const HArray<PartitionId>& owners,
    const PartitionId root,
    CommunicatorPtr comm ) 
{
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    SCAI_ASSERT_VALID_INDEX_ERROR( root, size, "illegal root processor specified" )

    IndexType globalSize[ 2 ];

    HArray<IndexType> localSizes;         // only defined by root
    HArray<IndexType> sortedIndexes;      // only defined by root

    if ( rank == root )
    {
        HArrayUtils::bucketSortSizes( localSizes, sortedIndexes, owners, size );

        globalSize[0] = owners.size();
        globalSize[1] = sortedIndexes.size();
    }

    comm->bcast( globalSize, 2, root );

    SCAI_ASSERT_EQ_ERROR( globalSize[0], globalSize[1], "illegal owners specified" )

    HArray<IndexType> myGlobalIndexes;

    comm->scatterVArray( myGlobalIndexes, root, sortedIndexes, localSizes );

    bool checkFlag = false;  // no check required

    return std::make_shared<GeneralDistribution>( globalSize[0], std::move( myGlobalIndexes ), checkFlag, comm );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<const GeneralDistribution> generalDistributionByNewOwners(
    const Distribution& dist,
    const hmemo::HArray<PartitionId>& newOwners )
{
    const IndexType globalSize = dist.getGlobalSize();

    CommunicatorPtr comm = dist.getCommunicatorPtr();

    const IndexType NP = comm->getSize();

    SCAI_ASSERT_DEBUG( HArrayUtils::validIndexes( newOwners, NP ), "illegal owners, #processors = " << NP )

    HArray<IndexType> myNewGlobalIndexes;  // new global indexes for this processor

    if ( NP == 1 )
    {
        SCAI_ASSERT_DEBUG( HArrayUtils::validIndexes( newOwners, NP ), "illegal owners, #partitions = " << NP )
        HArrayUtils::setOrder( myNewGlobalIndexes, globalSize );
    }
    else
    {
        // use a global exchange plan, returns exactly the new global indexes (unsorted) on each processor 

        HArray<IndexType> myOldGlobalIndexes;
        dist.getOwnedIndexes( myOldGlobalIndexes );
        dmemo::globalExchange( myNewGlobalIndexes, myOldGlobalIndexes, newOwners, comm );
    }

    bool checkFlag = false;  // cannot be illegal if new owners has only valid indexes

    return std::make_shared<GeneralDistribution>( globalSize, std::move( myNewGlobalIndexes ), checkFlag, comm );
}

/* ---------------------------------------------------------------------- */

GeneralDistribution::~GeneralDistribution()
{
    if ( static_cast<IndexType>( mGlobal2Local.size() ) != mLocal2Global.size() )
    {
        SCAI_LOG_ERROR( logger, "Oops: size mismatch: mGlobal2Local: " << mGlobal2Local.size()
                        << ", mLocal2Global : " << mLocal2Global.size() )
    }

    SCAI_LOG_DEBUG( logger, "~GeneralDistribution" )
}

/* ---------------------------------------------------------------------- */

bool GeneralDistribution::isLocal( const IndexType globalIndex ) const
{
    return mGlobal2Local.count( globalIndex ) > 0;
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getLocalSize() const
{
    return mLocal2Global.size();
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::fillIndexMap()
{
    SCAI_REGION( "Distribution.General.fillIndexMap" )
 
    IndexType nLocal = mLocal2Global.size();

    SCAI_LOG_INFO( logger, "fillHashMap, #local indexes = " << nLocal )

    // add my owned global indexes into unordered_map, [globalIndex : localIndex] 

    auto rIndexes = hostReadAccess( mLocal2Global );

    for ( IndexType i = 0; i < nLocal; i++ ) 
    {
        mGlobal2Local[rIndexes[i]] = i;
    }
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::local2Global( const IndexType localIndex ) const
{
    return mLocal2Global[localIndex];
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::global2Local( const IndexType globalIndex ) const
{
    // do a binary search in the array of global indexes for entries owned by this partition
    // return HArrayUtils::findPosInSortedIndexes( mLocal2Global, globalIndex );

    auto pos = mGlobal2Local.find( globalIndex );

    if ( pos == mGlobal2Local.end() )
    {
        return invalidIndex;
    } 
    else 
    {
        return pos->second;
    }
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::global2LocalV( hmemo::HArray<IndexType>& localIndexes, const hmemo::HArray<IndexType>& globalIndexes ) const
{
    SCAI_REGION( "Distribution.General.global2LocalV" )

    IndexType nnz = globalIndexes.size();

    ReadAccess<IndexType> rGlobal( globalIndexes );
    WriteOnlyAccess<IndexType> wLocal( localIndexes, nnz );

    #pragma omp parallel for

    for ( IndexType i = 0; i < nnz; ++i )
    {
        // avoid virtual call by calling method of this class directly
        wLocal[i] = GeneralDistribution::global2Local( rGlobal[i] );
    }

    // do a binary search in the array of global indexes for entries owned by this partition
    // HArrayUtils::findPosInSortedIndexesV( localIndexes, mLocal2Global, globalIndexes );
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getBlockDistributionSize() const
{
    // Note: we assume that the local indexes are descending and do not contain doubles

    IndexType localSize = mLocal2Global.size();

    bool isBlocked = true;

    ReadAccess<IndexType> rIndexes( mLocal2Global );

    for ( IndexType i = 0; i < localSize; ++i )
    {
        if ( rIndexes[i] - rIndexes[0] != i )
        {
            isBlocked = false;
        }
    }

    CommunicatorPtr comm = getCommunicatorPtr();

    isBlocked = comm->all( isBlocked );

    if ( !isBlocked )
    {
        return invalidIndex;
    }

    // Each processor has a contiguous part, but verify that it is in the same order

    auto genBlock = genBlockDistributionBySize( localSize, comm );

    IndexType lb = genBlock->lb();
    IndexType ub = genBlock->ub();

    isBlocked = true;

    if ( localSize > 0 )
    {
        isBlocked = ( rIndexes[0] == lb ) && ( rIndexes[localSize - 1] == ub - 1 );
    }

    isBlocked = comm->all( isBlocked );

    if ( !isBlocked )
    {
        return invalidIndex;
    }

    return localSize;
}

/* ---------------------------------------------------------------------- */

bool GeneralDistribution::isEqual( const Distribution& other ) const
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

    const GeneralDistribution& otherGen = reinterpret_cast<const GeneralDistribution&>( other );

    bool localSameSize = otherGen.getLocalSize() == getLocalSize();

    bool allSameSize = mCommunicator->all( localSameSize );

    SCAI_LOG_DEBUG( logger, *this << ": localSameSize = " << localSameSize << ", allSameSize = " << allSameSize )

    if ( !allSameSize )
    {
        return false;
    }

    // values will only be compared it sizes are same on all processors

    bool localSameVals = HArrayUtils::maxDiffNorm( mLocal2Global, otherGen.getMyIndexes() ) == 0;
    bool allSameVals   = mCommunicator->all( localSameVals );

    SCAI_LOG_DEBUG( logger, *this << ": localSameVals = " << localSameVals << ", allSameVals = " << allSameVals )

    return allSameVals;
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "GeneralDistribution( size = " << mLocal2Global.size() << " of " << mGlobalSize << ", comm = "
           << *mCommunicator << " )";
}

/* ---------------------------------------------------------------------- */

static void setOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes, const HArray<IndexType>& offsets )
{
    IndexType globalSize = indexes.size();
    IndexType nOwners    = offsets.size() - 1;

    WriteOnlyAccess<PartitionId> wOwners( owners, indexes.size() );
    ReadAccess<IndexType> rOffsets( offsets );
    ReadAccess<IndexType> rIndexes( indexes );

    // for testing: init

    for ( IndexType i = 0; i < globalSize; ++i )
    {
        wOwners[i] = invalidPartition;
    }

    for ( IndexType owner = 0; owner < nOwners; ++owner )
    {
        for ( IndexType j = rOffsets[owner]; j < rOffsets[owner + 1]; ++j )
        {
            wOwners[rIndexes[j]] = owner;
        }
    }
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::computeBlockDistributedOwners( HArray<IndexType>& blockDistributedOwners,
                                                         const IndexType globalSize, 
                                                         const HArray<IndexType>& ownedIndexes, 
                                                         CommunicatorPtr comm )
{
    SCAI_REGION( "Distribution.General.blockDistOwners" )

    // get my range for block distribution of size

    BlockDistribution blockDist( globalSize, comm );

    const IndexType localBlockSize = blockDist.getLocalSize();

    HArrayUtils::setSameValue( blockDistributedOwners, localBlockSize, invalidPartition );

    // blockDistributedOwners[ ownedIndexes[i] ] = rank, global scatter

    bool unique = true;  // no global index appears twice 'globally' seen

    auto plan = globalAddressingPlan ( blockDist, ownedIndexes, unique );

    plan.scatterOwner( blockDistributedOwners );

    // Sanitize check there must be no invalidPartition value for any global index

    bool okay = HArrayUtils::allScalar( blockDistributedOwners, common::CompareOp::NE, invalidPartition );

    okay = comm->all( okay );

    if ( !okay )
    {
        COMMON_THROWEXCEPTION( "Illegal general distribution, at least one index does not appear" )
    }
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::enableBlockDistributedOwners() const
{
    SCAI_LOG_INFO( logger, "enableBlockDistributedOwners" )

    if ( mBlockDistributedOwners.get() )
    {
        return;
    }

    mBlockDistributedOwners.reset( new HArray<PartitionId>() );
    computeBlockDistributedOwners( *mBlockDistributedOwners, getGlobalSize(), mLocal2Global, getCommunicatorPtr() );
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::disableBlockDistributedOwners() const
{
    SCAI_LOG_INFO( logger, "disableBlockDistributedOwners" )

    mBlockDistributedOwners.reset();
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::computeOwners(
    hmemo::HArray<PartitionId>& owners, 
    const hmemo::HArray<IndexType>& indexes ) const
{
    enableBlockDistributedOwners();

    SCAI_LOG_INFO( logger, "computeOwners for indexes = " << indexes )

    SCAI_REGION( "Distribution.General.computeOwners" )

    // initialize the array with invalidPartition, so we can easily identify non-available owners

    HArrayUtils::setSameValue( owners, indexes.size(), invalidPartition );

    BlockDistribution blockDist( getGlobalSize(), getCommunicatorPtr() );

    // owners[i] = blockDistribtedOwners[ indexes[i] ], global gather 

    auto plan = globalAddressingPlan( blockDist, indexes );
    
    plan.gather( owners, *mBlockDistributedOwners );
    
    SCAI_LOG_DEBUG( logger, "computeOwners done, owners = " << owners )
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::allOwners( HArray<PartitionId>& owners, const PartitionId root ) const
{
    ContextPtr ctx = Context::getHostPtr();

    PartitionId rank  = mCommunicator->getRank();
    PartitionId parts = mCommunicator->getSize();

    IndexType localSize = getLocalSize();

    // gather number of local rows

    HArray<IndexType> localSizes;

    {
        WriteOnlyAccess<IndexType> wLocalSizes( localSizes, ctx, parts );
        mCommunicator->gather( wLocalSizes.get(), 1, root, &localSize );
    }

    HArray<IndexType> offsets;
    offsets.reserve( ctx, parts + 1 );;

    if ( rank == root )
    {
        HArrayUtils::assign( offsets, localSizes, ctx );
        IndexType nTotal = HArrayUtils::scan1( offsets );
        SCAI_ASSERT( nTotal == mGlobalSize, "sum of local rows is not global size" )
    }

    // gather global indices of local rows

    HArray<IndexType> allIndexes( rank == root ? mGlobalSize : 2 );

    {
        WriteAccess<IndexType> wAllIndexes( allIndexes );
        ReadAccess<IndexType> rLocalSizes( localSizes );
        ReadAccess<IndexType> rIndexes( mLocal2Global );
        mCommunicator->gatherV( wAllIndexes.get(), localSize, root, rIndexes.get(), rLocalSizes.get() );
    }

    // now we ca scatter the owner ids in the owner array

    if ( rank == root )
    {
        setOwners( owners, allIndexes, offsets );
    }
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const
{
    HArrayUtils::assign( myGlobalIndexes, mLocal2Global );
}

/* ---------------------------------------------------------------------- */

bool GeneralDistribution::hasAnyAddressing() const
{
    return mAnyAddressing.get() != nullptr;
}

void GeneralDistribution::enableAnyAddressing() const
{
    if ( mAnyAddressing.get() ) 
    {
        return;   // already done
    }

    mAnyAddressing.reset( new AnyAddressing( *this ) );
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getAnyLocalSize( const PartitionId rank ) const
{
    SCAI_ASSERT( mAnyAddressing.get(), "any addressing not enabled" )

    return mAnyAddressing->localSize( rank );
}

/* ---------------------------------------------------------------------- */

PartitionId GeneralDistribution::getAnyOwner( const IndexType globalIndex ) const
{
    SCAI_ASSERT( mAnyAddressing.get(), "any addressing not enabled" )

    return mAnyAddressing->owner( globalIndex );
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const
{
    SCAI_ASSERT( mAnyAddressing.get(), "any addressing not enabled" )

    // here the owner is important as local index requires size offsets

    return mAnyAddressing->localIndex( globalIndex, owner );
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getAnyGlobalIndex( const IndexType localIndex, const PartitionId owner ) const
{
    SCAI_ASSERT( mAnyAddressing.get(), "any addressing not enabled" )

    return mAnyAddressing->globalIndex( localIndex, owner );
}

} /* end namespace dmemo */

} /* end namespace scai */
