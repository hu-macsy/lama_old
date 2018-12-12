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
#include <scai/dmemo/GlobalExchangePlan.hpp>

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
    const CommunicatorPtr comm ) :

    Distribution( globalSize, comm ),
    mLocal2Global( std::move( myIndexes ) )

{
    SCAI_LOG_INFO( logger, "GeneralDistribution( size = " << globalSize
                   << ", myIndexes = " << myIndexes.size() << ", comm = " << *comm )

    SCAI_ASSERT_DEBUG( HArrayUtils::validIndexes( myIndexes, globalSize ), "myIndexes contains illegal values" )

    // Note: the constructor is completely local, but make some consistency check now

    IndexType nLocal = mLocal2Global.size();

    SCAI_ASSERT_EQ_ERROR( mGlobalSize, comm->sum( nLocal ), "illegal general distribution" )

    // sort the array ascending, should be very fast if already sorted

    HArrayUtils::sort( NULL, &mLocal2Global, mLocal2Global, true );

    fillIndexMap();
    setBlockDistributedOwners();
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<GeneralDistribution> generalDistribution( 
    const IndexType globalSize,
    HArray<IndexType> myGlobalIndexes,
    const CommunicatorPtr comm )
{
    // here we could sort the global indexes in place

    return std::make_shared<GeneralDistribution>( globalSize, std::move( myGlobalIndexes ), comm );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<GeneralDistribution> generalDistributionByOwners(
    const HArray<PartitionId>& owners,
    const PartitionId root,
    CommunicatorPtr comm ) 
{
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    SCAI_ASSERT_VALID_INDEX_ERROR( root, size, "illegal root processor specified" )

    IndexType globalSize;

    if ( rank == root )
    {
        globalSize = owners.size();
    }

    comm->bcast( &globalSize, 1, root );

    HArray<IndexType> localSizes;         // only defined by root
    HArray<IndexType> sortedIndexes;      // only defined by root

    if ( rank == root )
    {
        HArrayUtils::bucketSortSizes( localSizes, sortedIndexes, owners, size );
    }

    HArray<IndexType> myGlobalIndexes;

    comm->scatterVArray( myGlobalIndexes, root, sortedIndexes, localSizes );

    return std::make_shared<GeneralDistribution>( globalSize, std::move( myGlobalIndexes ), comm );
}

/* ---------------------------------------------------------------------- */

std::shared_ptr<GeneralDistribution> generalDistributionNew(
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
        // use a global exchange pattern, very common for many operations

        HArray<IndexType> myOldGlobalIndexes;
        dist.getOwnedIndexes( myOldGlobalIndexes );
        dmemo::globalExchange( myNewGlobalIndexes, myOldGlobalIndexes, newOwners, *comm );
    }

    return std::make_shared<GeneralDistribution>( globalSize, std::move( myNewGlobalIndexes ), comm );
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

IndexType GeneralDistribution::local2global( const IndexType localIndex ) const
{
    return mLocal2Global[localIndex];
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::global2local( const IndexType globalIndex ) const
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

void GeneralDistribution::global2localV( hmemo::HArray<IndexType>& localIndexes, const hmemo::HArray<IndexType>& globalIndexes ) const
{
    SCAI_REGION( "Distribution.General.global2localV" )

    IndexType nnz = globalIndexes.size();

    ReadAccess<IndexType> rGlobal( globalIndexes );
    WriteOnlyAccess<IndexType> wLocal( localIndexes, nnz );

    #pragma omp parallel for

    for ( IndexType i = 0; i < nnz; ++i )
    {
        // avoid virtual call by calling method of this class directly
        wLocal[i] = GeneralDistribution::global2local( rGlobal[i] );
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

    auto genBlock = genBlockDistribution( localSize, comm );

    IndexType lb;
    IndexType ub;

    genBlock->getLocalRange( lb, ub );

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
    const IndexType lb             = blockDist.lb();

    HArrayUtils::setSameValue( blockDistributedOwners, localBlockSize, invalidPartition );

    // each processor queries the owners of its range lb, ..., ub-1

    HArray<PartitionId> blockOwners;

    blockDist.computeOwners( blockOwners, ownedIndexes );
 
    GlobalExchangePlan plan( blockOwners, *comm );

    HArray<IndexType> receivedIndexes;

    plan.exchange( receivedIndexes, ownedIndexes, *comm );

    HArray<PartitionId> sources;  // sources[i] is the onwer of receivedIndexes[i]

    plan.getSource( sources );  // sources[i] is the onwer of receivedIndexes[i]

    // receivedIndexes are all values in my block range lb, ..., ub

    HArrayUtils::compute( receivedIndexes, receivedIndexes, common::BinaryOp::SUB, lb );

    bool unique = true;  // there should be no double values in receivedIndexes

    HArrayUtils::scatter( blockDistributedOwners, receivedIndexes, unique, sources, common::BinaryOp::COPY );

    // Sanitize check there must be no invalidPartition value for any global index

    bool okay = HArrayUtils::allScalar( blockDistributedOwners, common::CompareOp::NE, invalidPartition );

    okay = comm->all( okay );

    if ( !okay )
    {
        COMMON_THROWEXCEPTION( "Illegal general distribution, at least one index does not appear" )
    }
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::setBlockDistributedOwners()
{
    SCAI_LOG_INFO( logger, "setBlockDistributedOwners" )
    computeBlockDistributedOwners( mBlockDistributedOwners, getGlobalSize(), mLocal2Global, getCommunicatorPtr() );
    SCAI_LOG_DEBUG( logger, "setBlockDistributedOwners done" )
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::computeOwners(
    hmemo::HArray<PartitionId>& owners, 
    const hmemo::HArray<IndexType>& indexes ) const
{
    SCAI_LOG_INFO( logger, "computeOwners for indexes = " << indexes )

    SCAI_REGION( "Distribution.General.computeOwners" )

    // initialize the array with invalidPartition, so we can easily identify non-available owners

    HArrayUtils::setSameValue( owners, indexes.size(), invalidPartition );

    BlockDistribution blockDist( getGlobalSize(), getCommunicatorPtr() );

    HArray<PartitionId> target;    // target processor that we will query for the owner

    // get for each queried index the processor that knows the owner

    blockDist.computeOwners( target, indexes );

    const Communicator& comm = getCommunicator();

    GlobalExchangePlan plan( target, comm );

    SCAI_LOG_DEBUG( logger, "computeOwners: plan for queries ready." )

    HArray<IndexType> queriedIndexes;  // that are queried indexes for which I know the owner

    plan.exchange( queriedIndexes, indexes, comm );

    SCAI_LOG_DEBUG( logger, "queries received: " << queriedIndexes )

    HArrayUtils::compute( queriedIndexes, queriedIndexes, common::BinaryOp::SUB, blockDist.lb() );

    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( queriedIndexes, blockDist.getLocalSize() ), "serious wrong index" )

    SCAI_LOG_DEBUG( logger, "computeOwners: now compute answers" )

    HArray<PartitionId> answerOwners;  // will take the owners of the queried indexes

    HArrayUtils::gather( answerOwners, mBlockDistributedOwners, queriedIndexes, common::BinaryOp::COPY );

    SCAI_LOG_DEBUG( logger, "computeOwners: answers = " << answerOwners << ", exchange back" )

    plan.exchangeBack( owners, answerOwners, comm );

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
    return mAllOwners.size() > 0;
}

void GeneralDistribution::enableAnyAddressing() const
{
    if ( mAllOwners.size() > 0 )
    {
        // already computed, but just verify correct sizes

        SCAI_ASSERT_EQ_DEBUG( mAllOwners.size(), getGlobalSize(), "serious mismatch" )
        SCAI_ASSERT_EQ_DEBUG( mAllLocalOffsets.size(), getCommunicator().getSize() + 1, "serious mismatch" )
        SCAI_ASSERT_EQ_DEBUG( mAllLocal2Global.size(), getGlobalSize(), "serious mismatch" )
        SCAI_ASSERT_EQ_DEBUG( mAllGlobal2Local.size(), getGlobalSize(), "serious mismatch" )

        return;   // already done
    }

    // compute mAllOwners

    HArray<IndexType> indexes;   // will contain all column indexes to get all owners

    HArrayUtils::setOrder( indexes, mGlobalSize );

    Distribution::computeOwners( mAllOwners, indexes );

    // bucket sort the owners, gives offsets and permutation to block values according to owners

    HArrayUtils::bucketSortOffsets( mAllLocalOffsets, mAllLocal2Global, mAllOwners, mCommunicator->getSize() );
    HArrayUtils::inversePerm( mAllGlobal2Local, mAllLocal2Global ); // global2local
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getAnyLocalSize( const PartitionId rank ) const
{
    SCAI_ASSERT( mAllLocalOffsets.size() > 0, "any addressing not enabled" )

    return mAllLocalOffsets[ rank + 1] - mAllLocalOffsets[rank];
}

/* ---------------------------------------------------------------------- */

PartitionId GeneralDistribution::getAnyOwner( const IndexType globalIndex ) const
{
    SCAI_ASSERT( mAllOwners.size() > 0, "any addressing not enabled" )

    return mAllOwners[ globalIndex ];
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const
{
    SCAI_ASSERT( mAllGlobal2Local.size() > 0, "any addressing not enabled" )

    // here the owner is important as local index  requires size offsets

    return mAllGlobal2Local[ globalIndex ] - mAllLocalOffsets[ owner ];
}

/* ---------------------------------------------------------------------- */

IndexType GeneralDistribution::getAnyGlobalIndex( const IndexType localIndex, const PartitionId owner ) const
{
    SCAI_ASSERT( mAllGlobal2Local.size() > 0, "any addressing not enabled" )

    // here the owner is important as local index  requires size offsets

    return mAllLocal2Global[ localIndex + mAllLocalOffsets[ owner ] ];
}

} /* end namespace dmemo */

} /* end namespace scai */
