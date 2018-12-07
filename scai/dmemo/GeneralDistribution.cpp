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

GeneralDistribution::GeneralDistribution(
    const Distribution& other,
    const HArray<PartitionId>& owners ) :

    Distribution( other.getGlobalSize(), other.getCommunicatorPtr() )

{
    SCAI_LOG_INFO( logger, "GeneralDistribution( dist = " << other << ", new local owners =  " << owners )

    SCAI_ASSERT_EQ_DEBUG( other.getLocalSize(), owners.size(), "serious size mismatch" )

    const IndexType nLocal = owners.size();

    const Communicator& comm = *mCommunicator;

    const PartitionId numPartitions = comm.getSize();

    if ( numPartitions == 1 )
    {
        // owners[i] == 0 for al i, mLocal2Global = { 0, 1, ..., globalSize - 1 }

        SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( owners, 1 ), "illegal owners, #partitions = " << numPartitions )

        HArrayUtils::setOrder( mLocal2Global, nLocal );

        fillIndexMap();
        setBlockDistributedOwners();

        return;
    }

    // make a bucket sort with owners

    HArray<IndexType> offsets;    // number of elements for each bucket/partition
    HArray<IndexType> perm;       // permutation to resort in buckets

    HArrayUtils::bucketSortOffsets( offsets, perm, owners, mCommunicator->getSize() );

    HArray<IndexType> sortedIndexes;    // current indexes resorted according the buckets

    if ( other.getKind() == GeneralDistribution::theCreateValue )
    {
        SCAI_ASSERT_DEBUG( dynamic_cast<const GeneralDistribution*>( &other ), "no general dist" )
        const GeneralDistribution& otherGen = reinterpret_cast<const GeneralDistribution&>( other );
        HArrayUtils::gather( sortedIndexes, otherGen.getMyIndexes(), perm, common::BinaryOp::COPY );
    }
    else
    {
        HArray<IndexType> currentIndexes;
        other.getOwnedIndexes( currentIndexes );
        HArrayUtils::gather( sortedIndexes, currentIndexes, perm, common::BinaryOp::COPY );
    }

    // make communication plans for sending data and receiving to exchange

    auto sendPlan = CommunicationPlan::buildByOffsets( hostReadAccess( offsets ).get(), numPartitions );
    auto recvPlan = comm.transpose( sendPlan );

    SCAI_LOG_INFO( logger, comm << ": send plan: " << sendPlan << ", rev plan: " << recvPlan );

    // we just receive all the values in mLocal2Global

    IndexType newLocalSize = recvPlan.totalQuantity();

    HArray<IndexType> myNewIndexes;

    {
        WriteOnlyAccess<IndexType> recvVals( myNewIndexes, newLocalSize );
        ReadAccess<IndexType> sendVals( sortedIndexes );
        comm.exchangeByPlan( recvVals.get(), recvPlan, sendVals.get(), sendPlan );
    }

    // Important: the new local indexes must be sorted
    // Note: actually it would be sufficient to have a mergesort

    HArrayUtils::sort( NULL, &mLocal2Global, myNewIndexes, true );

    fillIndexMap();
    setBlockDistributedOwners();
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

void GeneralDistribution::setBlockDistributedOwners()
{
    SCAI_REGION( "Distribution.General.blockDistOwners" )

    // get my range for block distribution of size

    BlockDistribution blockDist( getGlobalSize(), getCommunicatorPtr() );

    const Communicator& comm = getCommunicator();
    const PartitionId   np   = comm.getSize();

    hmemo::ContextPtr hostCtx = hmemo::Context::getHostPtr();

    const IndexType localBlockSize = blockDist.getLocalSize();
    const IndexType myLB           = blockDist.lb();

    HArrayUtils::setSameValue( mBlockDistributedOwners, localBlockSize, invalidPartition, hostCtx );

    // get Owners of myIndexes according to the block distribution

    HArray<PartitionId> blockOwners;
    blockDist.computeOwners( blockOwners, mLocal2Global );

    // bucketSort of the owners

    HArray<IndexType> perm;
    HArray<IndexType> offsets;

    HArrayUtils::bucketSortOffsets( offsets, perm, blockOwners, np );

    HArray<IndexType> sendIndexes;

    HArrayUtils::gather( sendIndexes, mLocal2Global, perm, common::BinaryOp::COPY );

    auto sendPlan = CommunicationPlan::buildByOffsets( hostReadAccess( offsets ).get(), np );
    auto recvPlan = comm.transpose( sendPlan );

    SCAI_LOG_DEBUG( logger, comm << ": query owners send plan: " << sendPlan )
    SCAI_LOG_DEBUG( logger, comm << ": query owners recv plan: " << recvPlan )

    HArray<IndexType> receivedIndexes;

    comm.exchangeByPlan( receivedIndexes, recvPlan, sendIndexes, sendPlan );

    {
        auto wLocalOwners = hmemo::hostWriteAccess( mBlockDistributedOwners );
        auto rGlobalIndexes = hmemo::hostReadAccess( receivedIndexes );

        for ( PartitionId k = 0; k < recvPlan.size(); k++ )
        {
            CommunicationPlan::Entry entry = recvPlan[k];

            for ( IndexType i = 0; i < entry.quantity; ++i )
            {
                IndexType localIndex = rGlobalIndexes[entry.offset + i] - myLB;
                SCAI_ASSERT_VALID_INDEX_DEBUG( localIndex, localBlockSize, "serious problem" )
                wLocalOwners[localIndex] = entry.partitionId;
            }
        }
    }

    // Sanitize check there must be no invalidPartition value for any global index

    bool okay = HArrayUtils::allScalar( mBlockDistributedOwners, common::CompareOp::NE, invalidPartition );

    okay = comm.all( okay );

    if ( !okay )
    {
        COMMON_THROWEXCEPTION( "Illegal general distribution, at least one index does not appear" )
    }
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::computeOwners(
    hmemo::HArray<PartitionId>& owners, 
    const hmemo::HArray<IndexType>& indexes ) const
{
    const bool localDebug = false;

    SCAI_REGION( "Distribution.General.computeOwners" )

    HArrayUtils::setSameValue( owners, indexes.size(), invalidPartition, Context::getHostPtr() );

    BlockDistribution blockDist( getGlobalSize(), getCommunicatorPtr() );

    const Communicator& comm = getCommunicator();
    const PartitionId   np   = comm.getSize();

    HArray<PartitionId> blockIndexOwners;

    // get for each queried index the processor that knows the owner

    blockDist.computeOwners( blockIndexOwners, indexes );

    // bucketSort of the owners

    HArray<IndexType> perm;
    HArray<IndexType> sizes;

    HArrayUtils::bucketSortSizes( sizes, perm, blockIndexOwners, np );

    HArray<IndexType> sendIndexes;

    HArrayUtils::gather( sendIndexes, indexes, perm, common::BinaryOp::COPY );

    if ( localDebug )
    {
        auto rIndexes = hostReadAccess( sendIndexes );
        for ( IndexType i = 0; i < sendIndexes.size(); ++i )
        {
            std::cout << comm << ": queryIndex[ " << i << " ] = " << rIndexes[i] << std::endl;
        }
    }

    CommunicationPlan sendPlan( hostReadAccess( sizes ) );
    auto recvPlan = comm.transpose( sendPlan );

    SCAI_LOG_DEBUG( logger, comm << ": send plan: " << sendPlan )
    SCAI_LOG_DEBUG( logger, comm << ": recv plan: " << recvPlan )

    HArray<IndexType> receivedIndexes;

    comm.exchangeByPlan( receivedIndexes, recvPlan, sendIndexes, sendPlan );

    if ( localDebug )
    {
        auto rIndexes = hostReadAccess( receivedIndexes );
        for ( IndexType i = 0; i < receivedIndexes.size(); ++i )
        {
            std::cout << comm << ": queriedIndex[ " << i << " ] = " << rIndexes[i] << std::endl;
        }
    }

    HArrayUtils::compute( receivedIndexes, receivedIndexes, common::BinaryOp::SUB, blockDist.lb() );

    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( receivedIndexes, blockDist.getLocalSize() ), "serious wrong index" )

    HArray<PartitionId> sendOwners;  // will take the owners of the queried indexes

    HArrayUtils::gather( sendOwners, mBlockDistributedOwners, receivedIndexes, common::BinaryOp::COPY );

    if ( localDebug )
    {
        auto rIndexes = hostReadAccess( receivedIndexes );
        auto rOwners  = hostReadAccess( sendOwners );

        for ( IndexType i = 0; i < receivedIndexes.size(); ++i )
        {
            std::cout << comm << ": queriedIndex[ " << i << " ] = " << rIndexes[i] 
                      << " has owner " << rOwners[i] << std::endl;
        }
    }

    HArray<PartitionId> recvOwners;  // will receive the owners for the queried indexes

    comm.exchangeByPlan( recvOwners, sendPlan, sendOwners, recvPlan );

    if ( localDebug )
    {
        auto rIndexes = hostReadAccess( sendIndexes );
        auto rOwners  = hostReadAccess( recvOwners );
        auto rPerm    = hostReadAccess( perm );
        for ( IndexType i = 0; i < sendIndexes.size(); ++i )
        {
            std::cout << comm << ": queriedIndex[ " << i << " ] = " << rIndexes[i] 
                      << " from pos " << rPerm[i] 
                      << " has owner " << rOwners[i] << std::endl;
        }
    }

    // now set the owner for each queried global index

    HArrayUtils::scatter( owners, perm, true, recvOwners, common::BinaryOp::COPY );
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
