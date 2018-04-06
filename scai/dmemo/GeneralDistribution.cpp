/**
 * @file GeneralDistribution.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief GeneralDistribution.cpp
 * @author brandes
 * @date 25.02.2011
 */

// hpp
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

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
    const HArray<IndexType>& myIndexes,
    const CommunicatorPtr communicator ) :

    Distribution( globalSize, communicator )

{
    SCAI_LOG_INFO( logger, "GeneralDistribution( size = " << globalSize
                   << ", myIndexes = " << myIndexes.size() << ", comm = " << *communicator )

    SCAI_ASSERT_DEBUG( HArrayUtils::validIndexes( myIndexes, globalSize ), "myIndexes contains illegal values" )

    // make a copy of the indexes and sort them ascending = true, perm not needed

    HArrayUtils::sort( NULL, &mLocal2Global, myIndexes, true );

    // Note: the constructor is completely local, but make some consistency check now

    IndexType nLocal = myIndexes.size();

    SCAI_ASSERT_EQ_ERROR( mGlobalSize, communicator->sum( nLocal ), "illegal general distribution" )

    fillIndexMap();
}

/* ---------------------------------------------------------------------- */

GeneralDistribution::GeneralDistribution(
    const HArray<PartitionId>& owners,
    const CommunicatorPtr communicator ) :

    Distribution( 0, communicator )

{
    using namespace hmemo;

    ContextPtr ctx = Context::getHostPtr();

    PartitionId rank = mCommunicator->getRank();
    PartitionId size = mCommunicator->getSize();

    if ( rank == MASTER )
    {
        mGlobalSize = owners.size();
    }

    mCommunicator->bcast( &mGlobalSize, 1, MASTER );

    SCAI_LOG_DEBUG( logger, *mCommunicator << ": global size = " << mGlobalSize )

    HArray<IndexType> localSizes;
    HArray<IndexType> localOffsets;

    // count in localSizes for each partition the owners
    // owners = [ 0, 1, 2, 1, 2, 1, 0 ] -> sizes = [2, 3, 2 ]
    // sizes.sum() == owners.size()

    if ( rank == MASTER )
    {
        HArrayUtils::bucketCount( localSizes, owners, size );
        IndexType lsum = HArrayUtils::sum( localSizes );
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": sum( localSizes ) = " << lsum << ", must be " << mGlobalSize );
    }
    else
    {
        // all other procs intialize localSizes with a dummy value to avoid read access of uninitialized array
        HArrayUtils::setOrder( localSizes, 1 );
    }

    IndexType localSize;

    // scatter partition sizes

    {
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": before scatter, localSizes = " << localSizes )
        ReadAccess<IndexType> rSizes( localSizes );
        mCommunicator->scatter( &localSize, 1, MASTER, rSizes.get() );
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": after scatter, localSize = " << localSize )
    }

    SCAI_LOG_DEBUG( logger, *mCommunicator << ": owns " << localSize << " of " << mGlobalSize << " elements" )

    if ( rank == MASTER )
    {
        localOffsets.reserve( ctx, size + 1 );
        localOffsets = localSizes;
        HArrayUtils::scan1( localOffsets );
        SCAI_LOG_DEBUG( logger, "scan done, sum = " << localOffsets[ size ] )
    }

    // Now resort 0, ..., n - 1 according to the owners

    HArray<IndexType> sortedIndexes;

    if ( rank == MASTER )
    {
        SCAI_LOG_DEBUG( logger, "reorder for indexes" )

        ContextPtr loc = Context::getHostPtr();

        static LAMAKernel<UtilKernelTrait::sortInBuckets<PartitionId> > sortInBuckets;

        sortInBuckets.getSupportedContext( loc );

        WriteOnlyAccess<IndexType> wIndexes( sortedIndexes, loc, mGlobalSize );
        WriteAccess<IndexType> wOffsets( localOffsets, loc );
        ReadAccess<PartitionId> rOwners( owners, loc );

        sortInBuckets[loc]( wIndexes, wOffsets, size, rOwners, mGlobalSize );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": initialize sortedIndexes " )

        // all other procs intialize sortedIndexes with a dummy value to avoid read access of uninitialized array
        HArrayUtils::setOrder( sortedIndexes, 1 );
    }

    WriteOnlyAccess<IndexType> wLocal2Global( mLocal2Global, localSize );

    {
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": before scatterV, sortedIndexes = " << sortedIndexes  )
        ReadAccess<IndexType> rIndexes( sortedIndexes );
        ReadAccess<IndexType> rSizes( localSizes );
        mCommunicator->scatterV( wLocal2Global.get(), localSize, MASTER, rIndexes.get(), rSizes.get() );
        SCAI_LOG_DEBUG( logger, *mCommunicator << ": after scatterV, sortedIndexes = " << sortedIndexes )
    }

    wLocal2Global.release();

    fillIndexMap();
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

        return;
    }

    // make a bucket sort with owners

    HArray<IndexType> offsets;    // number of elements for each bucket/partition
    HArray<IndexType> perm;       // permutation to resort in buckets

    HArrayUtils::bucketSort( offsets, perm, owners, mCommunicator->getSize() );

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
    auto recvPlan = sendPlan.transpose( comm );

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

    // add local indexes into unordered_map

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

    GenBlockDistribution genBlock( mGlobalSize, localSize, comm );

    IndexType lb;
    IndexType ub;

    genBlock.getLocalRange( lb, ub );

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

void GeneralDistribution::getBlockDistributedOwners( hmemo::HArray<PartitionId>& localOwners ) const
{
    // get my range for block distribution of size

    BlockDistribution blockDist( getGlobalSize(), getCommunicatorPtr() );

    const Communicator& comm = getCommunicator();
    const PartitionId   np   = comm.getSize();

    hmemo::ContextPtr hostCtx = hmemo::Context::getHostPtr();

    const IndexType localBlockSize = blockDist.getLocalSize();
    const IndexType myLB           = blockDist.lb();

    HArrayUtils::setSameValue( localOwners, localBlockSize, invalidPartition, hostCtx );

    // get Owners of myIndexes according to the block distribution

    HArray<PartitionId> blockOwners;
    blockDist.computeOwners( blockOwners, mLocal2Global );

    // bucketSort of the owners

    HArray<IndexType> perm;
    HArray<IndexType> offsets;

    HArrayUtils::bucketSort( offsets, perm, blockOwners, np );

    HArray<IndexType> sendIndexes;

    HArrayUtils::gather( sendIndexes, mLocal2Global, perm, common::BinaryOp::COPY );

    auto sendPlan = CommunicationPlan::buildByOffsets( hostReadAccess( offsets ).get(), np );
    auto recvPlan = sendPlan.transpose( comm );

    SCAI_LOG_ERROR( logger, comm << ": query owners send plan: " << sendPlan )
    SCAI_LOG_ERROR( logger, comm << ": query owners recv plan: " << recvPlan )

    HArray<IndexType> receivedIndexes;

    comm.exchangeByPlan( receivedIndexes, recvPlan, sendIndexes, sendPlan );

    {
        auto wLocalOwners = hmemo::hostWriteAccess( localOwners );
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
}

/* ---------------------------------------------------------------------- */

void GeneralDistribution::computeOwners(
    hmemo::HArray<PartitionId>& owners, 
    const hmemo::HArray<IndexType>& indexes ) const
{
    HArrayUtils::setSameValue( owners, indexes.size(), invalidPartition, Context::getHostPtr() );

    BlockDistribution blockDist( getGlobalSize(), getCommunicatorPtr() );

    const Communicator& comm = getCommunicator();
    const PartitionId   np   = comm.getSize();

    HArray<PartitionId> blockDistOwners;

    getBlockDistributedOwners( blockDistOwners );

    if ( true )
    {
        auto rOwners = hostReadAccess( blockDistOwners );

        IndexType lb = blockDist.lb();

        for ( IndexType i = 0; i < blockDistOwners.size(); ++i )
        {
            std::cout << comm << ": owner[ " << ( i + lb ) << " ] = " << rOwners[i] << std::endl;
        }
    }

    HArray<PartitionId> blockIndexOwners;

    // get the processors that know the owners of the queried indexes

    blockDist.computeOwners( blockIndexOwners, indexes );

    // bucketSort of the owners

    HArray<IndexType> perm;
    HArray<IndexType> offsets;

    HArrayUtils::bucketSort( offsets, perm, blockIndexOwners, np );

    HArray<IndexType> sendIndexes;

    HArrayUtils::gather( sendIndexes, indexes, perm, common::BinaryOp::COPY );

    if ( true )
    {
        auto rIndexes = hostReadAccess( sendIndexes );
        for ( IndexType i = 0; i < sendIndexes.size(); ++i )
        {
            std::cout << comm << ": queryIndex[ " << i << " ] = " << rIndexes[i] << std::endl;
        }
    }

    auto sendPlan = CommunicationPlan::buildByOffsets( hostReadAccess( offsets ).get(), np );
    auto recvPlan = sendPlan.transpose( comm );

    SCAI_LOG_DEBUG( logger, comm << ": send plan: " << sendPlan )
    SCAI_LOG_DEBUG( logger, comm << ": recv plan: " << recvPlan )

    HArray<IndexType> receivedIndexes;

    comm.exchangeByPlan( receivedIndexes, recvPlan, sendIndexes, sendPlan );

    if ( true )
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

    HArrayUtils::gather( sendOwners, blockDistOwners, receivedIndexes, common::BinaryOp::COPY );

    if ( true )
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

    if ( true )
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
        SCAI_ASSERT_EQ_DEBUG( mAllLocalOffsets.size(), getNumPartitions() + 1, "serious mismatch" )
        SCAI_ASSERT_EQ_DEBUG( mAllLocal2Global.size(), getGlobalSize(), "serious mismatch" )
        SCAI_ASSERT_EQ_DEBUG( mAllGlobal2Local.size(), getGlobalSize(), "serious mismatch" )

        return;   // already done
    }

    // compute mAllOwners

    HArray<IndexType> indexes;   // will contain all column indexes to get all owners

    HArrayUtils::setOrder( indexes, mGlobalSize );

    Distribution::computeOwners( mAllOwners, indexes );

    // bucket sort the owners, gives offsets and permutation to block values according to owners

    HArrayUtils::bucketSort( mAllLocalOffsets, mAllLocal2Global, mAllOwners, mCommunicator->getSize() );
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
