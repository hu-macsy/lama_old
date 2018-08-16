/**
 * @file HaloBuilder.cpp
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
 * @brief HaloBuilder.cpp
 * @author brandes
 * @date 24.03.2011
 */

// hpp
#include <scai/dmemo/HaloBuilder.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>

#include <vector>

using namespace scai::hmemo;

using scai::utilskernel::HArrayUtils;

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( HaloBuilder::logger, "Halo.Builder" )

void HaloBuilder::build( const Distribution& distribution, const HArray<IndexType>& requiredIndexes, Halo& halo )
{
    SCAI_REGION( "HaloBuilder.build" )
    halo.clear();

    const PartitionId noPartitions = distribution.getNumPartitions();
    const Communicator& communicator = distribution.getCommunicator();
    SCAI_LOG_INFO( logger,
                   communicator << ": building halo for " << noPartitions << " partitions, # requiredIndexes = " << requiredIndexes.size() )
    IndexType numIndexes = requiredIndexes.size();
    HArray<PartitionId> owners;

    {
        SCAI_REGION( "HaloBuilder.computeOwners" )
        distribution.computeOwners( owners, requiredIndexes );
    }

#ifdef SCAI_LOG_TRACE
    {
        ReadAccess<IndexType> rIndexes( requiredIndexes );
        ReadAccess<PartitionId> rOwners( owners );

        for ( IndexType i = 0; i < numIndexes; ++i )
        {
            SCAI_LOG_TRACE( logger, "Index " << rIndexes[i] << " belongs to " << rOwners[i] )
        }
    }
#endif

    CommunicationPlan& requiredPlan = halo.mRequiredPlan;

    //allocate Required plan with the nodes, where we get data from

    {
        ReadAccess<PartitionId> rOwners( owners );
        requiredPlan.allocateByOwners( noPartitions, rOwners.get(), owners.size() );
    }

    SCAI_LOG_INFO( logger,
                   communicator << ": allocated required plan for " << noPartitions << " partitions, size = " << requiredPlan.size() << ", total quantity = " << requiredPlan.totalQuantity() )

    // sort required indexes by the owner and define global->local mapping

    std::vector<IndexType> counts( noPartitions, invalidIndex );  // initialize with illegal index

    for ( PartitionId p = 0; p < requiredPlan.size(); ++p )
    {
        SCAI_LOG_TRACE( logger,
                        "requiredPlan[ " << p << "]: offset = " << requiredPlan[p].offset << ", pid = " << requiredPlan[p].partitionId << ", quantity = " << requiredPlan[p].quantity )
        counts[requiredPlan[p].partitionId] = requiredPlan[p].offset;
    }

    ContextPtr contextPtr = Context::getHostPtr();
    // allocate array for the required indexes sorted by owner
    WriteAccess<IndexType> requiredIndexesByOwner( halo.mRequiredIndexes, contextPtr );
    ReadAccess<PartitionId> rOwners( owners, contextPtr );
    ReadAccess<IndexType> rIndexes( requiredIndexes, contextPtr );
    requiredIndexesByOwner.resize( requiredPlan.totalQuantity() );

    // Note: size of requiredIndexesByOwner == requiredPlan.totalQuanitity()

    for ( IndexType jj = 0; jj < requiredPlan.totalQuantity(); ++jj )
    {
        PartitionId owner = rOwners[jj];
        SCAI_ASSERT( owner != invalidPartition, "No owner for required Index " << rIndexes[jj] )
        //The offset for the owner
        IndexType haloIndex = counts[owner];
        SCAI_ASSERT( haloIndex != invalidIndex, "No offset for owner " << owner << " of required Index " << rIndexes[jj] )
        requiredIndexesByOwner[haloIndex] = rIndexes[jj];
        // define the corresponding global->local mapping
        halo.setGlobal2Halo( rIndexes[jj], haloIndex );
        ++counts[owner];
    }

    requiredIndexesByOwner.release();
    // requiredPlan is ready, now we build the provides plan
    CommunicationPlan& providesPlan = halo.mProvidesPlan;
    providesPlan = requiredPlan.transpose( communicator );
    SCAI_LOG_DEBUG( logger, communicator << ": providesPlan = " << providesPlan )
    SCAI_LOG_TRACE( logger, "requiredPlan: " << requiredPlan )
    SCAI_LOG_TRACE( logger, "providesPlan: " << providesPlan )
    // communicate required indexes to other processors to get provideIndexes
    communicator.exchangeByPlan( halo.mProvidesIndexes, providesPlan, halo.mRequiredIndexes, requiredPlan );
    SCAI_LOG_INFO( logger, "exchanged plan indexes" )
#ifdef SCAI_LOG_TRACE
    {
        ReadAccess<IndexType> provide( halo.mProvidesIndexes, contextPtr );
        ReadAccess<IndexType> required( halo.mRequiredIndexes, contextPtr );

        for ( IndexType i = 0; i < provide.size(); ++i )
        {
            SCAI_LOG_TRACE( logger, "halo.mProvidesIndexes[" << i << "] " << provide[i] )
        }

        for ( IndexType i = 0; i < required.size(); ++i )
        {
            SCAI_LOG_TRACE( logger, "halo.mRequiredIndexes[" << i << "] " << required[i] )
        }
    }
#endif
    // localize the provides indexes that are still global from other processors
    WriteAccess<IndexType> providesIndexes( halo.mProvidesIndexes, contextPtr );

    for ( PartitionId p = 0; p < providesPlan.size(); ++p )
    {
        IndexType n = providesPlan[p].quantity;
        SCAI_LOG_TRACE( logger,
                        "Partition " << providesPlan[p].partitionId << " needs " << n << " entries from me(Partition: " << communicator.getRank() << "), offset = " << providesPlan[p].offset )
        IndexType* partitionumIndexes = providesIndexes.get() + providesPlan[p].offset;

        for ( IndexType i = 0; i < n; i++ )
        {
            IndexType localIndex = distribution.global2local( partitionumIndexes[i] );
            SCAI_ASSERT( localIndex != invalidIndex,
                         "global index " << partitionumIndexes[i] << " is not local on Rank " << communicator.getRank() )
            partitionumIndexes[i] = localIndex;
        }
    }
}

static void createRequiredGlobal2HaloMapping( std::map<IndexType, IndexType> & map, const HArray<IndexType> & requiredIndexes )
{
    ReadAccess<IndexType> rRequiredIndexes( requiredIndexes );

    for ( IndexType i = 0; i < rRequiredIndexes.size(); ++i )
    {
        const auto globalIndex = rRequiredIndexes[i];
        const auto haloIndex = i;
        map[globalIndex] = haloIndex;
    }
}

static HArray<IndexType> globalizeProvidedIndexes( const HArray<IndexType> & haloProvidedIndexes,
        const HArray<IndexType> & halo2global )
{
    HArray<IndexType> globalProvidedIndexes;
    HArrayUtils::gather( globalProvidedIndexes, halo2global, haloProvidedIndexes, common::BinaryOp::COPY );
    return globalProvidedIndexes;
}

void HaloBuilder::buildFromProvidedOwners( const Communicator& comm,
        const HArray<IndexType>& halo2global,
        const HArray<PartitionId>& ownersOfProvided,
        Halo& halo )
{
    SCAI_ASSERT_EQUAL_ERROR( halo2global.size(), ownersOfProvided.size() );
    SCAI_REGION( "HaloBuilder.buildFromProvidedOwners" )
    halo.clear();

    // TODO: Make context an argument (with default to Host/default context)
    const auto contextPtr = Context::getContextPtr();
    const auto numPartitions = comm.getSize();

    auto& requiredPlan = halo.mRequiredPlan;
    auto& providedPlan = halo.mProvidesPlan;
    auto& requiredIndexes = halo.mRequiredIndexes;
    auto& providedIndexes = halo.mProvidesIndexes;

    hmemo::HArray<IndexType> offsets( numPartitions + 1, contextPtr );
    providedIndexes.resize( ownersOfProvided.size() );

    HArrayUtils::bucketSort( offsets, providedIndexes, ownersOfProvided, numPartitions, contextPtr );

    {
        const auto rOffsets = hostReadAccess( offsets );
        providedPlan.allocateByOffsets( rOffsets.get(), numPartitions );
    }

    requiredPlan = providedPlan.transpose( comm );
    requiredIndexes.resize( requiredPlan.totalQuantity() );

    const auto globalProvidedIndexes = globalizeProvidedIndexes( providedIndexes, halo2global );
    comm.exchangeByPlan( requiredIndexes, requiredPlan, globalProvidedIndexes, providedPlan );

    createRequiredGlobal2HaloMapping( halo.mGlobal2Halo, requiredIndexes );
}

} /* end namespace dmemo */

} /* end namespace scai */
