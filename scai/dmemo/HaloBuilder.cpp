/**
 * @file HaloBuilder.cpp
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
 * @brief HaloBuilder.cpp
 * @author brandes
 * @date 24.03.2011
 */

// hpp
#include <scai/dmemo/HaloBuilder.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>

#include <set>

using namespace scai::hmemo;

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( HaloBuilder::logger, "Halo.Builder" )

void HaloBuilder::build( const Distribution& distribution, const HArray<IndexType>& requiredIndexes, Halo& halo )
{
    SCAI_REGION( "HaloBuilder.build" )
    const PartitionId noPartitions = distribution.getNumPartitions();
    const Communicator& communicator = distribution.getCommunicator();
    SCAI_LOG_INFO( logger,
                   communicator << ": building halo for " << noPartitions << " partitions, # requiredIndexes = " << requiredIndexes.size() )
    IndexType nIndexes = requiredIndexes.size();
    HArray<PartitionId> owners;

    {
        SCAI_REGION( "HaloBuilder.computeOwners" )
        distribution.computeOwners( owners, requiredIndexes );
    }

#ifdef SCAI_LOG_TRACE
    {
        ReadAccess<IndexType> rIndexes( requiredIndexes );
        ReadAccess<PartitionId> rOwners( owners );

        for ( IndexType i = 0; i < nIndexes; ++i )
        {
            SCAI_LOG_TRACE( logger, "Index " << rIndexes[i] << " belongs to " << rOwners[i] )
        }
    }
#endif

    CommunicationPlan& requiredPlan = halo.mRequiredPlan;

    //allocate Required plan with the nodes, where we get data from

    {
        ReadAccess<PartitionId> rOwners( owners );
        requiredPlan.allocate( noPartitions, rOwners.get(), owners.size() );
    }

    SCAI_LOG_INFO( logger,
                   communicator << ": allocated required plan for " << noPartitions << " partitions, size = " << requiredPlan.size() << ", total quantity = " << requiredPlan.totalQuantity() )

    // sort required indexes by the owner and define global->local mapping

    std::vector<IndexType> counts( noPartitions, nIndex );  // initialize with illegal index

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
        SCAI_ASSERT( owner != nPartition, "No owner for required Index " << rIndexes[jj] )
        //The offset for the owner
        IndexType haloIndex = counts[owner];
        SCAI_ASSERT( haloIndex != nIndex, "No offset for owner " << owner << " of required Index " << rIndexes[jj] )
        requiredIndexesByOwner[haloIndex] = rIndexes[jj];
        // define the corresponding global->local mapping
        halo.setGlobal2Halo( rIndexes[jj], haloIndex );
        ++counts[owner];
    }

    requiredIndexesByOwner.release();
    // requiredPlan is ready, now we build the provides plan
    CommunicationPlan& providesPlan = halo.mProvidesPlan;
    providesPlan.allocateTranspose( requiredPlan, communicator );
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
        IndexType* partitionIndexes = providesIndexes.get() + providesPlan[p].offset;

        for ( IndexType i = 0; i < n; i++ )
        {
            IndexType localIndex = distribution.global2local( partitionIndexes[i] );
            SCAI_ASSERT( localIndex != nIndex,
                         "global index " << partitionIndexes[i] << " is not local on Rank " << communicator.getRank() )
            partitionIndexes[i] = localIndex;
        }
    }
}

void HaloBuilder::build( const Distribution& distribution, const HArray<IndexType>& requiredIndexes, 
    const HArray<IndexType>& providedIndexes, const PartitionId partner, Halo& halo)
{
    SCAI_REGION( "HaloBuilder.buildWithPartner" )
    const PartitionId noPartitions = distribution.getNumPartitions();
    ContextPtr contextPtr = Context::getHostPtr();

    //skip computation of owners

    CommunicationPlan& requiredPlan = halo.mRequiredPlan;
    IndexType nIndexes = requiredIndexes.size();
    IndexType nProvidedIndexes = providedIndexes.size();
    ReadAccess<IndexType> rIndexes( requiredIndexes, contextPtr );

    HArray<PartitionId> owners;
    owners.init(partner, nIndexes);

    //allocate Required plan with the nodes, where we get data from
    {
        ReadAccess<PartitionId> rOwners( owners );
        requiredPlan.allocate( noPartitions, rOwners.get(), owners.size() );
    }
    
    std::vector<IndexType> counts( noPartitions, nIndex );  // initialize with illegal index
    for ( PartitionId p = 0; p < requiredPlan.size(); ++p )
    {
        counts[requiredPlan[p].partitionId] = requiredPlan[p].offset;
    }

    WriteAccess<IndexType> requiredIndexesByOwner( halo.mRequiredIndexes, contextPtr );
    requiredIndexesByOwner.resize( requiredPlan.totalQuantity() );

    for ( IndexType jj = 0; jj < requiredPlan.totalQuantity(); ++jj )
    {
        IndexType haloIndex = counts[partner];
        requiredIndexesByOwner[haloIndex] = rIndexes[jj];
        // define the corresponding global->local mapping
        halo.setGlobal2Halo( rIndexes[jj], haloIndex );
        ++counts[partner];
    }
    requiredIndexesByOwner.release();

    CommunicationPlan& providesPlan = halo.mProvidesPlan;
    HArray<PartitionId> targets;
    targets.init(partner, nProvidedIndexes);

    //allocate Required plan with the nodes, where we get data from
    {
        ReadAccess<PartitionId> rTargets( targets );
        providesPlan.allocate( noPartitions, rTargets.get(), targets.size() );
    }

    SCAI_ASSERT( providesPlan.totalQuantity() == nProvidedIndexes,
        "Provides plan has " << providesPlan.totalQuantity() << " entries, but " << nProvidedIndexes << " were given." )

    WriteAccess<IndexType> providesIndexesAccess( halo.mProvidesIndexes, contextPtr );
    ReadAccess<IndexType> pIndexes( providedIndexes, contextPtr );

    providesIndexesAccess.resize( providesPlan.totalQuantity() );

    for (IndexType i = 0; i < nProvidedIndexes; i++) {
        providesIndexesAccess[i] = distribution.global2local( pIndexes[i] );
    }
}

void HaloBuilder::coarsenHalo(const Distribution& coarseDistribution, const Halo& halo, const scai::hmemo::HArray<IndexType>& localFineToCoarse, const scai::hmemo::HArray<IndexType>& haloFineToCoarse, Halo& coarseHalo) {
    SCAI_REGION( "HaloBuilder.coarsenHalo" )
    scai::hmemo::ReadAccess<IndexType> providedIndices(halo.getProvidesIndexes());
    scai::hmemo::ReadAccess<IndexType> requiredIndices(halo.getRequiredIndexes());
    scai::dmemo::CommunicationPlan sendPlan = halo.getProvidesPlan();
    scai::dmemo::CommunicationPlan recvPlan = halo.getRequiredPlan();

    SCAI_ASSERT(providedIndices.size() == sendPlan.totalQuantity(), "Communication plan does not fit provided indices.");
    SCAI_ASSERT(requiredIndices.size() == recvPlan.totalQuantity(), "Communication plan does not fit required indices.");

    std::vector<IndexType> newProvidedIndices;
    std::vector<IndexType> sendQuantities;

    {
        scai::hmemo::ReadAccess<IndexType> rFineToCoarse(localFineToCoarse);
        //construct new send plan
        for (IndexType i = 0; i < sendPlan.size(); i++) {
            scai::dmemo::CommunicationPlan::Entry entry = sendPlan[i];
            if (IndexType(sendQuantities.size()) <= entry.partitionId) sendQuantities.resize(entry.partitionId+1);
            std::set<IndexType> sendSet;
            for (IndexType j = entry.offset; j < entry.offset + entry.quantity; j++) {
                SCAI_ASSERT(j < providedIndices.size(), "Communication plan does not fit provided indices.");
                IndexType provIndex = providedIndices[j];
                SCAI_ASSERT(provIndex < rFineToCoarse.size(), "Provided index " << provIndex << " seemingly not local.");
                sendSet.insert(coarseDistribution.global2local(rFineToCoarse[providedIndices[j]]));
            }
            newProvidedIndices.insert(newProvidedIndices.end(), sendSet.begin(), sendSet.end());
            sendQuantities[entry.partitionId] = sendSet.size();
        }
    }
    SCAI_ASSERT(IndexType(newProvidedIndices.size()) <= providedIndices.size(), "New index list is bigger than old one.");
    coarseHalo.mProvidesPlan.allocate(sendQuantities.data(), sendQuantities.size());

    SCAI_ASSERT( coarseHalo.mProvidesPlan.totalQuantity() == IndexType(newProvidedIndices.size()),
        "Send plan has " << coarseHalo.mProvidesPlan.totalQuantity() << " entries, but " << newProvidedIndices.size() << " were given." )

    std::vector<IndexType> newRequiredIndices;
    std::vector<IndexType> recvQuantities;

    {
        scai::hmemo::ReadAccess<IndexType> rFineToCoarse(haloFineToCoarse);
        //construct new recv plan
        for (IndexType i = 0; i < recvPlan.size(); i++) {
            scai::dmemo::CommunicationPlan::Entry entry = recvPlan[i];
            if (IndexType(recvQuantities.size()) <= entry.partitionId) recvQuantities.resize(entry.partitionId+1);
            std::set<IndexType> recvSet;
            for (IndexType j = entry.offset; j < entry.offset + entry.quantity; j++) {
                IndexType reqIndex = requiredIndices[j];
                SCAI_ASSERT(halo.global2halo(reqIndex) != nIndex, "Index" << reqIndex << " seemingly not in halo");
                SCAI_ASSERT(halo.global2halo(reqIndex) < rFineToCoarse.size(), "Index" << halo.global2halo(reqIndex) << " too big for halo data");
                recvSet.insert(rFineToCoarse[halo.global2halo(requiredIndices[j])]);
            }
            for (IndexType reqIndex : recvSet) {
                newRequiredIndices.push_back(reqIndex);
            }
            recvQuantities[entry.partitionId] = recvSet.size();
        }
    }
    SCAI_ASSERT(IndexType(newRequiredIndices.size()) <= requiredIndices.size(), "New index list is bigger than old one.");

    coarseHalo.mRequiredPlan.allocate(recvQuantities.data(), recvQuantities.size());

    SCAI_ASSERT( coarseHalo.mRequiredPlan.totalQuantity() == IndexType(newRequiredIndices.size()),
        "Provides plan has " << coarseHalo.mRequiredPlan.totalQuantity() << " entries, but " << newRequiredIndices.size() << " were given." )
    
    coarseHalo.mRequiredIndexes = HArray<IndexType>(newRequiredIndices.size(), newRequiredIndices.data());
    coarseHalo.mProvidesIndexes = HArray<IndexType>(newProvidedIndices.size(), newProvidedIndices.data());
    scai::hmemo::ReadAccess<IndexType> rRequired(coarseHalo.mRequiredIndexes);
    for (IndexType i = 0; i < rRequired.size(); i++) {
        coarseHalo.setGlobal2Halo(rRequired[i], i);
    }
}


} /* end namespace dmemo */

} /* end namespace scai */
