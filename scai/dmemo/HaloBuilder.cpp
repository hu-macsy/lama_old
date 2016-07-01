/**
 * @file HaloBuilder.cpp
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
 * @brief HaloBuilder.cpp
 * @author brandes
 * @date 24.03.2011
 */

// hpp
#include <scai/dmemo/HaloBuilder.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>

using namespace scai::hmemo;

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( HaloBuilder::logger, "Halo.Builder" )

void HaloBuilder::build( const Distribution& distribution, const std::vector<IndexType>& requiredIndexes, Halo& halo )
{
    SCAI_REGION( "HaloBuilder.build" )
    const PartitionId noPartitions = distribution.getNumPartitions();
    const Communicator& communicator = distribution.getCommunicator();
    SCAI_LOG_INFO( logger,
                   communicator << ": building halo for " << noPartitions << " partitions, # requiredIndexes = " << requiredIndexes.size() )
    IndexType nIndexes = requiredIndexes.size();
    std::vector<PartitionId> owners( nIndexes );
    {
        SCAI_REGION( "HaloBuilder.computeOwners" )
        communicator.computeOwners( &owners[0], distribution, &requiredIndexes[0], nIndexes );
    }
#ifdef SCAI_LOG_TRACE

    for ( unsigned int i = 0; i < requiredIndexes.size(); ++i )
    {
        SCAI_LOG_TRACE( logger, "Index " << requiredIndexes[i] << " belongs to " << owners[i] )
    }

#endif
    CommunicationPlan& requiredPlan = halo.mRequiredPlan;
    //allocate Required plan with the nodes, where we get data from
    requiredPlan.allocate( noPartitions, owners.data(), owners.size() );
    SCAI_LOG_INFO( logger,
                   communicator << ": allocated required plan for " << noPartitions << " partitions, size = " << requiredPlan.size() << ", total quantity = " << requiredPlan.totalQuantity() )
    // sort required indexes by the owner and define global->local mapping
    std::vector<IndexType> counts( noPartitions, -1 );

    for ( IndexType p = 0; p < requiredPlan.size(); ++p )
    {
        SCAI_LOG_TRACE( logger,
                        "requiredPlan[ " << p << "]: offset = " << requiredPlan[p].offset << ", pid = " << requiredPlan[p].partitionId << ", quantity = " << requiredPlan[p].quantity )
        counts[requiredPlan[p].partitionId] = requiredPlan[p].offset;
    }

    ContextPtr contextPtr = Context::getHostPtr();
    // allocate array for the required indexes sorted by owner
    WriteAccess<IndexType> requiredIndexesByOwner( halo.mRequiredIndexes, contextPtr );
    requiredIndexesByOwner.resize( requiredPlan.totalQuantity() );

    // Note: size of requiredIndexesByOwner == requiredPlan.totalQuanitity()

    for ( IndexType jj = 0; jj < requiredPlan.totalQuantity(); ++jj )
    {
        PartitionId owner = owners[jj];
        SCAI_ASSERT( owner >= 0 && owner < noPartitions, "No owner for required Index " << requiredIndexes[jj] )
        //The offset for the owner
        IndexType haloIndex = counts[owner];
        SCAI_ASSERT( haloIndex >= 0, "No offset for owner " << owner << " of required Index " << requiredIndexes[jj] )
        requiredIndexesByOwner[haloIndex] = requiredIndexes[jj];
        // define the corresponding global->local mapping
        halo.setGlobal2Halo( requiredIndexes[jj], haloIndex );
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

} /* end namespace dmemo */

} /* end namespace scai */
