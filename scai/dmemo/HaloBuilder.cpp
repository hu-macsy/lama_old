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

using scai::utilskernel::HArrayUtils;

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( HaloBuilder::logger, "Halo.Builder" )

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

void HaloBuilder::buildFromRequired( Halo& halo, const Distribution& distribution, const HArray<IndexType>& requiredIndexes )
{
    SCAI_REGION( "HaloBuilder.build" )
    halo.clear();

    const Communicator& communicator = distribution.getCommunicator();
    const PartitionId NP = communicator.getSize();

    SCAI_LOG_INFO( logger,
                   communicator << ": building halo for " << NP << " processors, # requiredIndexes = " << requiredIndexes.size() )

    HArray<PartitionId> owners;

    {
        SCAI_REGION( "HaloBuilder.computeOwners" )
        distribution.computeOwners( owners, requiredIndexes );
    }

    // global exchange of required indexes similiar to globalExchange, but here we keep the communication plans

    // sort the required indexes by the owners via bucket sort, bucket sizes will be used for communication plan

    hmemo::HArray<IndexType> perm;
    hmemo::HArray<IndexType> sizes;

    utilskernel::HArrayUtils::bucketSortSizes( sizes, perm, owners, NP );

    // bucket sorted required indexes become part of the halo

    utilskernel::HArrayUtils::gather( halo.mRequiredIndexes, requiredIndexes, perm, common::BinaryOp::COPY );

    // global indexes -> halo indexes mapping is also essential for building halo

    createRequiredGlobal2HaloMapping( halo.mGlobal2Halo, halo.mRequiredIndexes );

    halo.mRequiredPlan = CommunicationPlan( hostReadAccess( sizes ) );

    // requiredPlan is ready, now we build the provides plan
    halo.mProvidesPlan = communicator.transpose( halo.mRequiredPlan );

    // exchange required indexes (sorted in buckets for each processor)  with other processors to get provideIndexes

    communicator.exchangeByPlan( halo.mProvidesIndexes, halo.mProvidesPlan, halo.mRequiredIndexes, halo.mRequiredPlan );

    SCAI_LOG_INFO( logger, "exchanged non-local indexes" )

    // localize the provides indexes that are still global from other processors

    distribution.global2localV( halo.mProvidesIndexes, halo.mProvidesIndexes );
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

    hmemo::HArray<IndexType> sizes;

    providedIndexes.resize( ownersOfProvided.size() );

    HArrayUtils::bucketSortSizes( sizes, providedIndexes, ownersOfProvided, numPartitions, contextPtr );

    providedPlan = CommunicationPlan( hostReadAccess( sizes ) );
    requiredPlan = comm.transpose( providedPlan );

    requiredIndexes.resize( requiredPlan.totalQuantity() );

    const auto globalProvidedIndexes = globalizeProvidedIndexes( providedIndexes, halo2global );
    comm.exchangeByPlan( requiredIndexes, requiredPlan, globalProvidedIndexes, providedPlan );

    createRequiredGlobal2HaloMapping( halo.mGlobal2Halo, requiredIndexes );
}

} /* end namespace dmemo */

} /* end namespace scai */
