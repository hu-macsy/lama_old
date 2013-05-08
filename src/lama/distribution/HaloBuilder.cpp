/**
 * @file HaloBuilder.cpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief HaloBuilder.cpp
 * @author brandes
 * @date 24.03.2011
 * $Id$
 */

// hpp
#include <lama/distribution/HaloBuilder.hpp>

// others
#include <lama/exception/LAMAAssert.hpp>

// tracing
#include <lama/tracing.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( HaloBuilder::logger, "Halo.Builder" )

void HaloBuilder::build( const Distribution& distribution, const std::vector<IndexType>& requiredIndexes, Halo& halo )
{
    LAMA_REGION( "HaloBuilder.build" )

    const PartitionId noPartitions = distribution.getNumPartitions();

    const Communicator& communicator = distribution.getCommunicator();

    LAMA_LOG_INFO( logger,
                   communicator << ": building halo for " << noPartitions << " partitions, # requiredIndexes = " << requiredIndexes.size() )

    std::vector<PartitionId> owners;
    owners.reserve( requiredIndexes.size() );

    {
        LAMA_REGION( "HaloBuilder.computeOwners" )
        communicator.computeOwners( requiredIndexes, distribution, owners );
    }
#ifdef LAMA_LOG_TRACE
    for ( unsigned int i = 0; i < requiredIndexes.size(); ++i )
    {
        LAMA_LOG_TRACE( logger, "Index " << requiredIndexes[i] << " belongs to " << owners[i] )
    }
#endif

    CommunicationPlan& requiredPlan = halo.mRequiredPlan;

    //allocate Required plan with the nodes, where we get data from
    requiredPlan.allocate( noPartitions, owners );

    LAMA_LOG_INFO( logger,
                   communicator << ": allocated required plan for " << noPartitions << " partitions, size = " << requiredPlan.size() << ", total quantity = " << requiredPlan.totalQuantity() )

    // sort required indexes by the owner and define global->local mapping

    std::vector<IndexType> counts( noPartitions, -1 );

    for ( IndexType p = 0; p < requiredPlan.size(); ++p )
    {
        LAMA_LOG_TRACE( logger,
                        "requiredPlan[ " << p << "]: offset = " << requiredPlan[p].offset << ", pid = " << requiredPlan[p].partitionId << ", quantity = " << requiredPlan[p].quantity )
        counts[requiredPlan[p].partitionId] = requiredPlan[p].offset;
    }

    // allocate array for the required indexes sorted by owner

    HostWriteAccess<IndexType> requiredIndexesByOwner( halo.mRequiredIndexes );

    requiredIndexesByOwner.resize( requiredPlan.totalQuantity() );

    // Note: size of requiredIndexesByOnwer == requiredPlan.totalQuanitity()

    for ( IndexType jj = 0; jj < requiredPlan.totalQuantity(); ++jj )
    {
        PartitionId owner = owners[jj];

        LAMA_ASSERT( owner >= 0 && owner < noPartitions, "No owner for required Index " << requiredIndexes[jj] )

        //The offset for the owner
        IndexType haloIndex = counts[owner];

        LAMA_ASSERT( haloIndex >= 0, "No offset for owner "<<owner<< " of required Index "<<requiredIndexes[jj] )

        requiredIndexesByOwner[haloIndex] = requiredIndexes[jj];

        // define the corresponding global->local mapping

        halo.setGlobal2Halo( requiredIndexes[jj], haloIndex );

        ++counts[owner];
    }

    requiredIndexesByOwner.release();

    // requiredPlan is ready, now we build the provides plan

    CommunicationPlan& providesPlan = halo.mProvidesPlan;

    providesPlan.allocateTranspose( requiredPlan, communicator );

    LAMA_LOG_DEBUG( logger, communicator << ": providesPlan = " << providesPlan )

    LAMA_LOG_TRACE( logger, "requiredPlan: " << requiredPlan )
    LAMA_LOG_TRACE( logger, "providesPlan: " << providesPlan )

    // communicate required indexes to other processors to get provideIndexes

    communicator.exchangeByPlan( halo.mProvidesIndexes, providesPlan, halo.mRequiredIndexes, requiredPlan );

    LAMA_LOG_INFO( logger, "exchanged plan indexes" )
#ifdef LAMA_LOG_TRACE
    {
        HostReadAccess<IndexType> provide( halo.mProvidesIndexes );
        HostReadAccess<IndexType> required( halo.mRequiredIndexes );
        for ( int i = 0; i < provide.size(); ++i )
        {
            LAMA_LOG_TRACE( logger, "halo.mProvidesIndexes[" << i << "] " << provide[i] )
        }
        for ( int i = 0; i < required.size(); ++i )
        {
            LAMA_LOG_TRACE( logger, "halo.mRequiredIndexes[" << i << "] " << required[i] )
        }
    }
#endif

    // localize the provides indexes that are still global from other processors

    HostWriteAccess<IndexType> providesIndexes( halo.mProvidesIndexes );

    for ( PartitionId p = 0; p < providesPlan.size(); ++p )
    {
        IndexType n = providesPlan[p].quantity;

        LAMA_LOG_TRACE( logger,
                        "Partition " << providesPlan[p].partitionId << " needs " << n << " entries from me(Partition: " << communicator.getRank() << "), offset = " << providesPlan[p].offset )

        IndexType* partitionIndexes = providesIndexes.get() + providesPlan[p].offset;

        for ( IndexType i = 0; i < n; i++ )
        {
            IndexType localIndex = distribution.global2local( partitionIndexes[i] );
            LAMA_ASSERT( localIndex != nIndex,
                         "global index "<<partitionIndexes[i]<<" is not local on Rank " << communicator.getRank() )
            partitionIndexes[i] = localIndex;
        }
    }
}

} // namespace LAMA
