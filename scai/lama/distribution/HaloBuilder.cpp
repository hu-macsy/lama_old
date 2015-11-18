/**
 * @file HaloBuilder.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @since 1.0.0
 */

// hpp
#include <scai/lama/distribution/HaloBuilder.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>

using namespace scai::hmemo;

namespace scai
{

namespace lama
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

    for( unsigned int i = 0; i < requiredIndexes.size(); ++i )
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

    for( IndexType p = 0; p < requiredPlan.size(); ++p )
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

    for( IndexType jj = 0; jj < requiredPlan.totalQuantity(); ++jj )
    {
        PartitionId owner = owners[jj];

        SCAI_ASSERT( owner >= 0 && owner < noPartitions, "No owner for required Index " << requiredIndexes[jj] )

        //The offset for the owner
        IndexType haloIndex = counts[owner];

        SCAI_ASSERT( haloIndex >= 0, "No offset for owner "<<owner<< " of required Index "<<requiredIndexes[jj] )

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

        for( int i = 0; i < provide.size(); ++i )
        {
            SCAI_LOG_TRACE( logger, "halo.mProvidesIndexes[" << i << "] " << provide[i] )
        }

        for( int i = 0; i < required.size(); ++i )
        {
            SCAI_LOG_TRACE( logger, "halo.mRequiredIndexes[" << i << "] " << required[i] )
        }
    }
#endif

    // localize the provides indexes that are still global from other processors

    WriteAccess<IndexType> providesIndexes( halo.mProvidesIndexes, contextPtr );

    for( PartitionId p = 0; p < providesPlan.size(); ++p )
    {
        IndexType n = providesPlan[p].quantity;

        SCAI_LOG_TRACE( logger,
                        "Partition " << providesPlan[p].partitionId << " needs " << n << " entries from me(Partition: " << communicator.getRank() << "), offset = " << providesPlan[p].offset )

        IndexType* partitionIndexes = providesIndexes.get() + providesPlan[p].offset;

        for( IndexType i = 0; i < n; i++ )
        {
            IndexType localIndex = distribution.global2local( partitionIndexes[i] );
            SCAI_ASSERT( localIndex != nIndex,
                         "global index "<<partitionIndexes[i]<<" is not local on Rank " << communicator.getRank() )
            partitionIndexes[i] = localIndex;
        }
    }
}

} /* end namespace lama */

} /* end namespace scai */
