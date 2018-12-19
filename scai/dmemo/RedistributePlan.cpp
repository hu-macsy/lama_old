/**
 * @file RedistributePlan.cpp
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
 * @brief Implementation of methods for RedistributePlan class.
 * @author Thomas Brandes, Andreas Longva
 * @date 08.10.2011
 */

// hpp
#include <scai/dmemo/RedistributePlan.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

// local library
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/hmemo/HostReadAccess.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostWriteOnlyAccess.hpp>

#include <memory>
#include <algorithm>

using namespace scai::hmemo;

using scai::utilskernel::HArrayUtils;
using scai::dmemo::GeneralDistribution;

namespace scai
{

using std::unique_ptr;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( RedistributePlan::logger, "RedistributePlan" )

static HArray<IndexType> ownedGlobalIndexesForDist( const Distribution& dist )
{
    HArray<IndexType> indexes;
    dist.getOwnedIndexes( indexes );
    return indexes;
}

// Partitions local indexes of the source distribution into "keep"
// and "exchange", based on the new owner of each individual index.
static void partitionSourceIndexes( HArray<IndexType> & keepLocalIndexes,
                                    HArray<IndexType> & exchangeLocalIndexes,
                                    const HArray<PartitionId> & newOwnersOfLocalElements,
                                    const Distribution& sourceDist )
{
    SCAI_ASSERT_EQ_DEBUG( newOwnersOfLocalElements.size(), sourceDist.getLocalSize(), "sourceDist and newOwners must have same size" );
    const auto rank = sourceDist.getCommunicator().getRank();
    const auto numPartitions = sourceDist.getCommunicator().getSize();
    const auto sourceNumLocal = sourceDist.getLocalSize();

    const auto rNewOwners = hostReadAccess( newOwnersOfLocalElements );
    auto wKeep = hostWriteOnlyAccess( keepLocalIndexes, sourceNumLocal );
    auto wExchange = hostWriteOnlyAccess( exchangeLocalIndexes, sourceNumLocal );

    IndexType numKeep = 0;
    IndexType numExchange = 0;

    for ( IndexType localSourceIndex = 0; localSourceIndex < sourceDist.getLocalSize(); ++localSourceIndex )
    {
        const auto newOwner = rNewOwners[localSourceIndex];
        SCAI_ASSERT_VALID_INDEX( newOwner, numPartitions, "owner index out of range" );

        if ( rank == newOwner )
        {
            wKeep[numKeep++] = localSourceIndex;
        }
        else
        {
            wExchange[numExchange++] = localSourceIndex;
        }
    }

    wKeep.resize( numKeep );
    wExchange.resize( numExchange );
}

template <typename ValueType>
static HArray<ValueType> gather( const HArray<ValueType> & source, const HArray<IndexType> & indexes )
{
    HArray<ValueType> result;
    HArrayUtils::gather( result, source, indexes, common::BinaryOp::COPY );
    return result;
}

static HArray<IndexType> local2global( const HArray<IndexType> & localIndexes, const Distribution& dist )
{
    HArray<IndexType> globalIndexes;

    auto wGlobal = hostWriteOnlyAccess( globalIndexes, localIndexes.size() );
    auto rLocal = hostReadAccess( localIndexes );

    std::transform( rLocal.begin(), rLocal.end(), wGlobal.begin(),
                    [&dist] ( IndexType localIndex )
    {
        return dist.local2Global( localIndex );
    } );

    return globalIndexes;
}

RedistributePlan::RedistributePlan( DistributionPtr targetDistribution, DistributionPtr sourceDistribution ) : 

    mSourceDistribution( sourceDistribution ), 
    mTargetDistribution( targetDistribution )

{
    SCAI_ASSERT_ERROR( sourceDistribution, "source distribution is not allowed to be null" )
    SCAI_ASSERT_ERROR( targetDistribution, "target distribution is not allowed to be null" )
    SCAI_ASSERT_EQ_ERROR( sourceDistribution->getCommunicator(), targetDistribution->getCommunicator(),
                          "source and target distributions must have the same communicator" );

    HArray<PartitionId> targetOwners;

    const auto commSize = sourceDistribution->getCommunicator().getSize();

    if ( targetDistribution->hasAnyAddressing() || commSize <= 2 )
    {
        // Computing owners is cheap, so do so directly
        HArray<PartitionId> sourceGlobalIndexes;
        sourceDistribution->getOwnedIndexes( sourceGlobalIndexes );
        targetDistribution->computeOwners( targetOwners, sourceGlobalIndexes );
    }
    else
    {
        // Building the necessary data structures for a RedistributePlan usually relies
        // on determining where to send the data. For some distributions, computing owners is cheap,
        // whereas for others (such as general distributions), this is an expensive process.
        // In the case that computing owners directly is expensive, we can recover them in
        // an asymptotically speaking far cheaper way by going through an intermediate
        // distribution for which we *can* compute the owners cheaply (e.g. block, cyclic, ...).
        //
        // The approach below is attributed to Moritz von Looz-Corswarem, who
        // pointed out the optimization opportunity to us, and provided source code and experimental
        // results to show its efficacy. The code below is loosely based on his original code.

        const auto globalSize = sourceDistribution->getGlobalSize();
        const auto comm = sourceDistribution->getCommunicatorPtr();
        const auto rank = comm->getRank();
        const auto intermediateDist = std::make_shared<BlockDistribution>( globalSize, comm );

        SCAI_LOG_INFO( logger, "build RedistributePlan via intermediate dist = " << *intermediateDist )

        RedistributePlan targetToIntermediate( intermediateDist, targetDistribution );

        // Note: source to intermediate first, then reverse
        RedistributePlan intermediateToSource( intermediateDist, sourceDistribution );
        intermediateToSource.reverse();

        // In order to find out what the owners in target of the source global indexes are,
        // we work our way backwards from target. We know that all local elements in target
        // have owner equal to the rank, and since redistribution does not change the
        // associated global index of the elements, we can simply redistribute the owners
        // (which start out as all identical to rank) through the intermediate distribution
        // and finally to the source in order to recover the desired new owner for
        // each local source index.
        HArray<PartitionId> ownersInTarget( targetDistribution->getLocalSize(), rank );
        HArray<PartitionId> ownersInIntermediate;
        targetToIntermediate.redistribute( ownersInIntermediate, ownersInTarget );

        // Reuse storage in order to possibly avoid allocation (depending on relative sizes)
        auto ownersInSource = std::move( ownersInTarget );
        intermediateToSource.redistribute( ownersInSource, ownersInIntermediate );
        targetOwners = std::move ( ownersInSource );
    }


    const auto targetGlobalIndexes = initializeFromNewOwners( targetOwners, *sourceDistribution );

    SCAI_ASSERT_DEBUG(
        HArrayUtils::all ( targetGlobalIndexes, common::CompareOp::EQ, ownedGlobalIndexesForDist( *targetDistribution ) ),
        "Internal error: mismatch between expected global indexes and target distribution" );
}



RedistributePlan::RedistributePlan( const scai::hmemo::HArray< PartitionId >& newOwnersOfLocalElements, DistributionPtr sourceDistribution )
    :   mSourceDistribution( sourceDistribution )
{
    SCAI_ASSERT_ERROR( sourceDistribution, "source distribution is not allowed to be null" );
    SCAI_ASSERT_EQ_ERROR( newOwnersOfLocalElements.size(), sourceDistribution->getLocalSize(),
                          "size of new owners must be equal to local size of distribution" );

    const auto targetGlobalIndexes = initializeFromNewOwners( newOwnersOfLocalElements, *sourceDistribution );

    mTargetDistribution = generalDistributionUnchecked( sourceDistribution->getGlobalSize(),
                                                        targetGlobalIndexes,
                                                        sourceDistribution->getCommunicatorPtr() );
}

// Note: returns global target indexes
HArray<IndexType> RedistributePlan::initializeFromNewOwners( const hmemo::HArray<PartitionId> & newOwnersOfLocalElements, const Distribution& sourceDist )
{
    const auto sourceNumLocal = sourceDist.getLocalSize();

    SCAI_ASSERT_EQ_ERROR( sourceNumLocal, newOwnersOfLocalElements.size(),
                          "Array of owners must have size equal to number of local values in source distribution." );

    HArray<IndexType> providedSourceIndexes;
    partitionSourceIndexes( mKeepSourceIndexes, providedSourceIndexes, newOwnersOfLocalElements, sourceDist );

    const auto globalKeepIndexes = local2global ( mKeepSourceIndexes, sourceDist );
    const auto globalProvidedIndexes = local2global( providedSourceIndexes, sourceDist );
    const auto newOwnersOfProvided = gather( newOwnersOfLocalElements, providedSourceIndexes );

    // We only put the exchange indexes into the Halo, as this might work considerably better when
    // most elements are kept (present both in source and target dist). This means that the Halo is working
    // with the index set given by our exchange indexes rather than local indexes of the source distribution

    HArray<IndexType> sizes;
    HArray<IndexType> perm;

    const Communicator& comm = sourceDist.getCommunicator();
    HArrayUtils::bucketSortSizes( sizes, perm, newOwnersOfProvided, comm.getSize() );

    mExchangeSendPlan = CommunicationPlan( hostReadAccess( sizes ) );
    mExchangeReceivePlan = comm.transpose( mExchangeSendPlan );

    SCAI_LOG_INFO( logger, "RedistributePlan with new owners, send plan = " << mExchangeSendPlan )

    HArray<IndexType> sortedGlobalProvidedIndexes;
    HArrayUtils::gather( sortedGlobalProvidedIndexes, globalProvidedIndexes, perm, common::BinaryOp::COPY );

    HArray<IndexType> requiredIndexes;
    comm.exchangeByPlan( requiredIndexes, mExchangeReceivePlan, sortedGlobalProvidedIndexes, mExchangeSendPlan );

    HArray<IndexType> sortedRequiredIndexes;
    HArray<IndexType> sortPermutation;
    HArrayUtils::sort( &sortPermutation, &sortedRequiredIndexes, requiredIndexes, true );

    HArray<IndexType> targetGlobalIndexes;
    HArray<IndexType> mapFromExchangeToTarget;
    HArrayUtils::mergeAndMap( targetGlobalIndexes, mapFromExchangeToTarget, mKeepTargetIndexes, sortedRequiredIndexes, globalKeepIndexes );

    // Repurpose the storage of sortedRequiredIndexes (same size and type as inversePerm) to further additional memory allocation
    auto inversePerm = std::move( sortedRequiredIndexes );
    HArrayUtils::inversePerm( inversePerm, sortPermutation );
    mExchangeSourceIndexes = gather( providedSourceIndexes, perm );
    mExchangeTargetIndexes = gather( mapFromExchangeToTarget, inversePerm );

    return targetGlobalIndexes;
}

/* -------------------------------------------------------------------------- */

void RedistributePlan::writeAt( std::ostream& stream ) const
{
    stream << "RedistributePlan( ";
    stream << *mSourceDistribution << "->" << *mTargetDistribution;
    stream << ", " << getSourceLocalSize() << "->" << getTargetLocalSize();
    stream << ", local:" << getNumLocalValues();
    stream << ", source halo :" << getExchangeSourceSize();
    stream << ", target halo :" << getExchangeTargetSize();
    stream << ")";
}

/* -------------------------------------------------------------------------- */

DistributionPtr RedistributePlan::getSourceDistributionPtr() const
{
    return mSourceDistribution;
}

/* -------------------------------------------------------------------------- */

DistributionPtr RedistributePlan::getTargetDistributionPtr() const
{
    return mTargetDistribution;
}

void RedistributePlan::reverse()
{
    std::swap( mExchangeReceivePlan, mExchangeSendPlan );
    std::swap( mSourceDistribution, mTargetDistribution );
    std::swap( mKeepSourceIndexes, mKeepTargetIndexes );
    std::swap( mExchangeSourceIndexes, mExchangeTargetIndexes );
}

} /* end namespace dmemo */

} /* end namespace scai */
