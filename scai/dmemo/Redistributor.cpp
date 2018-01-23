/**
 * @file Redistributor.cpp
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
 * @brief Implementation of methods for Redistributor class.
 * @author Thomas Brandes
 * @date 08.10.2011
 */

// hpp
#include <scai/dmemo/Redistributor.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

// local library
#include <scai/dmemo/HaloBuilder.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>

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

static void partitionKeepAndRequired( HArray<IndexType> & keepSourceIndexes,
                                      HArray<IndexType> & keepTargetIndexes,
                                      HArray<IndexType> & globalRequiredIndexes,
                                      const Distribution & sourceDist,
                                      const Distribution & targetDist )
{
    const auto sourceNumLocal = sourceDist.getLocalSize();
    const auto targetNumLocal = targetDist.getLocalSize();

    // Determine the maximum sizes that keep/required can have. We allocate
    // for this size, and then later resize the arrays to their actual values.
    const auto maxKeepSize = std::min( sourceNumLocal, targetNumLocal );
    const auto maxRequiredSize = targetNumLocal;

    auto wKeepSource = hostWriteOnlyAccess( keepSourceIndexes, maxKeepSize );
    auto wKeepTarget = hostWriteOnlyAccess( keepTargetIndexes, maxKeepSize );
    auto wRequired = hostWriteOnlyAccess( globalRequiredIndexes, maxRequiredSize );

    IndexType numKeep = 0;
    IndexType numRequired = 0;

    for ( IndexType localTargetIndex = 0; localTargetIndex < targetDist.getLocalSize(); ++localTargetIndex )
    {
        const auto globalIndex = targetDist.local2global( localTargetIndex );
        const auto localSourceIndex = sourceDist.global2local(globalIndex);

        if ( localSourceIndex != invalidIndex )
        {
            wKeepSource[numKeep] = localSourceIndex;
            wKeepTarget[numKeep++] = localTargetIndex;
        }
        else
        {
            wRequired[numRequired++] = globalIndex;
        }
    }

    wKeepSource.resize( numKeep );
    wKeepTarget.resize( numKeep );
    wRequired.resize( numRequired );
}


SCAI_LOG_DEF_LOGGER( Redistributor::logger, "Redistributor" )

Redistributor::Redistributor( DistributionPtr targetDistribution, DistributionPtr sourceDistribution )

    : mSourceDistribution( sourceDistribution ), mTargetDistribution( targetDistribution )

{
    SCAI_ASSERT_ERROR( sourceDistribution, "NULL pointer for source distribution" )
    SCAI_ASSERT_ERROR( targetDistribution, "NULL pointer for target distribution" )

    const Distribution& sourceDist = *sourceDistribution;
    const Distribution& targetDist = *targetDistribution;
    SCAI_ASSERT_EQ_ERROR( sourceDist.getGlobalSize(), targetDist.getGlobalSize(), "source and target distributions must have the same global size" );
    SCAI_ASSERT_EQ_ERROR( sourceDist.getCommunicator(), targetDist.getCommunicator(), "source and target distributions must have the same communicator" );

    HArray<IndexType> globalRequiredIndexes;
    partitionKeepAndRequired( mKeepSourceIndexes, mKeepTargetIndexes, globalRequiredIndexes, sourceDist, targetDist);
    SCAI_ASSERT_EQ_ERROR( mKeepSourceIndexes.size(), mKeepTargetIndexes.size(), "keep indexes size mismatch");
    const auto numLocalValues = mKeepSourceIndexes.size();

    HaloBuilder::build( sourceDist, globalRequiredIndexes, mHalo );
    mExchangeSourceIndexes = mHalo.getProvidesIndexes();

    ContextPtr contextPtr = Context::getHostPtr();
    WriteAccess<IndexType> exchangeTargetIndexes( mExchangeTargetIndexes, contextPtr );

    // Now add the indexes where to scatter the halo into destination
    IndexType haloSize = mHalo.getHaloSize();
    SCAI_ASSERT_ERROR( numLocalValues + haloSize == getTargetLocalSize(), "size mismatch" )
    exchangeTargetIndexes.resize( haloSize );

    const auto rGlobalRequiredIndexes = hostReadAccess( globalRequiredIndexes );

    for ( IndexType i = 0; i < haloSize; i++ )
    {
        IndexType globalIndex = rGlobalRequiredIndexes[i];
        IndexType localIndex = targetDist.global2local( globalIndex );
        IndexType haloIndex = mHalo.global2halo( globalIndex );
        SCAI_LOG_TRACE( logger,
                        "saved mapping for target dist: local = " << localIndex << ", global = " << globalIndex << ", halo = " << haloIndex )
        exchangeTargetIndexes[haloIndex] = localIndex;
    }

    // all info needed for redistribution is now available
    SCAI_LOG_INFO( logger, "constructed " << *this )
}

// Partitions local indexes of the source distribution into "keep"
// and "exchange", based on the new owner of each individual index.
static void partitionLocalIndexes(HArray<IndexType> & keepLocalIndexes,
                                  HArray<IndexType> & exchangeLocalIndexes,
                                  const HArray<PartitionId> & newOwnersOfLocalElements,
                                  const Distribution & sourceDist)
{
    SCAI_ASSERT_EQ_DEBUG(newOwnersOfLocalElements.size(), sourceDist.getLocalSize(), "sourceDist and newOwners must have same size");
    const auto rank = sourceDist.getCommunicator().getRank();
    const auto numPartitions = sourceDist.getCommunicator().getSize();
    const auto sourceNumLocal = sourceDist.getLocalSize();

    const auto rNewOwners = hostReadAccess(newOwnersOfLocalElements);
    auto wKeep = hostWriteOnlyAccess(keepLocalIndexes, sourceNumLocal);
    auto wExchange = hostWriteOnlyAccess(exchangeLocalIndexes, sourceNumLocal);

    IndexType numKeep = 0;
    IndexType numExchange = 0;

    for (IndexType localSourceIndex = 0; localSourceIndex < sourceDist.getLocalSize(); ++localSourceIndex)
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

    wKeep.resize(numKeep);
    wExchange.resize(numExchange);
}

template <typename ValueType>
static HArray<ValueType> selectIndexes( const HArray<ValueType> & source, const HArray<IndexType> & indexes )
{
    HArray<ValueType> result;
    HArrayUtils::gather( result, source, indexes, common::BinaryOp::COPY );
    return result;
}

// TODO: Make this a method of Distribution (with name local2global)? (can give default impl, but
// allow subclasses to override it for a more efficient implementation)
static void local2globalInto(HArray<IndexType> & indexes, const Distribution & dist)
{
    for (auto & x : hostWriteAccess(indexes))
    {
        x = dist.local2global(x);
    }
}

static HArray<IndexType> local2global( const HArray<IndexType> & localIndexes, const Distribution & dist )
{
    auto result = localIndexes;
    local2globalInto( result, dist );
    return result;
}

Redistributor::Redistributor( const scai::hmemo::HArray< PartitionId >& newOwnersOfLocalElements, DistributionPtr sourceDistribution )
    :   mSourceDistribution(sourceDistribution)
{
    const auto & sourceDist = *sourceDistribution;
    const auto sourceNumLocal = sourceDistribution->getLocalSize();

    SCAI_ASSERT_EQ_ERROR( sourceNumLocal, newOwnersOfLocalElements.size(),
                          "Array of owners must have size equal to number of local values in source distribution." );

    HArray<IndexType> providedSourceIndexes;
    partitionLocalIndexes( mKeepSourceIndexes, providedSourceIndexes, newOwnersOfLocalElements, *sourceDistribution );

    const auto globalKeepIndexes = local2global (mKeepSourceIndexes, sourceDist );
    const auto globalProvidedIndexes = local2global( providedSourceIndexes, sourceDist );
    const auto newOwnersOfProvided = selectIndexes( newOwnersOfLocalElements, providedSourceIndexes );

    // We only put the exchange indexes into the Halo, as this might work considerably better when
    // most elements are kept (present both in source and target dist). This means that the Halo is working
    // with the index set given by our exchange indexes rather than local indexes of the source distribution
    HaloBuilder::buildFromProvidedOwners( sourceDistribution->getCommunicator(), globalProvidedIndexes, newOwnersOfProvided, mHalo );

    HArray<IndexType> sortedRequiredIndexes;
    HArray<IndexType> sortPermutation;
    HArrayUtils::sort( &sortPermutation, &sortedRequiredIndexes, mHalo.getRequiredIndexes(), true );

    HArray<IndexType> targetGlobalIndexes;
    HArray<IndexType> mapFromExchangeToTarget;
    HArrayUtils::mergeAndMap( targetGlobalIndexes, mapFromExchangeToTarget, mKeepTargetIndexes, sortedRequiredIndexes, globalKeepIndexes );

    // Repurpose the storage of sortedRequiredIndexes (same size and type as inversePerm) to further additional memory allocation
    auto inversePerm = std::move( sortedRequiredIndexes );
    HArrayUtils::inversePerm(inversePerm, sortPermutation);
    mExchangeSourceIndexes = selectIndexes( providedSourceIndexes, mHalo.getProvidesIndexes() );
    mExchangeTargetIndexes = selectIndexes( mapFromExchangeToTarget, inversePerm );
    mTargetDistribution = DistributionPtr ( new GeneralDistribution( sourceDistribution->getGlobalSize(),
                                                                     targetGlobalIndexes,
                                                                     sourceDistribution->getCommunicatorPtr() ) );
}

/* -------------------------------------------------------------------------- */

void Redistributor::writeAt( std::ostream& stream ) const
{
    stream << "Redistributor( ";
    stream << *mSourceDistribution << "->" << *mTargetDistribution;
    stream << ", " << getSourceLocalSize() << "->" << getTargetLocalSize();
    stream << ", local:" << getNumLocalValues();
    stream << ", source halo :" << getExchangeSourceSize();
    stream << ", target halo :" << getExchangeTargetSize();
    stream << ")";
}

/* -------------------------------------------------------------------------- */

DistributionPtr Redistributor::getSourceDistributionPtr() const
{
    return mSourceDistribution;
}

/* -------------------------------------------------------------------------- */

DistributionPtr Redistributor::getTargetDistributionPtr() const
{
    return mTargetDistribution;
}

/* -------------------------------------------------------------------------- */

void Redistributor::buildVPlans( const IndexType haloSourceSizes[], const IndexType haloTargetSizes[] ) const
{
    const IndexType numProvides = mHalo.getProvidesPlan().totalQuantity();
    const IndexType numRequired = mHalo.getRequiredPlan().totalQuantity();
    // calculate number of provided and required values by summing up the corresponding quantities
    unique_ptr<IndexType[]> provideQuantities( new IndexType[numProvides] );
    unique_ptr<IndexType[]> requiredQuantities( new IndexType[numRequired] );

    // For building the new schedule we need the sizes, can be calculated by the offsets

    for ( IndexType i = 0; i < numProvides; i++ )
    {
        IndexType size = haloSourceSizes[i];
        provideQuantities[i] = size;
        SCAI_LOG_DEBUG( logger, "provides[" << i << "] = " << size )
    }

    for ( IndexType i = 0; i < numRequired; i++ )
    {
        IndexType size = haloTargetSizes[i];
        requiredQuantities[i] = size;
        SCAI_LOG_DEBUG( logger, "required[" << i << "] = " << size )
    }

    mProvidesPlan.reset( new CommunicationPlan( mHalo.getProvidesPlan(), provideQuantities.get() ) );
    mRequiredPlan.reset( new CommunicationPlan( mHalo.getRequiredPlan(), requiredQuantities.get() ) );
    SCAI_LOG_INFO( logger, "providesPlan = " << *mProvidesPlan )
    SCAI_LOG_INFO( logger, "requiredPlan = " << *mRequiredPlan )
}

/* -------------------------------------------------------------------------- */

void Redistributor::buildRowPlans(
    const HArray<IndexType>& targetSizes,
    const HArray<IndexType>& sourceSizes ) const
{
    const IndexType numProvides = mHalo.getProvidesPlan().totalQuantity();
    const IndexType numRequired = mHalo.getRequiredPlan().totalQuantity();
    // calculate number of provided and required values by summing up the corresponding quantities
    unique_ptr<IndexType[]> provideQuantities( new IndexType[numProvides] );
    unique_ptr<IndexType[]> requiredQuantities( new IndexType[numRequired] );
    ContextPtr contextPtr = Context::getHostPtr();
    // For building the new schedule we need the sizes, can be calculated by the offsets
    {
        ReadAccess<IndexType> indexes( mExchangeSourceIndexes, contextPtr );
        ReadAccess<IndexType> sizes( sourceSizes, contextPtr );

        for ( IndexType i = 0; i < numProvides; i++ )
        {
            IndexType size = sizes[indexes[i]];
            provideQuantities[i] = size;
            SCAI_LOG_DEBUG( logger, "provides[" << i << "] = " << size )
        }
    }
    {
        ReadAccess<IndexType> indexes( mExchangeTargetIndexes, contextPtr );
        ReadAccess<IndexType> sizes( targetSizes, contextPtr );

        for ( IndexType i = 0; i < numRequired; i++ )
        {
            IndexType size = sizes[indexes[i]];
            requiredQuantities[i] = size;
            SCAI_LOG_DEBUG( logger, "required[" << i << "] = " << size )
        }
    }
    mProvidesPlan.reset( new CommunicationPlan( mHalo.getProvidesPlan(), provideQuantities.get() ) );
    mRequiredPlan.reset( new CommunicationPlan( mHalo.getRequiredPlan(), requiredQuantities.get() ) );
    SCAI_LOG_INFO( logger, "providesPlan = " << *mProvidesPlan )
    SCAI_LOG_INFO( logger, "requiredPlan = " << *mRequiredPlan )
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

/*
 template COMMON_DLL_IMPORTEXPORT
 void Redistributor::redistributeN ( HArray<float>& targetArray,
 const HArray<float>& sourceArray,
 IndexType n ) const;

 template COMMON_DLL_IMPORTEXPORT
 void Redistributor::redistributeN ( HArray<double>& targetArray,
 const HArray<double>& sourceArray,
 IndexType n ) const;
 */

} /* end namespace dmemo */

} /* end namespace scai */
