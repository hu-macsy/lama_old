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
#include <scai/dmemo/GlobalExchangePlan.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/hmemo/HostReadAccess.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostWriteOnlyAccess.hpp>

#include <memory>
#include <algorithm>

namespace scai
{

using namespace hmemo;

using utilskernel::HArrayUtils;
using dmemo::GeneralDistribution;

using std::unique_ptr;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( RedistributePlan::logger, "RedistributePlan" )

/* -------------------------------------------------------------------------- */

RedistributePlan::RedistributePlan( DistributionPtr targetDistribution, DistributionPtr sourceDistribution ) : 

    mSourceDistribution( sourceDistribution ), 
    mTargetDistribution( targetDistribution )

{
    SCAI_ASSERT_ERROR( sourceDistribution, "source distribution is not allowed to be null" )
    SCAI_ASSERT_ERROR( targetDistribution, "target distribution is not allowed to be null" )
    SCAI_ASSERT_EQ_ERROR( sourceDistribution->getCommunicator(), targetDistribution->getCommunicator(),
                          "source and target distributions must have the same communicator" );

    // Each processor computes the new owners of owned indexes from source distribution

    auto targetOwners = mTargetDistribution->owner( mSourceDistribution->ownedGlobalIndexes() );

    const auto targetGlobalIndexes = initializeFromNewOwners( targetOwners, *sourceDistribution );

    SCAI_ASSERT_ERROR(
        HArrayUtils::all ( targetGlobalIndexes, common::CompareOp::EQ, targetDistribution->ownedGlobalIndexes() ),
        "Internal error: mismatch between expected global indexes and target distribution" );
}

RedistributePlan::RedistributePlan( 
    const HArray< PartitionId >& newOwnersOfLocalElements, 
    DistributionPtr sourceDistribution ) :

    mSourceDistribution( sourceDistribution )

{
    SCAI_ASSERT_ERROR( sourceDistribution, "source distribution is not allowed to be null" );
    SCAI_ASSERT_EQ_ERROR( newOwnersOfLocalElements.size(), sourceDistribution->getLocalSize(),
                          "size of new owners must be equal to local size of distribution" );

    const auto targetGlobalIndexes = initializeFromNewOwners( newOwnersOfLocalElements, *sourceDistribution );

    mTargetDistribution = generalDistributionUnchecked( sourceDistribution->getGlobalSize(),
                                                        std::move( targetGlobalIndexes ),
                                                        sourceDistribution->getCommunicatorPtr() );
}

/* -------------------------------------------------------------------------- */

/** Help routine to extract 'local' indexes of a communication plan */

static void splitSelf( HArray<IndexType>& local, 
                       CommunicationPlan& plan, 
                       HArray<IndexType>& permutation,
                       const IndexType rank )
{
    IndexType size;
    IndexType offset;

    plan.removeEntry( size, offset, rank );

    const IndexType n = permutation.size();

    {
        auto wLocal = hostWriteOnlyAccess( local, size );
        auto wPermutation = hostWriteAccess( permutation );
 
        for ( IndexType i = 0; i < size; ++i )
        {
            wLocal[i] = wPermutation[i + offset];
        }

        for ( IndexType i = offset; i + size < n; ++i )
        {
            wPermutation[i] = wPermutation[i + size];
        }

        wPermutation.resize( n - size );
    }

    SCAI_ASSERT_EQ_DEBUG( permutation.size(), plan.totalQuantity(), "serious mismatch" );
}

// Note: returns global target indexes

HArray<IndexType> RedistributePlan::initializeFromNewOwners( 
    const hmemo::HArray<PartitionId>& newOwners, 
    const Distribution& sourceDist )
{
    HArray<PartitionId> sourceGlobalIndexes = sourceDist.ownedGlobalIndexes();

    SCAI_ASSERT_EQ_ERROR( sourceGlobalIndexes.size(), newOwners.size(),
                          "Array of owners must have size equal to number of local values in source distribution." );

    GlobalExchangePlan plan( newOwners, sourceDist.getCommunicatorPtr() );

    HArray<IndexType> inGlobalIndexes;
    plan.exchange( inGlobalIndexes, sourceGlobalIndexes );

    // targetGlobalIndexes = sort( inGlobalIndexes ), keep permutation

    HArray<IndexType> sortPerm;
    HArray<IndexType> targetGlobalIndexes;
    HArrayUtils::sort( &sortPerm, &targetGlobalIndexes, inGlobalIndexes, true );

    // the inverse permutation tells us exactly where each of the incoming elements go
    HArrayUtils::inversePerm( mExchangeTargetIndexes, sortPerm );

    // move member variables of the exchange plan to here
    plan.splitUp( mExchangeSourceIndexes, mExchangeSendPlan, mExchangeReceivePlan );

    const PartitionId rank = sourceDist.getCommunicator().getRank();

    // extract self communication, use own arrays for it

    splitSelf( mKeepSourceIndexes, mExchangeSendPlan, mExchangeSourceIndexes, rank );
    splitSelf( mKeepTargetIndexes, mExchangeReceivePlan, mExchangeTargetIndexes, rank );

    SCAI_ASSERT_EQ_ERROR( mKeepSourceIndexes.size(), mKeepTargetIndexes.size(), "serious" )

    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mKeepSourceIndexes, sourceDist.getLocalSize() ), "serious" )
    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mExchangeSourceIndexes, sourceDist.getLocalSize() ), "serious" )

    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mKeepTargetIndexes, targetGlobalIndexes.size() ), "serious" )
    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mExchangeTargetIndexes, targetGlobalIndexes.size() ), "serious" )

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

/* -------------------------------------------------------------------------- */

void RedistributePlan::reverse()
{
    std::swap( mExchangeReceivePlan, mExchangeSendPlan );
    std::swap( mSourceDistribution, mTargetDistribution );
    std::swap( mKeepSourceIndexes, mKeepTargetIndexes );
    std::swap( mExchangeSourceIndexes, mExchangeTargetIndexes );
}

} /* end namespace dmemo */

} /* end namespace scai */
