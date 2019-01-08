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

/** Help routine to extract 'local' indexes of a communication plan */

void RedistributePlan::splitSelf( 
    HArray<IndexType>& local, 
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

/* -------------------------------------------------------------------------- */

RedistributePlan::RedistributePlan( 
    DistributionPtr targetDistribution,   
    hmemo::HArray<IndexType> unpackTargetPerm, 
    CommunicationPlan recvTargetPlan, 
    DistributionPtr sourceDistribution,
    hmemo::HArray<IndexType> packSourcePerm,  
    CommunicationPlan sendSourcePlan ) :

    mSourceDistribution( sourceDistribution ),
    mTargetDistribution( targetDistribution ),
  
    mExchangeSourceIndexes( std::move( packSourcePerm ) ),
    mExchangeTargetIndexes( std::move( unpackTargetPerm ) ),

    mExchangeSendPlan( std::move( sendSourcePlan ) ),
    mExchangeReceivePlan( std::move( recvTargetPlan ) )
{
    // make some checks right at the beginning

    SCAI_ASSERT_EQ_ERROR( mSourceDistribution->getGlobalSize(), mTargetDistribution->getGlobalSize(), 
                          "redistribute only possible with same sizes" )

    SCAI_ASSERT_EQ_ERROR( mSourceDistribution->getCommunicator(), mTargetDistribution->getCommunicator(), 
                          "redistribute only for distributions with same communicator" )

    SCAI_ASSERT_EQ_ERROR( mExchangeTargetIndexes.size(), mTargetDistribution->getLocalSize(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( mExchangeReceivePlan.totalQuantity(), mExchangeTargetIndexes.size(), "serious mismatch" )

    SCAI_ASSERT_EQ_ERROR( mExchangeSourceIndexes.size(), mSourceDistribution->getLocalSize(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( mExchangeSendPlan.totalQuantity(), mExchangeSourceIndexes.size(), "serious mismatch" )

    // just split up the local part, i.e. the data that is kept by this processor

    const PartitionId rank = mSourceDistribution->getCommunicator().getRank();

    // extract self communication, use own arrays for it

    splitSelf( mKeepSourceIndexes, mExchangeSendPlan, mExchangeSourceIndexes, rank );
    splitSelf( mKeepTargetIndexes, mExchangeReceivePlan, mExchangeTargetIndexes, rank );

    SCAI_ASSERT_EQ_DEBUG( mExchangeReceivePlan.totalQuantity(), mExchangeTargetIndexes.size(), "serious mismatch" )
    SCAI_ASSERT_EQ_DEBUG( mExchangeSendPlan.totalQuantity(), mExchangeSourceIndexes.size(), "serious mismatch" )

    SCAI_ASSERT_EQ_DEBUG( mKeepSourceIndexes.size(), mKeepTargetIndexes.size(), "serious mismatch" )

    SCAI_ASSERT_EQ_DEBUG( mExchangeSourceIndexes.size() + mKeepSourceIndexes.size(), 
                          mSourceDistribution->getLocalSize(), "serious mismatch" )

    SCAI_ASSERT_EQ_DEBUG( mExchangeTargetIndexes.size() + mKeepTargetIndexes.size(), 
                          mTargetDistribution->getLocalSize(), "serious mismatch" )
}

/* -------------------------------------------------------------------------- */

RedistributePlan redistributePlanByNewDistribution( 

    DistributionPtr targetDistribution, 
    DistributionPtr sourceDistribution ) 

{
    SCAI_ASSERT_ERROR( sourceDistribution, "source distribution is not allowed to be null" )
    SCAI_ASSERT_ERROR( targetDistribution, "target distribution is not allowed to be null" )

    // get the new owners of owned indexes from source distribution

    auto newOwners = targetDistribution->owner( sourceDistribution->ownedGlobalIndexes() );

    auto plan = redistributePlanByNewOwners( newOwners, sourceDistribution );

    // the new target distribution in the plan is a GeneralDistribution that must be the same
    // but we replace it with the original one that might be simpler 

    plan.resetTargetDistribution( targetDistribution );

    return plan;
}

/* -------------------------------------------------------------------------- */

void RedistributePlan::resetTargetDistribution( DistributionPtr targetDistribution )
{
    SCAI_ASSERT_ERROR(
        HArrayUtils::all ( targetDistribution->ownedGlobalIndexes(), 
                           common::CompareOp::EQ, 
                           mTargetDistribution->ownedGlobalIndexes() ),

        "mismatch of new target distribution" )

    mTargetDistribution = targetDistribution;
}

/* -------------------------------------------------------------------------- */

RedistributePlan redistributePlanByNewOwners( 
    const HArray< PartitionId >& newOwners,
    DistributionPtr sourceDistribution ) 
{
    SCAI_ASSERT_ERROR( sourceDistribution, "source distribution is not allowed to be null" );

    SCAI_ASSERT_EQ_ERROR( newOwners.size(), sourceDistribution->getLocalSize(),
                          "new owner for each owned index required" );

    auto sourceGlobalIndexes = sourceDistribution->ownedGlobalIndexes();

    SCAI_ASSERT_EQ_ERROR( sourceGlobalIndexes.size(), newOwners.size(),
                          "Array of owners must have size equal to number of local values in source distribution." );

    auto plan = globalExchangePlan( newOwners, sourceDistribution->getCommunicatorPtr() );

    auto inGlobalIndexes = plan.exchangeF( sourceGlobalIndexes );

    // targetGlobalIndexes = sort( inGlobalIndexes ), keep permutation

    HArray<IndexType> sortPerm;
    HArray<IndexType> targetGlobalIndexes;
    HArrayUtils::sort( &sortPerm, &targetGlobalIndexes, inGlobalIndexes, true );

    // the sorted incoming global indexes give the new 'GeneralDistribution', no further checks required

    auto targetDistribution = generalDistributionUnchecked( 
        sourceDistribution->getGlobalSize(),
        std::move( targetGlobalIndexes ),
        sourceDistribution->getCommunicatorPtr() );

    // the inverse permutation tells us exactly where each of the incoming elements go, use it for plan

    HArray<IndexType> unpackTargetPerm;
    HArrayUtils::inversePerm( unpackTargetPerm, sortPerm );

    // move member variables of the exchange plan to here

    CommunicationPlan sendPlan;
    CommunicationPlan recvPlan;
    HArray<IndexType> sourcePackPerm;

    plan.splitUp( sourcePackPerm, sendPlan, recvPlan );

    return RedistributePlan( targetDistribution, std::move( unpackTargetPerm ), std::move( recvPlan ),
                             sourceDistribution, std::move( sourcePackPerm ), std::move( sendPlan ) );
}

/* -------------------------------------------------------------------------- */

void RedistributePlan::writeAt( std::ostream& stream ) const
{
    stream << "RedistributePlan( ";
    stream << *mSourceDistribution << "->" << *mTargetDistribution;
    stream << ", keeps:" << mKeepSourceIndexes.size();
    stream << ", source exchg :" << mExchangeSourceIndexes.size();
    stream << ", target exchg :" << mExchangeTargetIndexes.size();
    stream << ")";
}

/* -------------------------------------------------------------------------- */

DistributionPtr RedistributePlan::getTargetDistributionPtr() const
{
    return mTargetDistribution;
}

DistributionPtr RedistributePlan::getSourceDistributionPtr() const
{
    return mSourceDistribution;
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
