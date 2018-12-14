/**
 * @file HaloPlan.cpp
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
 * @brief HaloPlan.cpp
 * @author Thomas Brandes
 * @date 23.02.2011
 */

// hpp
#include <scai/dmemo/HaloPlan.hpp>

#include <scai/hmemo/HostReadAccess.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( HaloPlan::logger, "HaloPlan" )

/* ---------------------------------------------------------------------- */

HaloPlan::HaloPlan()
{
}

/* ---------------------------------------------------------------------- */

HaloPlan::~HaloPlan()
{
}

/* ---------------------------------------------------------------------- */

void HaloPlan::swap( HaloPlan& other )
{
    mRequiredIndexes.swap( other.mRequiredIndexes );
    mProvidesIndexes.swap( other.mProvidesIndexes );
    mRequiredPlan.swap( other.mRequiredPlan );
    mProvidesPlan.swap( other.mProvidesPlan );
    std::swap( mGlobal2Halo, other.mGlobal2Halo );
}

/* ---------------------------------------------------------------------- */

HaloPlan::HaloPlan( 
    HArray<IndexType> requiredIndexes,
    HArray<IndexType> providesIndexes,
    CommunicationPlan requiredPlan,
    CommunicationPlan providesPlan ) :

    mRequiredIndexes( std::move( requiredIndexes ) ),
    mProvidesIndexes( std::move( providesIndexes ) ),
    mRequiredPlan( std::move( requiredPlan ) ),
    mProvidesPlan( std::move( providesPlan ) )

{
    SCAI_ASSERT_EQ_ERROR( mRequiredIndexes.size(), mRequiredPlan.totalQuantity(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( mProvidesIndexes.size(), mProvidesPlan.totalQuantity(), "serious mismatch" )

    // requiredPlan must be the tranpose of the providesPlan and vice versa, not checked here
    
    // now build the map from required indexes to halo indexes

    IndexType haloIndex = 0;

    for ( IndexType globalIndex : hostReadAccess( mRequiredIndexes ) )
    {
        mGlobal2Halo[globalIndex] = haloIndex++;
    }

    SCAI_ASSERT_EQ_ERROR( mRequiredIndexes.size(), static_cast<IndexType>( mGlobal2Halo.size() ), "double global indexes" )
}

/* ---------------------------------------------------------------------- */

void HaloPlan::splitUp(
    HArray<IndexType>& requiredIndexes,
    HArray<IndexType>& providesIndexes,
    CommunicationPlan& requiredPlan,
    CommunicationPlan& providesPlan )
{
    requiredIndexes = std::move( mRequiredIndexes );
    providesIndexes = std::move( mProvidesIndexes );
    requiredPlan = std::move( mRequiredPlan );
    providesPlan = std::move( mProvidesPlan );
}

/* ---------------------------------------------------------------------- */

void HaloPlan::clear()
{
    mRequiredPlan.clear();
    mProvidesPlan.clear();
    mRequiredIndexes.clear();
    mProvidesIndexes.clear();
    mGlobal2Halo.clear();
}

/* ---------------------------------------------------------------------- */

void HaloPlan::purge()
{
    mRequiredPlan.purge();
    mProvidesPlan.purge();
    mRequiredIndexes.purge();
    mProvidesIndexes.purge();
    // free memory of map by reallocation
    std::map<IndexType, IndexType>().swap( mGlobal2Halo );
}

/* ---------------------------------------------------------------------- */

void HaloPlan::halo2GlobalV( HArray<IndexType>& globalIndexes, const HArray<IndexType>& haloIndexes ) const
{
    const IndexType N        = haloIndexes.size();
    const IndexType haloSize = mRequiredIndexes.size();

    auto halo2global = hostReadAccess( mRequiredIndexes );
    auto rHalo       = hostReadAccess( haloIndexes );
    auto wGlobal     = hostWriteOnlyAccess( globalIndexes, N );

    #pragma omp parallel for
    for ( IndexType i = 0; i < N; ++i )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( rHalo[i], haloSize, "illegal halo index at pos " << i )
        wGlobal[i] = halo2global[ rHalo[i] ];
    }
}

/* ---------------------------------------------------------------------- */

void HaloPlan::global2HaloV( HArray<IndexType>& haloIndexes, const HArray<IndexType>& globalIndexes ) const
{
    const IndexType N = globalIndexes.size();

    auto rGlobal  = hostReadAccess( globalIndexes );
    auto wHalo    = hostWriteOnlyAccess( haloIndexes, N );

    #pragma omp parallel for
    for ( IndexType i = 0; i < N; ++i )
    {
        auto elem = mGlobal2Halo.find( rGlobal[i] );
        SCAI_ASSERT_ERROR( elem != mGlobal2Halo.end(), "global Index " << rGlobal[i] << " no halo index, never required" )
        wHalo[i] = elem->second;
    }
}

/* ---------------------------------------------------------------------- */

void HaloPlan::writeAt( std::ostream& stream ) const
{
    // write info this object
    stream << "HaloPlan( size = " << getHaloSize() << ", required plan = " << mRequiredPlan << ", provides plan = "
           << mProvidesPlan << ")";
}

/* ---------------------------------------------------------------------- */

HaloPlan HaloPlan::constructByRequiredIndexes( const HArray<IndexType>& requiredIndexes, const Distribution& distribution )
{
    const Communicator& communicator = distribution.getCommunicator();
    const PartitionId NP = communicator.getSize();

    HArray<PartitionId> owners;

    distribution.computeOwners( owners, requiredIndexes );

    // sort the required indexes by the owners via bucket sort, bucket sizes will be used for communication plan

    HArray<IndexType> perm;
    HArray<IndexType> sizes;

    utilskernel::HArrayUtils::bucketSortSizes( sizes, perm, owners, NP );

    SCAI_ASSERT_EQ_ERROR( perm.size(), owners.size(), "illegal owners." )

    // bucket sorted required indexes become part of the halo to map halo indexes back to global indexs

    HArray<IndexType> sortedRequiredIndexes;

    utilskernel::HArrayUtils::gather( sortedRequiredIndexes, requiredIndexes, perm, common::BinaryOp::COPY );

    CommunicationPlan requiredPlan = CommunicationPlan( hostReadAccess( sizes ) );
    CommunicationPlan providesPlan = communicator.transpose( requiredPlan );

    // exchange required indexes (sorted in buckets for each processor)  with other processors to get provideIndexes

    HArray<IndexType> providesIndexes;

    communicator.exchangeByPlan( providesIndexes, providesPlan, requiredIndexes, requiredPlan );

    SCAI_LOG_INFO( logger, "exchanged non-local indexes" )

    // localize the provides indexes that are still global from other processors

    distribution.global2LocalV( providesIndexes, providesIndexes );

    /*
    for ( IndexType localIndex : hostReadAccess( providesIndexes ) )
    {
        SCAI_ASSERT_NE_ERROR( localIndex, invalidIndex, "illegal local index, dist = " << distribution <<
                              ", halo plan = " << requiredPlan << ", local plan = " << providesPlan 
                              << ", providesIndexes = " << providesIndexes << ", required = " << requiredIndexes )
    }
    */

    return HaloPlan( std::move( sortedRequiredIndexes ), 
                     std::move( providesIndexes ),
                     std::move( requiredPlan ),
                     std::move( providesPlan ) );
}

template<typename ValueType>
void HaloPlan::updateHalo(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& sourceArray,
    const Communicator& comm ) const
{
    HArray<ValueType> sendValues;
    utilskernel::HArrayUtils::gather( sendValues, sourceArray, mProvidesIndexes, common::BinaryOp::COPY );
    comm.exchangeByPlan( haloArray, mRequiredPlan, sendValues, mProvidesPlan );
}

template<typename ValueType>
void HaloPlan::updateHalo(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& sourceArray,
    const Communicator& comm,
    HArray<ValueType>& tmpSendValues ) const
{
    SCAI_LOG_DEBUG( logger, "gather, halo = " << haloArray << ", source = " << sourceArray
                            << ", local Indexes = " << mProvidesIndexes
                            << " max local index = " << utilskernel::HArrayUtils::max( mProvidesIndexes ) )

    utilskernel::HArrayUtils::gather( tmpSendValues, sourceArray, mProvidesIndexes, common::BinaryOp::COPY );
    comm.exchangeByPlan( haloArray, mRequiredPlan, tmpSendValues, mProvidesPlan );
}

template<typename ValueType>
void HaloPlan::updateHaloDirect(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& sendArray,
    const Communicator& comm ) const
{
    comm.exchangeByPlan( haloArray, mRequiredPlan, sendArray, mProvidesPlan );
}

static void releaseArray( std::shared_ptr<_HArray> array )
{
    array->clear();
}

template<typename ValueType>
tasking::SyncToken* HaloPlan::updateHaloAsync(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& sourceArray,
    const Communicator& comm ) const
{
    SCAI_LOG_INFO( logger, comm << ": updateHaloAsync, source = " << sourceArray << " with this plan: " << *this )
    auto sendValues = std::make_shared<HArray<ValueType>>();
    utilskernel::HArrayUtils::gather( *sendValues, sourceArray, mProvidesIndexes, common::BinaryOp::COPY );
    SCAI_LOG_DEBUG( logger, comm << ": updateHaloAsync, send = " << *sendValues )
    auto token = comm.exchangeByPlanAsync( haloArray, mRequiredPlan, *sendValues, mProvidesPlan );
    token->pushRoutine( std::bind( releaseArray, sendValues ) );
    return token;
}

template<typename ValueType>
void HaloPlan::updateByHalo(
    HArray<ValueType>& sourceArray,
    const HArray<ValueType>& haloArray,
    common::BinaryOp op,
    const Communicator& comm ) const
{
    // receive the value at the same location where the soure array is actually valid

    HArray<ValueType> recvValues;

    comm.exchangeByPlan( recvValues, mProvidesPlan, haloArray, mRequiredPlan );
    bool unique = false;
    utilskernel::HArrayUtils::scatter( sourceArray, mProvidesIndexes, unique, recvValues, op );
}


#define SCAI_HALO_PLAN_INSTANTIATIONS( _type )                      \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void HaloPlan::updateHalo(                                      \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm ) const;                           \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void HaloPlan::updateHalo(                                      \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm,                                   \
        HArray<_type>& ) const;                                     \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void HaloPlan::updateHaloDirect(                                \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm ) const;                           \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    tasking::SyncToken* HaloPlan::updateHaloAsync(                  \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm ) const;                           \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void HaloPlan::updateByHalo(                                    \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        common::BinaryOp,                                           \
        const Communicator& comm ) const;                           \


SCAI_COMMON_LOOP( SCAI_HALO_PLAN_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_HALO_PLAN_INSTANTIATIONS


} /* end namespace dmemo */

} /* end namespace scai */
