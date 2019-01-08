/**
 * @file HaloExchangePlan.cpp
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
 * @brief HaloExchangePlan.cpp
 * @author Thomas Brandes
 * @date 17.12.2018
 */

// hpp
#include <scai/dmemo/HaloExchangePlan.hpp>

#include <scai/hmemo/HostReadAccess.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/TransferUtils.hpp>

#include <scai/tracing.hpp>

#include <set>

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( HaloExchangePlan::logger, "HaloExchangePlan" )

/* ---------------------------------------------------------------------- */

HaloExchangePlan::HaloExchangePlan()
{
}

/* ---------------------------------------------------------------------- */

HaloExchangePlan::~HaloExchangePlan()
{
}

/* ---------------------------------------------------------------------- */

void HaloExchangePlan::swap( HaloExchangePlan& other )
{
    mHalo2GlobalIndexes.swap( other.mHalo2GlobalIndexes );
    mLocalIndexes.swap( other.mLocalIndexes );
    mHaloCommPlan.swap( other.mHaloCommPlan );
    mLocalCommPlan.swap( other.mLocalCommPlan );
    std::swap( mGlobal2Halo, other.mGlobal2Halo );
}

/* ---------------------------------------------------------------------- */

HaloExchangePlan::HaloExchangePlan( 
    HArray<IndexType> halo2GlobalIndexes,
    HArray<IndexType> localIndexes,
    CommunicationPlan haloCommPlan,
    CommunicationPlan localCommPlan,
    std::map<IndexType, IndexType> global2Halo ) :

    mHalo2GlobalIndexes( std::move( halo2GlobalIndexes ) ),
    mLocalIndexes( std::move( localIndexes ) ),
    mHaloCommPlan( std::move( haloCommPlan ) ),
    mLocalCommPlan( std::move( localCommPlan ) ),
    mGlobal2Halo( std::move( global2Halo ) )

{
    SCAI_ASSERT_EQ_ERROR( mHalo2GlobalIndexes.size(), mHaloCommPlan.totalQuantity(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( mLocalIndexes.size(), mLocalCommPlan.totalQuantity(), "serious mismatch" )

    // haloCommPlan must be the tranpose of the localCommPlan and vice versa, not checked here
    
    const IndexType mapSize = mGlobal2Halo.size();

    if ( mapSize != mHalo2GlobalIndexes.size() )
    {
        SCAI_LOG_WARN( logger, "double global indexes: " << mapSize << " of " 
                                << mHalo2GlobalIndexes.size() << " entries are only unique" )
    }
}

/* ---------------------------------------------------------------------- */

HaloExchangePlan::HaloExchangePlan( 
    HArray<IndexType> halo2GlobalIndexes,
    HArray<IndexType> localIndexes,
    CommunicationPlan haloCommPlan,
    CommunicationPlan localCommPlan ) :

    mHalo2GlobalIndexes( std::move( halo2GlobalIndexes ) ),
    mLocalIndexes( std::move( localIndexes ) ),
    mHaloCommPlan( std::move( haloCommPlan ) ),
    mLocalCommPlan( std::move( localCommPlan ) )

{
    SCAI_ASSERT_EQ_ERROR( mHalo2GlobalIndexes.size(), mHaloCommPlan.totalQuantity(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( mLocalIndexes.size(), mLocalCommPlan.totalQuantity(), "serious mismatch" )
    
    updateMap( mGlobal2Halo, mHalo2GlobalIndexes );
}

/* ---------------------------------------------------------------------- */

void HaloExchangePlan::splitUp(
    HArray<IndexType>& requiredIndexes,
    HArray<IndexType>& localIndexes,
    CommunicationPlan& haloCommPlan,
    CommunicationPlan& localCommPlan )
{
    requiredIndexes = std::move( mHalo2GlobalIndexes );
    localIndexes = std::move( mLocalIndexes );
    haloCommPlan = std::move( mHaloCommPlan );
    localCommPlan = std::move( mLocalCommPlan );
}

/* ---------------------------------------------------------------------- */

void HaloExchangePlan::clear()
{
    mHaloCommPlan.clear();
    mLocalCommPlan.clear();
    mHalo2GlobalIndexes.clear();
    mLocalIndexes.clear();
    mGlobal2Halo.clear();
}

/* ---------------------------------------------------------------------- */

void HaloExchangePlan::purge()
{
    mHaloCommPlan.purge();
    mLocalCommPlan.purge();
    mHalo2GlobalIndexes.purge();
    mLocalIndexes.purge();
    // free memory of map by reallocation
    std::map<IndexType, IndexType>().swap( mGlobal2Halo );
}

/* ---------------------------------------------------------------------- */

void HaloExchangePlan::halo2GlobalV( HArray<IndexType>& globalIndexes, const HArray<IndexType>& haloIndexes ) const
{
    const IndexType N        = haloIndexes.size();
    const IndexType haloSize = mHalo2GlobalIndexes.size();

    auto halo2global = hostReadAccess( mHalo2GlobalIndexes );
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

void HaloExchangePlan::global2HaloV( HArray<IndexType>& haloIndexes, const HArray<IndexType>& globalIndexes ) const
{
    SCAI_REGION( "HaloExchangePlan.global2HaloV" )

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

void HaloExchangePlan::writeAt( std::ostream& stream ) const
{
    // write info this object
    stream << "HaloExchangePlan( size = " << getHaloSize() 
           << ", halo comm = " << mHaloCommPlan 
           << ", local comm = " << mLocalCommPlan << ")";
}

/* ---------------------------------------------------------------------- */

void HaloExchangePlan::buildMap( 
    std::map<IndexType, IndexType>& globalMap, 
    HArray<PartitionId>& owners,
    const HArray<IndexType>& globalIndexes,
    const PartitionId NP )
{
    IndexType localIndex = 0;

    IndexType countDoubles = 0;

    auto rGlobalIndexes = hostReadAccess( globalIndexes );

    for ( auto& owner : hostWriteAccess( owners ) )
    {
        const IndexType globalIndex = rGlobalIndexes[localIndex];

        SCAI_ASSERT_VALID_INDEX_ERROR( owner, NP, "global index " << globalIndex << " has illegal onwer" )

        if ( globalMap.find( globalIndex ) == globalMap.end() )
        {
            globalMap[globalIndex] = localIndex;   // new entry
        }
        else
        {
            owner = invalidPartition;         // double entry, elim by setting invalid owner
            countDoubles++;
        } 

        localIndex++;
    }

    IndexType mapSize = globalMap.size();

    SCAI_ASSERT_EQ_ERROR( mapSize + countDoubles, globalIndexes.size(), "serious mismatch" )

    SCAI_LOG_INFO( logger, "build map, has " << mapSize << " entries" 
                           << ", " << countDoubles << " were double" )
}

/* ---------------------------------------------------------------------- */

void HaloExchangePlan::updateMap( 
    std::map<IndexType,IndexType>& globalMap, 
    const HArray<IndexType>& globalIndexes )
{
    IndexType localIndex = 0;

    for ( auto globalIndex : hostReadAccess( globalIndexes ) )
    {
        globalMap[globalIndex] = localIndex++;
    }
}

/* ---------------------------------------------------------------------- */

HaloExchangePlan HaloExchangePlan::haloExchangePlan( 
    const Distribution& distribution,
    const HArray<IndexType>& globalIndexes, 
    const bool elimDouble )
{
    SCAI_REGION( "HaloExchangePlan.construct" )

    const Communicator& communicator = distribution.getCommunicator();
    const PartitionId NP = communicator.getSize();

    HArray<PartitionId> owners;

    distribution.computeOwners( owners, globalIndexes );

    std::map<IndexType,IndexType> globalMap;  // used for mapping global indexes to halo indexes

    if ( elimDouble )
    {
        // eliminate double indexes by setting their owner to invalid, uses std::map
        buildMap( globalMap, owners, globalIndexes, NP ); 
    }

    // sort the required indexes by the owners via bucket sort, bucket sizes will be used for communication plan

    HArray<IndexType> perm;
    HArray<IndexType> sizes;

    utilskernel::HArrayUtils::bucketSortSizes( sizes, perm, owners, NP );

    // bucket sorted required indexes become part of the halo to map halo indexes back to global indexs

    HArray<IndexType> sortedRequiredIndexes;

    utilskernel::HArrayUtils::gather( sortedRequiredIndexes, globalIndexes, perm, common::BinaryOp::COPY );

    CommunicationPlan haloCommPlan = CommunicationPlan( hostReadAccess( sizes ) );
    CommunicationPlan localCommPlan = communicator.transpose( haloCommPlan );

    // exchange required indexes (sorted in buckets for each processor)  with other processors to get provideIndexes

    HArray<IndexType> localIndexes;

    communicator.exchangeByPlan( localIndexes, localCommPlan, sortedRequiredIndexes, haloCommPlan );

    SCAI_LOG_INFO( logger, "exchanged non-local indexes" )

    // localize the provides indexes that are still global from other processors

    distribution.global2LocalV( localIndexes, localIndexes );

    updateMap( globalMap, sortedRequiredIndexes );   // map must be updated now with sorted global indexes

    return HaloExchangePlan( std::move( sortedRequiredIndexes ), 
                             std::move( localIndexes ),
                             std::move( haloCommPlan ),
                             std::move( localCommPlan ),
                             std::move( globalMap ) );
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
void HaloExchangePlan::updateHalo(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& localArray,
    const Communicator& comm ) const
{
    HArray<ValueType> sendValues;
    utilskernel::HArrayUtils::gather( sendValues, localArray, mLocalIndexes, common::BinaryOp::COPY );
    comm.exchangeByPlan( haloArray, mHaloCommPlan, sendValues, mLocalCommPlan );
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
void HaloExchangePlan::updateHalo(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& localArray,
    const Communicator& comm,
    HArray<ValueType>& tmpSendValues ) const
{
    SCAI_LOG_DEBUG( logger, "gather, halo = " << haloArray << ", source = " << localArray
                            << ", local Indexes = " << mLocalIndexes
                            << " max local index = " << utilskernel::HArrayUtils::max( mLocalIndexes ) )

    utilskernel::HArrayUtils::gather( tmpSendValues, localArray, mLocalIndexes, common::BinaryOp::COPY );
    comm.exchangeByPlan( haloArray, mHaloCommPlan, tmpSendValues, mLocalCommPlan );
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
HArray<ValueType> HaloExchangePlan::updateHaloF(
    const HArray<ValueType>& localArray,
    const Communicator& comm ) const
{
    auto sendValues = utilskernel::HArrayUtils::gatherF( localArray, mLocalIndexes );
    return comm.exchangeByPlanF( mHaloCommPlan, sendValues, mLocalCommPlan );
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
void HaloExchangePlan::updateHaloDirect(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& sendArray,
    const Communicator& comm ) const
{
    comm.exchangeByPlan( haloArray, mHaloCommPlan, sendArray, mLocalCommPlan );
}

/* ---------------------------------------------------------------------- */

static void releaseArray( std::shared_ptr<_HArray> array )
{
    array->clear();
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* HaloExchangePlan::updateHaloAsync(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& localArray,
    const Communicator& comm ) const
{
    SCAI_LOG_INFO( logger, comm << ": updateHaloAsync, source = " << localArray << " with this plan: " << *this )
    auto sendValues = std::make_shared<HArray<ValueType>>();
    utilskernel::HArrayUtils::gather( *sendValues, localArray, mLocalIndexes, common::BinaryOp::COPY );
    SCAI_LOG_DEBUG( logger, comm << ": updateHaloAsync, send = " << *sendValues )
    auto token = comm.exchangeByPlanAsync( haloArray, mHaloCommPlan, *sendValues, mLocalCommPlan );
    token->pushRoutine( std::bind( releaseArray, sendValues ) );
    return token;
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* HaloExchangePlan::updateHaloDirectAsync(
    HArray<ValueType>& haloArray,
    const HArray<ValueType>& sendArray,
    const Communicator& comm ) const
{
    return comm.exchangeByPlanAsync( haloArray, mHaloCommPlan, sendArray, mLocalCommPlan );
}

/* ---------------------------------------------------------------------- */

template<typename ValueType>
void HaloExchangePlan::updateByHalo(
    HArray<ValueType>& localArray,
    const HArray<ValueType>& haloArray,
    common::BinaryOp op,
    const Communicator& comm ) const
{
    HArray<ValueType> recvValues;   // receiving halo values from other processors

    comm.exchangeByPlan( recvValues, mLocalCommPlan, haloArray, mHaloCommPlan );
    bool unique = false;
    utilskernel::HArrayUtils::scatter( localArray, mLocalIndexes, unique, recvValues, op );
}

/* ---------------------------------------------------------------------- */

#define SCAI_HALO_PLAN_INSTANTIATIONS( _type )                      \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void HaloExchangePlan::updateHalo(                              \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm ) const;                           \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void HaloExchangePlan::updateHalo(                              \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm,                                   \
        HArray<_type>& ) const;                                     \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void HaloExchangePlan::updateHaloDirect(                        \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm ) const;                           \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    tasking::SyncToken* HaloExchangePlan::updateHaloAsync(          \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm ) const;                           \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    tasking::SyncToken* HaloExchangePlan::updateHaloDirectAsync(    \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        const Communicator& comm ) const;                           \
                                                                    \
    template COMMON_DLL_IMPORTEXPORT                                \
    void HaloExchangePlan::updateByHalo(                            \
        HArray<_type>&,                                             \
        const HArray<_type>&,                                       \
        common::BinaryOp,                                           \
        const Communicator& comm ) const;                           \


SCAI_COMMON_LOOP( SCAI_HALO_PLAN_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_HALO_PLAN_INSTANTIATIONS


} /* end namespace dmemo */

} /* end namespace scai */
