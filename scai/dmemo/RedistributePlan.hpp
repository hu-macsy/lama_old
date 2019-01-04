/**
 * @file RedistributePlan.hpp
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
 * @brief Implementation of a class that handles redistribution of matrices and vectors
 * @author Thomas Brandes
 * @date 08.10.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// local library
#include <scai/dmemo/Distribution.hpp>

// internal scai libraries
#include <scai/hmemo/HArray.hpp>
#include <scai/utilskernel/TransferUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/tracing.hpp>

#include <memory>

namespace scai
{

namespace dmemo
{

/** This class allows to create objects that handle redistributions of vector and
 *  matrices from one distribution into another distribution.
 *
 *  Once created, it has built internal data structures like communicaton plans
 *  that restrict the redistribution just to the transfer of the corresponding data.
 */

class COMMON_DLL_IMPORTEXPORT RedistributePlan: public common::Printable
{
public:

    /** Build an object by the source and target distribution.
     *
     *  @param[in] targetDistribution  the new distribution
     *  @param[in] sourceDistribution  the old distribution
     *
     *  The global size of both distributions must be the same.
     */
    RedistributePlan( DistributionPtr targetDistribution, DistributionPtr sourceDistribution );

    /**
     * Build a RedistributePlan by a source distribution and a list of owners.
     *
     * By supplying the new owners instead of the target distribution, the expensive operation
     * of computing owners for elements can be avoided. The target distribution is instead built
     * as part of the construction of the RedistributePlan.
     *
     * This is only useful in some cases - for example when integrating with external redistribution software,
     * which might output exactly such a list of new owners. In these cases, the necessary communication required
     * to create a RedistributePlan may sometimes be significantly reduced.
     *
     * @param[in] newOwners array with new owner for currently owned global indexes
     * @param[in] sourceDistribution  The source distribution from which to redistribute elements.
     *
     *  *  newOwners.size() == sourceDistribution.getLocalSize()
     *  *  newOwners[i] is new owner of sourceDistribution.local2Global( i ), 
     */
    RedistributePlan( const hmemo::HArray<PartitionId>& newOwners, DistributionPtr sourceDistribution );

    /** Getter needed for distributions */

    DistributionPtr getTargetDistributionPtr() const;
    DistributionPtr getSourceDistributionPtr() const;

    IndexType getTargetLocalSize() const
    {
        return mTargetDistribution->getLocalSize();
    }

    /** Redistribution of a distributed vector as HArrays.
     *
     *  @param[out] targetArray  vector in target distribution
     *  @param[in]  sourceArray  vector in source distribution
     *
     *  Size of sourceArray must be at least sourceDistribution->getLocalSize()
     *  Size of targetArray must be at least targetDistribution->getLocalSize()
     *
     *  This routine is a simple straightforward solution. More
     *  efficient implementations are possible: overlapping local copies  and
     *  communication.
     */
    template<typename ValueType>
    void redistribute( hmemo::HArray<ValueType>& targetArray, const hmemo::HArray<ValueType>& sourceArray ) const;

    /** Redistribution of a distributed vector as HArrays.
     *
     *  @param[out] targetArray  vector in target distribution
     *  @param[in]  sourceArray  vector in source distribution
     *  @param[in]  n            TODO[doxy] Complete Description.
     *
     *  Size of sourceArray must be at least sourceDistribution->getLocalSize()
     *  Size of targetArray must be at least targetDistribution->getLocalSize()
     *
     *  This code only demonstrates how the redistributor schedule works. More
     *  efficient implementations are possible: overlapping local copies  and
     *  communication.
     */
    template<typename ValueType>
    void redistributeN( hmemo::HArray<ValueType>& targetArray, const hmemo::HArray<ValueType>& sourceArray, IndexType n ) const;

    /** Redistribution of ragged arrays. */

    template<typename ValueType>
    void redistributeV(
        hmemo::HArray<ValueType>& targetArray,
        const hmemo::HArray<IndexType>& targetOffsets,
        const hmemo::HArray<ValueType>& sourceArray,
        const hmemo::HArray<IndexType>& sourceOffsets ) const;

    const hmemo::HArray<IndexType>& getKeepSourceIndexes() const
    {
        return mKeepSourceIndexes;
    }

    const hmemo::HArray<IndexType>& getKeepTargetIndexes() const
    {
        return mKeepTargetIndexes;
    }

    const hmemo::HArray<IndexType>& getExchangeSourceIndexes() const
    {
        return mExchangeSourceIndexes;
    }

    const hmemo::HArray<IndexType>& getExchangeTargetIndexes() const
    {
        return mExchangeTargetIndexes;
    }

    /** Reverse the RedistributePlan.
     *
     * Has the same effect as redistributing with RedistributePlan( source, target ), except
     * that the reverse operation can be performed significantly cheaper than constructing
     * a new RedistributePlan in this way.
     */
    void reverse();

private:

    /**
     *  @param[in] sourceDist source distribution is the current distribution 
     *  @param[in] newOwners array with new owner for currently owned indexes
     */
    hmemo::HArray<IndexType> initializeFromNewOwners( 
        const Distribution& sourceDist,
        const hmemo::HArray< PartitionId >& newOwners );

    virtual void writeAt( std::ostream& stream ) const;

    DistributionPtr mSourceDistribution;
    DistributionPtr mTargetDistribution;

    // sorted local indices in source dist which remain local
    // (i.e. the corresponding global index is also local in target)
    hmemo::HArray<IndexType> mKeepSourceIndexes;

    // sorted local indices in target dist which remain local
    // (i.e. the corresponding global index is also local in source)
    hmemo::HArray<IndexType> mKeepTargetIndexes;

    // local indices in source/target, sorted according to the send/receive plans
    hmemo::HArray<IndexType> mExchangeSourceIndexes;
    hmemo::HArray<IndexType> mExchangeTargetIndexes;

    CommunicationPlan mExchangeSendPlan;
    CommunicationPlan mExchangeReceivePlan;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void RedistributePlan::redistribute( hmemo::HArray<ValueType>& targetArray, const hmemo::HArray<ValueType>& sourceArray ) const
{
    SCAI_REGION( "RedistributePlan.redistribute" )

    targetArray.clear();                        // just invalidate all valid data
    targetArray.resize( getTargetLocalSize() ); // resize it correctly

    const Communicator& comm = mSourceDistribution->getCommunicator();

    // allocate memory for source (provides) and target (required) that is exchanged

    hmemo::HArray<ValueType> sourceHalo( mExchangeSourceIndexes.size() );
    hmemo::HArray<ValueType> targetHalo( mExchangeTargetIndexes.size() );

    utilskernel::HArrayUtils::gather( sourceHalo, sourceArray, mExchangeSourceIndexes, common::BinaryOp::COPY );
    utilskernel::TransferUtils::copy( targetArray, mKeepTargetIndexes, sourceArray, mKeepSourceIndexes );
    comm.exchangeByPlan( targetHalo, mExchangeReceivePlan, sourceHalo, mExchangeSendPlan );

    const bool uniqueIndexes = true;   // redistribution is always permutation without double indexes

    utilskernel::HArrayUtils::scatter( targetArray, mExchangeTargetIndexes, uniqueIndexes, targetHalo, common::BinaryOp::COPY );
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void RedistributePlan::redistributeN(
    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<ValueType>& sourceArray,
    IndexType n ) const
{
    using utilskernel::TransferUtils;

    SCAI_REGION( "RedistributePlan.redistributeN" )

    const Communicator& comm = mSourceDistribution->getCommunicator();

    targetArray.clear();                            // just invalidate all valid data
    targetArray.resize( getTargetLocalSize() * n ); // resize it correctly

    // allocate memory for source (provides) and target (required) halo

    hmemo::HArray<ValueType> sourceHalo( n * mExchangeSourceIndexes.size() );
    hmemo::HArray<ValueType> targetHalo( n * mExchangeTargetIndexes.size() );

    SCAI_LOG_DEBUG( logger, "gather: sourceHalo " << mExchangeSourceIndexes.size() << " * " << n << " values" )
    TransferUtils::gatherN( sourceHalo, sourceArray, mExchangeSourceIndexes, n );
    SCAI_LOG_DEBUG( logger, "copy: source -> target " << mKeepTargetIndexes.size() << " * " << n << " values" )
    TransferUtils::copyN( targetArray, mKeepTargetIndexes, sourceArray, mKeepSourceIndexes, n );

    // Communication plans are built by multiplication with n

    auto recvPlanN = mExchangeReceivePlan.constructN( n );
    auto sendPlanN = mExchangeSendPlan.constructN( n );

    comm.exchangeByPlan( targetHalo, recvPlanN, sourceHalo, sendPlanN );

    SCAI_LOG_DEBUG( logger, "scatter: targetHalo " << mExchangeTargetIndexes.size() << " * " << n << " values" )
    TransferUtils::scatterN( targetArray, mExchangeTargetIndexes, targetHalo, n );
}

/* ------------------------------------------------------------------------------- */

template<typename ValueType>
void RedistributePlan::redistributeV(
    hmemo::HArray<ValueType>& targetArray,
    const hmemo::HArray<IndexType>& targetOffsets,
    const hmemo::HArray<ValueType>& sourceArray,
    const hmemo::HArray<IndexType>& sourceOffsets ) const
{
    using utilskernel::TransferUtils;

    SCAI_REGION( "RedistributePlan.redistributeV" )

    const Communicator& comm = mSourceDistribution->getCommunicator();

    {
        auto t = hmemo::hostWriteOnlyAccess<ValueType>( targetArray, targetOffsets[getTargetLocalSize()] );
    }

    hmemo::HArray<IndexType> sourceSendSizes;
    hmemo::HArray<IndexType> receiveTargetSizes;

    hmemo::ContextPtr hostCtx = hmemo::Context::getHostPtr();

    sparsekernel::CSRUtils::gatherSizes( sourceSendSizes, sourceOffsets, mExchangeSourceIndexes, hostCtx );
    sparsekernel::CSRUtils::gatherSizes( receiveTargetSizes, targetOffsets, mExchangeTargetIndexes, hostCtx );

    auto sourcePlanV = mExchangeSendPlan.constructRaggedBySizes( sourceSendSizes );
    auto targetPlanV = mExchangeReceivePlan.constructRaggedBySizes( receiveTargetSizes );

    hmemo::HArray<ValueType> sourceHalo( sourcePlanV.totalQuantity() );
    hmemo::HArray<ValueType> targetHalo( targetPlanV.totalQuantity() );

    TransferUtils::gatherV( sourceHalo, sourceArray, sourceOffsets, mExchangeSourceIndexes );

    TransferUtils::copyV( targetArray, targetOffsets, mKeepTargetIndexes,
                          sourceArray, sourceOffsets, mKeepSourceIndexes );

    comm.exchangeByPlan( targetHalo, targetPlanV, sourceHalo, sourcePlanV );

    TransferUtils::scatterV( targetArray, targetOffsets, mExchangeTargetIndexes, targetHalo );
}


} /* end namespace dmemo */

} /* end namespace scai */
