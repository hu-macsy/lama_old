/**
 * @file GlobalAddressingPlan.hpp
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
 * @brief Some template functions for typical global communication patterns
 * @author Thomas Brandes
 * @date 10.12.2018
 */
#pragma once

// internal scai libraris
#include <scai/dmemo/GlobalExchangePlan.hpp>

//#include <cmath>

namespace scai
{

namespace dmemo
{

/**
 *  A global addressing plan contains the full communication pattern for exchanging
 *  'global' data as specified by 'global' indexes. The plan might be used either
 *  for gathering non-local values or for scattering data into non-local parts.
 * 
 *  \code
 *      // A = B [ globalIndexes ]; B[globalIndexes] = A;
 *      GlobalAddressingPlan plan( globalIndexes, distB );
 *      plan.gather( A, B, comm );
 *      plan.scatter( B, A, common::BinaryOp::COPY, comm );
 *  \endcode
 */
class GlobalAddressingPlan : public GlobalExchangePlan
{
public:

    /** Construct a global addressing plan
     *
     *  @param[in] dist           is the distribution that specifies the mapping of the global indexes to local indexes
     *  @param[in] globalIndexes  is a 'distributed' array of global indexes
     *  @param[in] unique         might be set true if there is no double value among 'all' global indexes
     *
     *  Be careful that dist is not the distribution of the array globalIndexes.
     *  Even if globalindexes is also a distributed array its distribution does not matter here because
     *  only the local parts are needed.
     *
     *  Note:
     *  
     *   - global indexes out of range (dist.getGlobalSize()) will be ignored
     *   - if size() < globalIndexes.size() we know that values were out of range
     *
     *   \code
     *       auto plan = globalAddressingPlan( myGlobalIndexes, dist );
     *       const Communicator& comm = dist.getCommunicator();
     *       SCAI_ASSERT_EQ_ERROR( comm.sum( myGlobalIndexes.size(), comm.sum( plan.size() ),
     *                             "illegal indexes" )
     *   \endcode
     */
    static GlobalAddressingPlan globalAddressingPlan( 
        const Distribution& dist, 
        const hmemo::HArray<IndexType>& globalIndexes, 
        const bool unique = false );
    
    /** Constructor of a plan by its member variables 
     *
     *  @param[in] plan is a global exchange plan
     *  @param[in] localIndexes give the indexes for the local part of a distributed array
     *  @param[in] unique if true there are no double values in the array localIndexes
     *
     *  The flag unique is only used for optimization in the scatter operation. It might be set true
     *  if is guaranteed that among 'all' global indexes no entry appears twice.
     */
    GlobalAddressingPlan( GlobalExchangePlan plan, hmemo::HArray<IndexType> localIndexes, bool unique = false );

    /**
     * Default copy constructor is enabled.
     */
    GlobalAddressingPlan( const GlobalAddressingPlan& other ) = default;

    /**
     * Default move constructor is enabled.
     */
    GlobalAddressingPlan( GlobalAddressingPlan&& other ) = default;

    /**
     * Default assignment copy is enabled.
     */
    GlobalAddressingPlan& operator=( const GlobalAddressingPlan& other ) = default;

    /**
     * Default move assignment is enabled.
     */
    GlobalAddressingPlan& operator=( GlobalAddressingPlan&& other ) = default;

    /**
     *  Apply the global addressing plan for a gathering of remote data
     *
     *  @param[in,out] localArray is my local part of a 'distributed' array, size must be same to globalIndexes
     *  @param[in] remoteArray is my local part of the array that is remotely accessed
     *  @param[in] op specifies how to combine the gathered values with the existing values.
     *
     *  Note: also for op == BinaryOp::COPY the localArray should have been allocated with 
     *        sufficient size as the plan might have been built with illegal global indexes.
     */
    template<typename ValueType>
    void gather( 
        hmemo::HArray<ValueType>& localArray, 
        const hmemo::HArray<ValueType>& remoteArray, 
        const common::BinaryOp op = common::BinaryOp::COPY ) const;

    template<typename ValueType>
    void scatter( 
        hmemo::HArray<ValueType>& remoteArray, 
        const hmemo::HArray<ValueType>& localArray, 
        const common::BinaryOp op = common::BinaryOp::COPY ) const;

    /** 
     *  Optimized version for scatter with sourceArray = rank (same value), op = COPY
     */
    void scatterOwner( hmemo::HArray<PartitionId>& targetArray ) const;

    /** 
     *  Split up for member variables, so the data might be used for other purposes
     */
    void splitUp( 
        hmemo::HArray<IndexType>& sendIndexes, 
        CommunicationPlan& sendPlan, 
        CommunicationPlan& recvPlan,
        hmemo::HArray<IndexType>& recvIndexes );

private:

    hmemo::HArray<IndexType> mLocalIndexes;  // localized global addresses from other processors

    bool mUnique;    // if true there are no double entries in mLocalIndexes
};

/** 
 *   Provide GlobalAddressingPlan::globalAddressingPlan as free function
 */
inline GlobalAddressingPlan globalAddressingPlan( 
    const Distribution& dist, 
    const hmemo::HArray<IndexType>& globalIndexes, 
    const bool unique = false )
{
    return GlobalAddressingPlan::globalAddressingPlan( dist, globalIndexes, unique );
}

/* ------------------------------------------------------------------------- */
/*  Inline/template methods for GlobalAddressingPlan                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void GlobalAddressingPlan::gather( 
    hmemo::HArray<ValueType>& localArray, 
    const hmemo::HArray<ValueType>& remoteArray, 
    const common::BinaryOp op ) const
{
    // gather the values as required by global indexes used to build this plan

    auto sendValues = utilskernel::HArrayUtils::gatherF( remoteArray, mLocalIndexes );

    // send back via communication plan

    if ( op != common::BinaryOp::COPY )
    {
        SCAI_ASSERT_GE_ERROR( localArray.size(), sendSize(), "local array as target array too small" )
    }

    exchangeBack( localArray, sendValues, op );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void GlobalAddressingPlan::scatter( 
    hmemo::HArray<ValueType>& remoteArray, 
    const hmemo::HArray<ValueType>& localArray, 
    const common::BinaryOp op ) const
{
    hmemo::HArray<ValueType> recvData;
    exchange( recvData, localArray );
    utilskernel::HArrayUtils::scatter( remoteArray, mLocalIndexes, mUnique, recvData, op );
}

}

}
