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
     *  @param[in] indexes is a 'distributed' array of global indexes
     *  @param[in] dist    is the distribution that specifies the mapping of the global indexes to local indexes
     *
     *  Even if indexes is also a distributed array its distribution does not matter here because
     *  only the local parts are needed.
     *
     *  ToDo: explain what happens if there are illegal global indexes
     *  
     *   - global indexes out of range will be ignored
     *   - if size() < globalIndexes.size() we know that values were out of range
     *
     *   \code
     *       GlobalAddressingPlan plan( myIndexes, dist );
     *       const Communicator& comm = dist.getCommunicator();
     *       SCAI_ASSERT_EQ_ERROR( comm.sum( myIndexes.size(), comm.sum( plan.size() ),
     *                             "illegal indexes" )
     *   \endcode
     */
    GlobalAddressingPlan( const hmemo::HArray<IndexType>& globalIndexes, const Distribution& dist );

    GlobalAddressingPlan( const GlobalAddressingPlan& other ) = default;

    GlobalAddressingPlan( GlobalAddressingPlan&& other ) = default;

    /**
     *  Apply the global addressing plan for a gathering of remote data
     *
     *  @param[out] localArray is my local part of a 'distributed' array, size must be same to globalIndexes
     *  @param[in] remoteArray is my local part of the array that is remotely accessed
     */
    template<typename ValueType>
    void gather( 
        hmemo::HArray<ValueType>& localArray, 
        const hmemo::HArray<ValueType>& remoteArray, 
        const common::BinaryOp op = common::BinaryOp::COPY );

    template<typename ValueType>
    void scatter( 
        hmemo::HArray<ValueType>& remoteArray, 
        const hmemo::HArray<ValueType>& localArray, 
        const common::BinaryOp op = common::BinaryOp::COPY );

    /** 
     *  Optimized version for scatter with sourceArray = rank (same value), op = COPY
     */
    void scatterOwner( hmemo::HArray<PartitionId>& targetArray );

private:

    hmemo::HArray<IndexType> mLocalIndexes;  // localized global addresses from other processors
};

template<typename ValueType>
void GlobalAddressingPlan::gather( 
    hmemo::HArray<ValueType>& localArray, 
    const hmemo::HArray<ValueType>& remoteArray, 
    const common::BinaryOp op )
{
    hmemo::HArray<ValueType> sendValues;

    // gather the values as required by global indexes used to build this plan

    utilskernel::HArrayUtils::gather( sendValues, remoteArray, mLocalIndexes, common::BinaryOp::COPY );

    // send back via communication plan

    SCAI_ASSERT_GE_ERROR( localArray.size(), size(), "target array too small" )

    exchangeBack( localArray, sendValues, op );
}

template<typename ValueType>
void GlobalAddressingPlan::scatter( 
    hmemo::HArray<ValueType>& remoteArray, 
    const hmemo::HArray<ValueType>& localArray, 
    const common::BinaryOp op )
{
    hmemo::HArray<ValueType> recvData;
    exchange( recvData, localArray );
    utilskernel::HArrayUtils::scatter( remoteArray, mLocalIndexes, false, recvData, op );
}

}

}
