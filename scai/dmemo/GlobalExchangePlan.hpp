/**
 * @file GlobalExchangePlan.hpp
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
#include <scai/dmemo/Communicator.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

//#include <cmath>

namespace scai
{

namespace dmemo
{

/** Plan for global communication where for each entry the target processor is known. 
 *
 *  A global exchange plan between a set of processors is built by an array of processor ids 
 *  that specifies where each entry will go. When the plan is applied for an exchange of an
 *  array the values will be sent and received corresponding to this plan.
 */
class GlobalExchangePlan
{
public:

    GlobalExchangePlan( const hmemo::HArray<PartitionId>& target, const Communicator& comm );

    /** This method exchanges now data for which the plan was computed. 
     *
     *  @param[in] sendArray is the data to send, where sendArray[i] goes to target[i] 
     *  @param[in] comm      is the corresponding communicator
     *  @param[out] recvArray is the data that was received from the other processors.
     *
     *  Note: getSource might be used to find from which processor the received data came. 
     */
    template<typename ValueType>
    void exchange( hmemo::HArray<ValueType>& recvArray, const hmemo::HArray<ValueType>& sendArray, const Communicator& comm )
    {
        if ( sendArray.size() != mSendPerm.size() )
        {
            SCAI_LOG_WARN( logger, "exchange send array ( size = " << sendArray.size() 
                                   << " ) mismatches plan size = " << mSendPerm.size() )
        }

        hmemo::HArray<ValueType> sortedSendArray;
        utilskernel::HArrayUtils::gather( sortedSendArray, sendArray, mSendPerm, common::BinaryOp::COPY );
        comm.exchangeByPlan( recvArray, mRecvPlan, sortedSendArray, mSendPlan );
    } 

    /**
     *  If the plan was used to query values from other processors this routine is helpful to
     *  answer this query by using the same communication plans in the other direction.
     *  The received values match exactly the original values.
     */
    template<typename ValueType>
    void exchangeBack( hmemo::HArray<ValueType>& recvArray, const hmemo::HArray<ValueType>& sendArray, const Communicator& comm )
    {
        SCAI_ASSERT_EQ_ERROR( sendArray.size(), mRecvPlan.totalQuantity(), "serious mismatch" )

        hmemo::HArray<ValueType> recvArrayByProcs;

        comm.exchangeByPlan( recvArrayByProcs, mSendPlan, sendArray, mRecvPlan );

        utilskernel::HArrayUtils::scatter( recvArray, mSendPerm, true, recvArrayByProcs, common::BinaryOp::COPY );
    } 

    /** Get an array that contains for received values the processor id where it comes from
     *
     *  Usually this array is only needed in very special situations.
     */
    void getSource( hmemo::HArray<PartitionId>& source );

    /** 
     *  Split up for member variables, so the data might be used for other purposes
     */
    void splitUp( hmemo::HArray<IndexType>& perm, CommunicationPlan& sendPlan, CommunicationPlan& recvPlan );

private:

    hmemo::HArray<IndexType> mSendPerm;   // permutation to sort the entries according to the targets

    CommunicationPlan mSendPlan;   // plan for sending the data to exchange
    CommunicationPlan mRecvPlan;   // matching plan for receiving the data to exchange

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/**
 *  @brief Global exchange of data between processors
 *
 *  @param[out] recvValues values received from other processors
 *  @param[in]  sendValues values that will be sent to other processors
 *  @param[in]  target     same size as sendValues, contains for each entry the target processor
 *  @param[in]  comm       is the communicator used for exchanging data                  
 *
 *  Note: alias of recvValues and sendValues is allowed
 */
template<typename ValueType>
void globalExchange(
    hmemo::HArray<ValueType>& recvValues,
    const hmemo::HArray<ValueType>& sendValues,
    const hmemo::HArray<PartitionId>& target,
    const Communicator& comm )
{
    if ( comm.getSize() < 2 )
    {
        recvValues = sendValues;
        return;
    }

    GlobalExchangePlan plan( target, comm );

    plan.exchange( recvValues, sendValues, comm );
}

/**
 *  @brief Global exchange of arrays between processors, but this time with two arrays
 */
template<typename ValueType1, typename ValueType2>
void globalExchange(
    hmemo::HArray<ValueType1>& recvValues1,
    hmemo::HArray<ValueType2>& recvValues2,
    const hmemo::HArray<ValueType1>& sendValues1,
    const hmemo::HArray<ValueType2>& sendValues2,
    const hmemo::HArray<PartitionId>& target,
    const Communicator& comm ) 
{
    if ( comm.getSize() < 2 )
    {
        recvValues1 = sendValues1;
        recvValues2 = sendValues2;
        return;
    }

    GlobalExchangePlan plan( target, comm );

    plan.exchange( recvValues1, sendValues1, comm );
    plan.exchange( recvValues2, sendValues2, comm );
}

/**
 *  @brief Global exchange of arrays between processors, but this time with three arrays
 */
template<typename ValueType1, typename ValueType2, typename ValueType3>
void globalExchange(
    hmemo::HArray<ValueType1>& recvValues1,
    hmemo::HArray<ValueType2>& recvValues2,
    hmemo::HArray<ValueType3>& recvValues3,
    const hmemo::HArray<ValueType1>& sendValues1,
    const hmemo::HArray<ValueType2>& sendValues2,
    const hmemo::HArray<ValueType3>& sendValues3,
    const hmemo::HArray<PartitionId>& target,
    const Communicator& comm )
{
    if ( comm.getSize() < 2 )
    {
        recvValues1 = sendValues1;
        recvValues2 = sendValues2;
        recvValues3 = sendValues3;
        return;
    }

    GlobalExchangePlan plan( target, comm );

    plan.exchange( recvValues1, sendValues1, comm );
    plan.exchange( recvValues2, sendValues2, comm );
    plan.exchange( recvValues3, sendValues3, comm );
}

}

}

