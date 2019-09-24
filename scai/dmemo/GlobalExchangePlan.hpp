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
#include <scai/tracing.hpp>

//#include <cmath>

namespace scai
{

namespace dmemo
{

/** Plan for global communication where for each entry of the array the target processor is specified.
 *
 *  A global exchange plan between a set of processors is built by an array of processor ids 
 *  that specifies where each entry will go. When the plan is applied for an exchange of an
 *  array the values will be sent and received corresponding to this plan.
 *
 *  One central concept of using a global exchange plan is to pack the data for the processors
 *  into contiguous sections. 
 */
class GlobalExchangePlan
{
public:

    /**
     *  Construct a plan by specifying a target processor for each element.
     *
     *  @param[in] target    array of processor ids to which sending elements will go
     *  @param[in] comm      is the corresponding communicator
     *
     *  Some notes:
     *
     *  - all processors of the communicator comm must call this method as it implies also
     *    an all-to-all communication for building communication plans (matching receive).
     *  - each processor has its indiviual array of targets                                 
     *  - if there is an invalid target specified it will be ignored.
     */
    static GlobalExchangePlan globalExchangePlan( 
        const hmemo::HArray<PartitionId>& target, 
        CommunicatorPtr comm = Communicator::getCommunicatorPtr() );

    /**
     *  @brief Constructor for a global exchange plan by the corresponding member variables
     *
     *  @param[in] sendPermutation specifies how to pack the elements into contiguous parts for each processor
     *  @param[in] sendPlan specifies the sizes to send to each other processor
     *  @param[in] recvPlan is the counterpart to sendPlan (i.e. the transpose)
     *  @param[in] comm specifies the processor set for the data exchange.
     */
    GlobalExchangePlan( 
        hmemo::HArray<IndexType> sendPermutation, 
        CommunicationPlan sendPlan,
        CommunicationPlan recvPlan,
        CommunicatorPtr comm );

    /** Default copy constructor works fine as all member variables have a copy constructor. */

    GlobalExchangePlan( const GlobalExchangePlan& ) = default;

    /** Default move constructor works fine as all member variables have a move constructor. */

    GlobalExchangePlan( GlobalExchangePlan&& ) = default;

    /** Default copy assignment works fine as all member variables support it. */

    GlobalExchangePlan& operator=( const GlobalExchangePlan& ) = default;

    /** Default move assignment works fine as all member variables support it. */

    GlobalExchangePlan& operator=( GlobalExchangePlan&& ) = default;

    /**
     *  Query the number of values for that this processor will send
     */
    inline IndexType sendSize() const;

    /**
     *  Query the number of elements that specified this processor (me) as target
     */
    inline IndexType recvSize() const;

    /** This method exchanges now data for which the plan was computed. 
     *
     *  @param[in] sendArray is the data to send, where sendArray[i] goes to target[i] 
     *  @param[out] recvArray is the data that was received from the other processors.
     *
     *  Note: getSource might be used to find from which processor the received data came. 
     */
    template<typename ValueType>
    void exchange( hmemo::HArray<ValueType>& recvArray, const hmemo::HArray<ValueType>& sendArray ) const;

    /** Provide exchange as a function for convenience. */

    template<typename ValueType>
    hmemo::HArray<ValueType> exchangeF( const hmemo::HArray<ValueType>& sendArray ) const;

    /**
     *  If the plan was used to query values from other processors this routine is helpful to
     *  answer this query by using the same communication plans in the other direction.
     *  The received values match exactly the original values.
     */
    template<typename ValueType>
    void exchangeBack( 
        hmemo::HArray<ValueType>& recvArray, 
        const hmemo::HArray<ValueType>& sendArray, 
        const common::BinaryOp op = common::BinaryOp::COPY ) const;

    /** Get an array that contains for received values the processor id where it comes from
     *
     *  Usually this array is only needed in very special situations.
     */
    void getSource( hmemo::HArray<PartitionId>& source ) const;

    /** 
     *  Split up for member variables, so the data might be used for other purposes
     */
    void splitUp( hmemo::HArray<IndexType>& perm, CommunicationPlan& sendPlan, CommunicationPlan& recvPlan );

private:

    CommunicatorPtr mComm;                //!< keep the communicator for convenience

    hmemo::HArray<IndexType> mSendPerm;   //!<  permutation to sort the entries according to the targets

    CommunicationPlan mSendPlan;   //!<  plan for sending the (contiguous) data to exchange
    CommunicationPlan mRecvPlan;   //!<  matching plan for receiving the (contiguous) data to exchange

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* ------------------------------------------------------------------------- */
/* Free function to construct a global exchange plan                         */
/* ------------------------------------------------------------------------- */

/**
 *  Provide GlobalExchangePlan::globalExchangePlan as free function.
 */
inline GlobalExchangePlan globalExchangePlan( 
    const hmemo::HArray<PartitionId>& target, 
    CommunicatorPtr comm = Communicator::getCommunicatorPtr() )
{
    SCAI_REGION( "GlobalExchangePlan.build" )
    return GlobalExchangePlan::globalExchangePlan( target, comm );
}

/* ------------------------------------------------------------------------- */
/* Template methods for global exchange of data using a global exchange plan */
/* ------------------------------------------------------------------------- */

/**
 *  @brief Global exchange of data between processors
 *
 *  Each processor calls this methods with its own values that will be sent
 *  to other processors. The values sendValues[i] will be sent to processor
 *  target[i]. Each processor gets in recvValues the values received from
 *  the other processors. 
 *
 *  @param[out] recvValues values received from other processors
 *  @param[in]  sendValues values that will be sent to other processors
 *  @param[in]  target     same size as sendValues, contains for each entry the target processor
 *  @param[in]  comm       communicator used for the data exchange
 *
 *  Note: alias of recvValues and sendValues is allowed
 */
template<typename ValueType>
void globalExchange(
    hmemo::HArray<ValueType>& recvValues,
    const hmemo::HArray<ValueType>& sendValues,
    const hmemo::HArray<PartitionId>& target,
    dmemo::CommunicatorPtr comm )
{
    SCAI_REGION( "GlobalExchangePlan.globalExchange" )

    if ( comm->getSize() < 2 )
    {
        recvValues = sendValues;
        return;
    }

    auto plan = globalExchangePlan( target, comm );
    plan.exchange( recvValues, sendValues );
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
    CommunicatorPtr comm ) 
{
    SCAI_REGION( "GlobalExchangePlan.globalExchange2" )

    if ( comm->getSize() < 2 )
    {
        recvValues1 = sendValues1;
        recvValues2 = sendValues2;
        return;
    }

    auto plan = globalExchangePlan( target, comm );

    plan.exchange( recvValues1, sendValues1 );
    plan.exchange( recvValues2, sendValues2 );
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
    CommunicatorPtr comm )
{
    SCAI_REGION( "GlobalExchangePlan.globalExchange3" )
    if ( comm->getSize() < 2 )
    {
        recvValues1 = sendValues1;
        recvValues2 = sendValues2;
        recvValues3 = sendValues3;
        return;
    }

    auto plan = globalExchangePlan( target, comm );

    plan.exchange( recvValues1, sendValues1 );
    plan.exchange( recvValues2, sendValues2 );
    plan.exchange( recvValues3, sendValues3 );
}

/* ------------------------------------------------------------------------- */
/*  Inline/template methods for GlobalExchangePlan                           */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void GlobalExchangePlan::exchange( hmemo::HArray<ValueType>& recvArray, const hmemo::HArray<ValueType>& sendArray ) const
{
    SCAI_REGION( "GlobalExchangePlan.exchange" )

    if ( sendArray.size() != mSendPerm.size() )
    {
        SCAI_LOG_WARN( logger, "exchange send array ( size = " << sendArray.size() 
                               << " ) mismatches plan size = " << mSendPerm.size() )
    }

    hmemo::HArray<ValueType> sortedSendArray;
    utilskernel::HArrayUtils::gather( sortedSendArray, sendArray, mSendPerm, common::BinaryOp::COPY );
    mComm->exchangeByPlan( recvArray, mRecvPlan, sortedSendArray, mSendPlan );
} 

template<typename ValueType>
hmemo::HArray<ValueType> GlobalExchangePlan::exchangeF( const hmemo::HArray<ValueType>& sendArray ) const
{
    hmemo::HArray<ValueType> recvArray;
    exchange( recvArray, sendArray );
    return recvArray;
}

template<typename ValueType>
void GlobalExchangePlan::exchangeBack( 
    hmemo::HArray<ValueType>& recvArray, 
    const hmemo::HArray<ValueType>& sendArray, 
    const common::BinaryOp op ) const
{
    SCAI_REGION( "GlobalExchangePlan.exchangeBack" )

    SCAI_ASSERT_EQ_ERROR( sendArray.size(), mRecvPlan.totalQuantity(), "serious mismatch" )

    hmemo::HArray<ValueType> recvArrayByProcs;

    mComm->exchangeByPlan( recvArrayByProcs, mSendPlan, sendArray, mRecvPlan );

    bool unique = true;   // as mSendPerm is a permutation

    utilskernel::HArrayUtils::scatter( recvArray, mSendPerm, unique, recvArrayByProcs, op );
} 

IndexType GlobalExchangePlan::sendSize() const
{
    return mSendPlan.totalQuantity();
}

IndexType GlobalExchangePlan::recvSize() const
{
    return mRecvPlan.totalQuantity();
}

}

}

