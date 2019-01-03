/**
 * @file GlobalExchangePlan.cpp
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

#include <scai/dmemo/GlobalExchangePlan.hpp>

namespace scai
{

using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GlobalExchangePlan::logger, "GlobalExchangePlan" )

/* --------------------------------------------------------------------------- */

GlobalExchangePlan::GlobalExchangePlan(

    HArray<IndexType> sendPermutation,
    CommunicationPlan sendPlan,
    CommunicationPlan recvPlan,
    CommunicatorPtr   comm ) :

    mComm( comm ),
    mSendPerm( std::move( sendPermutation ) ),
    mSendPlan( std::move( sendPlan ) ),
    mRecvPlan( std::move( recvPlan ) )
{
    SCAI_ASSERT_EQ_ERROR( mSendPerm.size(), mSendPlan.totalQuantity(), "serious mismatch" )

    // no verification that receive plan is transpose of send plan
}

/* --------------------------------------------------------------------------- */

GlobalExchangePlan GlobalExchangePlan::globalExchangePlan( const HArray<PartitionId>& target, CommunicatorPtr comm )
{
    HArray<IndexType> sendSizes;
    HArray<IndexType> sendPerm;

    utilskernel::HArrayUtils::bucketSortSizes( sendSizes, sendPerm, target, comm->getSize() );

    auto sendPlan = CommunicationPlan( hostReadAccess( sendSizes ) );
    auto recvPlan = comm->transpose( sendPlan );

    // instead of validIndexes( target, comm.getSize() )

    if ( sendPlan.totalQuantity() < target.size() )
    {
        SCAI_LOG_WARN( logger, *comm << ": some values in target array are out-of-range" )
    } 

    return GlobalExchangePlan( std::move( sendPerm ), std::move( sendPlan ), std::move( recvPlan ), comm );
}

/* --------------------------------------------------------------------------- */

void GlobalExchangePlan::getSource( HArray<PartitionId>& source )
{
    const IndexType N = mRecvPlan.totalQuantity();  // size of received array

    auto wSource = hostWriteOnlyAccess( source, N );

    for ( PartitionId k = 0; k < mRecvPlan.size(); k++ )
    {
        const CommunicationPlan::Entry& entry = mRecvPlan[k];

        for ( IndexType i = 0; i < entry.quantity; ++i )
        {
            wSource[entry.offset + i] = entry.partitionId;
        }
    }
};

/* --------------------------------------------------------------------------- */

void GlobalExchangePlan::splitUp( HArray<IndexType>& perm, CommunicationPlan& sendPlan, CommunicationPlan& recvPlan )
{
    perm = std::move( mSendPerm );
    sendPlan = std::move( mSendPlan );
    recvPlan = std::move( mRecvPlan );
}

}

}

