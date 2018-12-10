/**
 * @file GlobalCommunication.hpp
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

/**
 *  @brief Global exchange of data between processors
 *
 *  @param[out] recvValues values received from other processors
 *  @param[in]  sendValues values that will be sent to other processors
 *  @param[in]  owners     same size as sendValues, contains for each entry the new owner
 *  @param[in]  comm       is the communicator used for exchanging data                  
 *
 *  Note: alias of recvValues and sendValues is legal
 */
template<typename ValueType>
void globalExchange(
    hmemo::HArray<ValueType>& recvValues,
    const hmemo::HArray<ValueType>& sendValues,
    const hmemo::HArray<PartitionId>& owners,
    const Communicator& comm )
{
    PartitionId np = comm.getSize();

    if ( np == 1 )
    {
        recvValues = sendValues;
        return;
    }

    // Sort the send values by their new owners (bucket sort)

    hmemo::HArray<IndexType> perm;
    hmemo::HArray<IndexType> sizes;

    utilskernel::HArrayUtils::bucketSortSizes( sizes, perm, owners, np );

    hmemo::HArray<ValueType> sortedSendValues;

    utilskernel::HArrayUtils::gather( sortedSendValues, sendValues, perm, common::BinaryOp::COPY );

    // use the sizes of the np buckets to build the send plan

    dmemo::CommunicationPlan sendPlan( hmemo::hostReadAccess( sizes ) );
    auto recvPlan = comm.transpose( sendPlan );

    comm.exchangeByPlan( recvValues, recvPlan, sortedSendValues, sendPlan );
}

template<typename ValueType1, typename ValueType2>
void globalExchange(
    hmemo::HArray<ValueType1>& recvValues1,
    hmemo::HArray<ValueType2>& recvValues2,
    const hmemo::HArray<ValueType1>& sendValues1,
    const hmemo::HArray<ValueType2>& sendValues2,
    const hmemo::HArray<PartitionId>& owners,
    const Communicator& comm ) 
{
    PartitionId np = comm.getSize();

    if ( np == 1 )
    {
        recvValues1 = sendValues1;
        recvValues2 = sendValues2;
        return;
    }

    hmemo::HArray<IndexType> perm;
    hmemo::HArray<IndexType> sizes;

    utilskernel::HArrayUtils::bucketSortSizes( sizes, perm, owners, np );

    hmemo::HArray<ValueType1> sortedSendValues1;
    hmemo::HArray<ValueType2> sortedSendValues2;

    utilskernel::HArrayUtils::gather( sortedSendValues1, sendValues1, perm, common::BinaryOp::COPY );
    utilskernel::HArrayUtils::gather( sortedSendValues2, sendValues2, perm, common::BinaryOp::COPY );

    dmemo::CommunicationPlan sendPlan( hmemo::hostReadAccess( sizes ) );
    auto recvPlan = comm.transpose( sendPlan );

    comm.exchangeByPlan( recvValues1, recvPlan, sortedSendValues1, sendPlan );
    comm.exchangeByPlan( recvValues2, recvPlan, sortedSendValues2, sendPlan );
}

template<typename ValueType1, typename ValueType2, typename ValueType3>
void globalExchange(
    hmemo::HArray<ValueType1>& recvValues1,
    hmemo::HArray<ValueType2>& recvValues2,
    hmemo::HArray<ValueType3>& recvValues3,
    const hmemo::HArray<ValueType1>& sendValues1,
    const hmemo::HArray<ValueType2>& sendValues2,
    const hmemo::HArray<ValueType3>& sendValues3,
    const hmemo::HArray<PartitionId>& owners,
    const Communicator& comm )
{
    PartitionId np = comm.getSize();

    hmemo::HArray<IndexType> perm;
    hmemo::HArray<IndexType> sizes;

    utilskernel::HArrayUtils::bucketSortSizes( sizes, perm, owners, np );

    hmemo::HArray<ValueType1> sortedSendValues1;
    hmemo::HArray<ValueType2> sortedSendValues2;
    hmemo::HArray<ValueType3> sortedSendValues3;

    utilskernel::HArrayUtils::gather( sortedSendValues1, sendValues1, perm, common::BinaryOp::COPY );
    utilskernel::HArrayUtils::gather( sortedSendValues2, sendValues2, perm, common::BinaryOp::COPY );
    utilskernel::HArrayUtils::gather( sortedSendValues3, sendValues3, perm, common::BinaryOp::COPY );

    dmemo::CommunicationPlan sendPlan( hmemo::hostReadAccess( sizes ) );
    auto recvPlan = comm.transpose( sendPlan );

    comm.exchangeByPlan( recvValues1, recvPlan, sortedSendValues1, sendPlan );
    comm.exchangeByPlan( recvValues2, recvPlan, sortedSendValues2, sendPlan );
    comm.exchangeByPlan( recvValues3, recvPlan, sortedSendValues3, sendPlan );
}

}

}

