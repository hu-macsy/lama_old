/**
 * @file GlobalAddressingPlan.cpp
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

// internal scai libraris
#include <scai/dmemo/GlobalAddressingPlan.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

namespace scai
{

using namespace hmemo;
using utilskernel::HArrayUtils;

namespace dmemo
{

GlobalAddressingPlan::GlobalAddressingPlan( GlobalExchangePlan plan, HArray<IndexType> localIndexes, const bool unique ) :

    GlobalExchangePlan( std::move( plan ) ),
    mLocalIndexes( std::move( localIndexes ) ),
    mUnique( unique )
{
    SCAI_ASSERT_EQ_ERROR( mLocalIndexes.size(), recvSize(), "serious mismatch" )
}

GlobalAddressingPlan GlobalAddressingPlan::globalAddressingPlan( 
    const Distribution& dist,
    const hmemo::HArray<IndexType>& globalIndexes, 
    const bool unique )
{
    auto exchangePlan = globalExchangePlan( dist.owner( globalIndexes ), dist.getCommunicatorPtr() );

    HArray<IndexType> localIndexes;

    exchangePlan.exchange( localIndexes, globalIndexes );

    dist.global2LocalV( localIndexes, localIndexes );

    // ToDo: verify in debug mode that if unique is true there are no double entries in localIndexes

    return GlobalAddressingPlan( std::move( exchangePlan ), std::move( localIndexes ), unique );
}

void GlobalAddressingPlan::scatterOwner( HArray<PartitionId>& targetArray )
{
    HArray<PartitionId> sources;  // will be same size as local indexes

    // here we save globalExchange operation of sourceArray = rank

    getSource( sources );   // sources[i] is the processor that required localIndexes[i]

    bool unique = true;

    HArrayUtils::scatter( targetArray, mLocalIndexes, unique, sources, common::BinaryOp::COPY );
}

void GlobalAddressingPlan::splitUp(
    HArray<IndexType>& sendIndexes,
    CommunicationPlan& sendPlan,
    CommunicationPlan& recvPlan,
    HArray<IndexType>& recvIndexes )
{
    // split up for member variables of base class

    GlobalExchangePlan::splitUp( sendIndexes, sendPlan, recvPlan );

    // and move the new member variable

    recvIndexes = std::move( mLocalIndexes );
}

}

}
