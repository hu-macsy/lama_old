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

GlobalAddressingPlan::GlobalAddressingPlan( 

    const HArray<IndexType>& globalIndexes, 
    const Distribution& dist ) :

    GlobalExchangePlan( dist.owner( globalIndexes ), dist.getCommunicatorPtr() )

{
    // use the built exchange plan to send the required global indexes to the processors
    exchange( mLocalIndexes, globalIndexes );

    // the global indexes required from the other processors will now be localized
    dist.global2LocalV( mLocalIndexes, mLocalIndexes );
}

void GlobalAddressingPlan::scatterOwner( HArray<PartitionId>& targetArray )
{
    HArray<PartitionId> sources;  // will be same size as local indexes

    // here we save globalExchange operation of sourceArray = rank

    getSource( sources );   // sources[i] is the processor that required localIndexes[i]

    bool unique = true;

    HArrayUtils::scatter( targetArray, mLocalIndexes, unique, sources, common::BinaryOp::COPY );
}

}

}
