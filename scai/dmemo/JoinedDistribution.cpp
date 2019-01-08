/**
 * @file JoinedDistribution.cpp
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
 * @brief Distribution class that stands for the concatenation of two distributions
 * @author Thomas Brandes
 * @date 27.07.2017
 */

#include <scai/dmemo/JoinedDistribution.hpp>

namespace scai

{

namespace dmemo

{

SCAI_LOG_DEF_LOGGER( JoinedDistribution::logger, "Distribution.JoinedDistribution" )

JoinedDistribution::JoinedDistribution ( DistributionPtr d1, DistributionPtr d2 ) :

    Distribution( d1->getGlobalSize() + d2->getGlobalSize(), d1->getCommunicatorPtr() ),
    mD1( d1 ),
    mD2( d2 )

{
    SCAI_ASSERT_EQ_ERROR( d1->getCommunicator(), d2->getCommunicator(), "different communicators" )
}

JoinedDistribution::~JoinedDistribution()
{
}

const char* JoinedDistribution::getKind() const
{
    return "JOINED";
}

IndexType JoinedDistribution::getBlockDistributionSize() const
{
    const PartitionId numPartitions = getCommunicator().getSize();

    if ( numPartitions == 1 )
    {
        return mGlobalSize;
    }

    if ( mD1->getGlobalSize() == 0 )
    {
        return mD2->getBlockDistributionSize();
    }

    if ( mD2->getGlobalSize() == 0 )
    {
        return mD1->getBlockDistributionSize();
    }

    // usually we never have it for a joined distribution

    return invalidIndex;
}

void JoinedDistribution::writeAt( std::ostream& stream ) const
{
    stream << "JoinedDistribution( " << *mD1 << ", " << *mD2 << " )";
}

}

}
