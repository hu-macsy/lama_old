/**
 * @file BlockPartitioning.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief BlockPartitioning.cpp
 * @author ThomasBrandes
 * @date 17.08.2017
 */

// hpp
#include <scai/partitioning/BlockPartitioning.hpp>

// local library
#include <scai/dmemo/GenBlockDistribution.hpp>

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/tracing.hpp>

namespace scai
{

using namespace lama;
using namespace dmemo;

namespace partitioning
{

SCAI_LOG_DEF_LOGGER( BlockPartitioning::logger, "Partitioning.BlockPartitioning" )

BlockPartitioning::BlockPartitioning()
{
}

BlockPartitioning::~BlockPartitioning()
{
    SCAI_LOG_INFO( logger, "~BlockPartitioning" )
}

void BlockPartitioning::writeAt( std::ostream& stream ) const
{
    stream << "BlockPartitioning";
}

DistributionPtr BlockPartitioning::partitionIt( const CommunicatorPtr comm, const Matrix& matrix, float weight ) const
{
    IndexType globalSize = matrix.getRowDistribution().getGlobalSize();

    // Block partitioning : just create a 'general' block distribution
    // Note: this does not take the connections into account 

    return DistributionPtr( new GenBlockDistribution( globalSize, weight, comm ) );
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in Partitioning factory )    *
 * ---------------------------------------------------------------------------------*/

const char BlockPartitioning::theCreateValue[] = "BLOCK";

std::string BlockPartitioning::createValue()
{
    return theCreateValue;
}

PartitioningPtr BlockPartitioning::create()
{
    return PartitioningPtr( new BlockPartitioning() );
}

} /* end namespace partitioning */

} /* end namespace scai */
