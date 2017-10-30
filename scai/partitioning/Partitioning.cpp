/**
 * @file Partitioning.cpp
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
 * @brief Implementation of general methods for partitioning classes.
 * @author Thomas Brandes
 * @date 18.07.2017
 */

#include <scai/partitioning/Partitioning.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>

#include <scai/tracing.hpp>

namespace scai
{

using namespace dmemo;

namespace partitioning
{

/* ------  Static class variables --------------------------------------- */

SCAI_LOG_DEF_LOGGER( Partitioning::logger, "Partitioning" )

/* ------  Static methods ----------------------------------------------- */

void Partitioning::gatherWeights( std::vector<float>& weights, const float weight, const Communicator& comm )
{
    const PartitionId MASTER = 0;

    PartitionId numPartitions = comm.getSize();

    weights.resize( numPartitions );

    SCAI_LOG_INFO( logger, "#weights = " << weights.size() )

    comm.gather( &weights[0], 1, MASTER, &weight );
    comm.bcast( &weights[0], numPartitions, MASTER );
}

/* ---------------------------------------------------------------------- */

void Partitioning::normWeights( std::vector<float>& weights )
{
    const size_t numPartitions = weights.size();

    // now each partition norms it

    float sum = 0;

    for ( size_t i = 0; i < numPartitions; ++i )
    {
        sum += weights[i];
    }

    float sumNorm = 0.0f;

    for ( size_t i = 0; i < numPartitions - 1; ++i )
    {
        weights[i] /= sum;
        sumNorm += weights[i];
    }

    weights[numPartitions - 1] = 1.0f - sumNorm;
}

/* ---------------------------------------------------------------------- */

void Partitioning::rectangularRedistribute( lama::_Matrix& matrix, const float weight ) const
{
    SCAI_REGION( "partitioning.rect" )

    CommunicatorPtr comm = matrix.getRowDistribution().getCommunicatorPtr();

    if ( comm->getSize() < 2 )
    {
        // no repartitioning for a single processor
        return;
    }

    hmemo::HArray<PartitionId> rowOwners;
    hmemo::HArray<PartitionId> colOwners;

    std::vector<float> weightVector;

    gatherWeights( weightVector, weight, *comm );
    normWeights( weightVector );

    hmemo::HArrayRef<float> weights( weightVector );

    // call the 'virtual' partitioning method, i.e. the one of the derived class

    {
        SCAI_REGION( "partitioning.rect_get" )
        rectangularPartitioning( rowOwners, colOwners, matrix, weights );
    }

    {
        SCAI_REGION( "partitioning.rect_redist" )
        DistributionPtr rowDist( new GeneralDistribution( rowOwners, comm ) );
        DistributionPtr colDist( new GeneralDistribution( colOwners, comm ) );

        matrix.redistribute( rowDist, colDist );
    }
}

/* ------  Constructors  ------------------------------------------------ */

Partitioning::Partitioning()
{
}

Partitioning::~Partitioning()
{
}

void Partitioning::writeAt( std::ostream& stream ) const
{
    stream << "Partitioning( any )";
}

} /* end namespace partitioning */

} /* end namespace scai */
