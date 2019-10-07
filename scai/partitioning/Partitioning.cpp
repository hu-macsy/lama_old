/**
 * @file Partitioning.cpp
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
 * @brief Implementation of general methods for partitioning classes.
 * @author Thomas Brandes
 * @date 18.07.2017
 */

#include <scai/partitioning/Partitioning.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/RedistributePlan.hpp>

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
    PartitionId numPartitions = comm.getSize();

    weights.resize( numPartitions );

    gatherAll( &weights[0], weight, comm );
}

/* ---------------------------------------------------------------------- */

void Partitioning::gatherWeights( hmemo::HArray<float>& weights, const float weight, const Communicator& comm )
{
    PartitionId numPartitions = comm.getSize();

    hmemo::WriteOnlyAccess<float> writeWeights( weights, numPartitions );

    gatherAll( writeWeights.get(), weight, comm );
}

/* ---------------------------------------------------------------------- */

void Partitioning::rectangularRedistribute( lama::_Matrix& matrix, const float weight ) const
{
    SCAI_REGION( "partitioning.rect" )

    if ( matrix.getRowDistribution().isReplicated() )
    {
        // no repartitioning for a single processor

        return;
    }

    CommunicatorPtr comm = matrix.getRowDistribution().getCommunicatorPtr();

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

        const PartitionId root = 0;

        auto rowDist = generalDistributionBySingleOwners( rowOwners, root, comm );
        auto colDist = generalDistributionBySingleOwners( colOwners, root, comm );

        matrix.redistribute( rowDist, colDist );
    }
}

/* ---------------------------------------------------------------------------------*/

void Partitioning::squarePartitioning(
    hmemo::HArray<PartitionId>& newLocalOwners,
    const lama::_Matrix& matrix,
    const float weight ) const
{
    const Communicator& comm = matrix.getRowDistribution().getCommunicator();

    hmemo::HArray<float> processorWeights;

    Partitioning::gatherWeights( processorWeights, weight, comm );

    // check that all weights are greater 0 

    SCAI_LOG_DEBUG( logger, "weights gathered = " << processorWeights );

    normWeights( processorWeights );

    SCAI_LOG_DEBUG( logger, "now call squarePartitioning with all processor weights = " << processorWeights );

    squarePartitioning( newLocalOwners, matrix, processorWeights );
}

/* ---------------------------------------------------------------------- */

void Partitioning::squarePartitioningW(
    hmemo::HArray<PartitionId>& newLocalOwners,
    const lama::_Matrix& matrix,
    const hmemo::HArray<float>& vertexWeights,
    const float weight ) const
{
    const Communicator& comm = matrix.getRowDistribution().getCommunicator();

    hmemo::HArray<float> processorWeights;

    Partitioning::gatherWeights( processorWeights, weight, comm );

    // check that all weights are greater 0 

    SCAI_LOG_DEBUG( logger, "weights gathered = " << processorWeights );

    normWeights( processorWeights );

    SCAI_LOG_DEBUG( logger, "now call squarePartitioning with all processor weights = " << processorWeights );

    squarePartitioningW( newLocalOwners, matrix, vertexWeights, processorWeights );
}

/* ---------------------------------------------------------------------- */

void Partitioning::squarePartitioningW(
    hmemo::HArray<PartitionId>& newLocalOwners,
    const lama::_Matrix& matrix,
    const hmemo::HArray<float>&,
    const hmemo::HArray<float>& weight ) const
{
    SCAI_LOG_ERROR( logger, "vertex weights are ignored" )
    squarePartitioning( newLocalOwners, matrix, weight );
}

/* ---------------------------------------------------------------------- */

dmemo::DistributionPtr Partitioning::partitionIt( 
    const dmemo::CommunicatorPtr comm, 
    const lama::_Matrix& matrix, 
    float weight ) const
{
    hmemo::HArray<float> processorWeights;   // replicated array with all weights needed

    gatherWeights( processorWeights, weight, *comm );
    normWeights( processorWeights );

    hmemo::HArray<IndexType> newLocalOwners;

    squarePartitioning( newLocalOwners, matrix, processorWeights );

    // build a new 'general' distribution by the new owners

    return dmemo::generalDistributionByNewOwners( matrix.getRowDistribution(), newLocalOwners );
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

/* ------  Instantiations ------------------------------------------------ */


} /* end namespace partitioning */

} /* end namespace scai */
