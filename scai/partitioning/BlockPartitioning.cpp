/**
 * @file BlockPartitioning.cpp
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
using namespace hmemo;

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

/* ---------------------------------------------------------------------------------*/

void BlockPartitioning::rectangularPartitioning( 
    HArray<PartitionId>& rowMapping,
    HArray<PartitionId>& colMapping,
    const _Matrix& matrix,
    const HArray<float>& processorWeights ) const
{
    IndexType npart = processorWeights.size();

    IndexType numRows = matrix.getNumRows();
    IndexType numCols = matrix.getNumColumns();

    IndexType rowBlockSize = ( numRows + npart - 1 ) / npart;
    IndexType colBlockSize = ( numCols + npart - 1 ) / npart;

    // distribute numRows according to the weights on the processors
    // distribute numColumns according to the weights on the processors
 
    const Distribution& rowDist = matrix.getRowDistribution();
    const Distribution& colDist = matrix.getColDistribution();

    IndexType numLocalRows = rowDist.getLocalSize();
    IndexType numLocalCols = colDist.getLocalSize();

    WriteOnlyAccess<PartitionId> wRowMapping( rowMapping, numLocalRows );

    for ( IndexType i = 0; i < numLocalRows; ++i )
    {
        IndexType globalI = rowDist.local2global( i );
        wRowMapping[i] = globalI / rowBlockSize;
    }

    WriteOnlyAccess<PartitionId> wColMapping( colMapping, numLocalCols );

    for ( IndexType j = 0; j < numLocalCols; ++j )
    {
        IndexType globalJ = colDist.local2global( j );
        wColMapping[j] = globalJ / colBlockSize;
    }
}

/* ---------------------------------------------------------------------------------*/

void BlockPartitioning::squarePartitioning(
    HArray<PartitionId>& newLocalOwners,
    const _Matrix& matrix,
    const HArray<float>& processorWeights ) const
{
    IndexType npart = processorWeights.size();

    SCAI_ASSERT_EQ_ERROR( matrix.getNumRows(), matrix.getNumColumns(), "square partitioning only for square matrices" )

    IndexType n = matrix.getNumRows();

    IndexType nb = ( n + npart - 1 ) / npart;  // block size

    // weights is not taken into account yet
    // ToDo: make the mapping like it is done for a general block distribution with weights

    const Distribution& dist = matrix.getRowDistribution();

    IndexType nLocal = dist.getLocalSize();

    HArray<IndexType> ownedIndexes;
    dist.getOwnedIndexes( ownedIndexes );

    WriteOnlyAccess<PartitionId> wMapping( newLocalOwners, nLocal );
    ReadAccess<IndexType> rOwners( ownedIndexes );

    for ( IndexType i = 0; i < nLocal; ++i )
    {
        wMapping[i] = rOwners[i] / nb;
    }
}

/* ---------------------------------------------------------------------------------*/

void BlockPartitioning::rectangularRedistribute( _Matrix& matrix, const float weight ) const
{
    CommunicatorPtr comm = matrix.getRowDistribution().getCommunicatorPtr();

    if ( comm->getSize() < 2 )
    {
        // no repartitioning for a single processor
        return;
    }

    IndexType numRows    = matrix.getRowDistribution().getGlobalSize();
    IndexType numColumns = matrix.getColDistribution().getGlobalSize();

    // Block partitioning : just create a 'general' block distribution
    // Note: this does not take the connections into account 

    DistributionPtr rowDist( new GenBlockDistribution( numRows, weight, comm ) );
    DistributionPtr colDist( new GenBlockDistribution( numColumns, weight, comm ) );

    matrix.redistribute( rowDist, colDist );
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
