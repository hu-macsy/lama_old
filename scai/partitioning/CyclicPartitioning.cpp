/**
 * @file CyclicPartitioning.cpp
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
 * @brief CyclicPartitioning.cpp
 * @author ThomasBrandes
 * @date 17.08.2017
 */

// hpp
#include <scai/partitioning/CyclicPartitioning.hpp>

// local library
#include <scai/dmemo/CyclicDistribution.hpp>

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

SCAI_LOG_DEF_LOGGER( CyclicPartitioning::logger, "Partitioning.BlockPartitioning" )

static inline IndexType getCyclicBlockSize( IndexType globalSize, IndexType np )
{
    // choose block size in such a way that each processor has roughly 10 blocks

    IndexType nb = static_cast<IndexType>( globalSize / ( 10.0 * np ) + 0.5 );

    if ( nb < 1 )
    {
        nb = 1;
    }

    return nb;
}

CyclicPartitioning::CyclicPartitioning()
{
}

CyclicPartitioning::~CyclicPartitioning()
{
    SCAI_LOG_INFO( logger, "~CyclicPartitioning" )
}

void CyclicPartitioning::writeAt( std::ostream& stream ) const
{
    stream << "CyclicPartitioning";
}

DistributionPtr CyclicPartitioning::partitionIt( const CommunicatorPtr comm, const _Matrix& matrix, float ) const
{
    IndexType globalSize = matrix.getRowDistribution().getGlobalSize();

    // Cyclic partitioning : just create a 'cyclic' block distribution

    IndexType np = comm->getSize();

    // Each processor should get roughly 10 chunks

    IndexType nb = getCyclicBlockSize( globalSize, np );

    return DistributionPtr( new CyclicDistribution( globalSize, nb, comm ) );
}

/* ---------------------------------------------------------------------------------*/

void CyclicPartitioning::rectangularPartitioning( 
    HArray<PartitionId>& rowMapping,
    HArray<PartitionId>& colMapping,
    const _Matrix& matrix,
    const HArray<float>& processorWeights ) const
{
    IndexType npart = processorWeights.size();

    IndexType numRows = matrix.getNumRows();
    IndexType numCols = matrix.getNumColumns();

    IndexType nbRows = getCyclicBlockSize( numRows, npart );
    IndexType nbCols = getCyclicBlockSize( numCols, npart );

    // weights cannot be taken into account for cyclic distributions
 
    const Distribution& rowDist = matrix.getRowDistribution();
    const Distribution& colDist = matrix.getColDistribution();

    IndexType numLocalRows = rowDist.getLocalSize();
    IndexType numLocalCols = colDist.getLocalSize();

    WriteOnlyAccess<PartitionId> wRowMapping( rowMapping, numLocalRows );

    for ( IndexType i = 0; i < numLocalRows; ++i )
    {
        IndexType globalI = rowDist.local2global( i );
        IndexType globalChunk = globalI / nbRows;
        wRowMapping[i] = globalChunk % npart;
    }

    WriteOnlyAccess<PartitionId> wColMapping( colMapping, numLocalCols );

    for ( IndexType j = 0; j < numLocalCols; ++j )
    {
        IndexType globalJ = colDist.local2global( j );
        IndexType globalChunk = globalJ / nbCols;
        wColMapping[j] = globalChunk % npart;
    }
}

/* ---------------------------------------------------------------------------------*/

void CyclicPartitioning::rectangularRedistribute( _Matrix& matrix, const float ) const
{
    CommunicatorPtr comm = matrix.getRowDistribution().getCommunicatorPtr();

    IndexType np = comm->getSize();

    if ( np < 2 )
    {
        // no repartitioning for a single processor
        return;
    }

    IndexType numRows    = matrix.getRowDistribution().getGlobalSize();
    IndexType numColumns = matrix.getColDistribution().getGlobalSize();

    // Cyclic partitioning : just create a 'cyclic' block distribution

    IndexType nbRows = getCyclicBlockSize( numRows, np );
    IndexType nbCols = getCyclicBlockSize( numColumns, np );

    DistributionPtr rowDist( new CyclicDistribution( numRows, nbRows, comm ) );
    DistributionPtr colDist( new CyclicDistribution( numColumns, nbCols, comm ) );

    matrix.redistribute( rowDist, colDist );
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in Partitioning factory )    *
 * ---------------------------------------------------------------------------------*/

const char CyclicPartitioning::theCreateValue[] = "CYCLIC";

std::string CyclicPartitioning::createValue()
{
    return theCreateValue;
}

PartitioningPtr CyclicPartitioning::create()
{
    return PartitioningPtr( new CyclicPartitioning() );
}

} /* end namespace partitioning */

} /* end namespace scai */
