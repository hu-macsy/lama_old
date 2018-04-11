/**
 * @file MetisPartitioning.cpp
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
 * @brief MetisPartitioning.cpp
 * @author ThomasBrandes
 * @date 17.08.2017
 */

// hpp
#include <scai/partitioning/MetisPartitioning.hpp>
#include <scai/partitioning/CSRGraph.hpp>
#include <scai/partitioning/CSRGraph2.hpp>

// local library

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/SingleDistribution.hpp>

// internal scai libraries

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/tracing.hpp>

#include <metis.h>

namespace scai
{

using namespace lama;
using namespace dmemo;
using namespace hmemo;

namespace partitioning
{

SCAI_LOG_DEF_LOGGER( MetisPartitioning::logger, "Partitioning.MetisPartitioning" )

#define MASTER 0

MetisPartitioning::MetisPartitioning()
{
}

void MetisPartitioning::writeAt( std::ostream& stream ) const
{
    stream << "MetisPartitioning";
}

MetisPartitioning::~MetisPartitioning()
{
    SCAI_LOG_INFO( logger, "~MetisPartitioning" )
}

/* ---------------------------------------------------------------------------------*/

void MetisPartitioning::squarePartitioning(
    HArray<PartitionId>& newLocalOwners,
    const lama::_Matrix& matrix,
    const HArray<float>& processorWeights ) const
{
    // The processor weight must be normed that they sum up to 1.0

    HArray<float> normedProcessorWeights( processorWeights );
    Partitioning::normWeights( normedProcessorWeights ); 

    const Distribution& dist = matrix.getRowDistribution();
    const Communicator& comm = dist.getCommunicator();

    IndexType nLocal  = dist.getLocalSize();
    IndexType nGlobal = dist.getGlobalSize();

    DistributionPtr singleDist( new SingleDistribution( nGlobal, dist.getCommunicatorPtr(), MASTER ) );

    CSRSparseMatrix<DefaultReal> csrMatrix;
    csrMatrix.assign( matrix );
    csrMatrix.redistribute( singleDist, singleDist );

    HArray<PartitionId> newOwners;

    SCAI_LOG_INFO( logger, comm << ", serial csrMatrix = " << csrMatrix )

    if ( comm.getRank() == MASTER )
    {
        SCAI_LOG_INFO( logger, "Master process does the partitioning" )
        squarePartitioningMaster( newOwners, csrMatrix, normedProcessorWeights );
    }
    else
    {
        SCAI_LOG_INFO( logger, comm << ": wait for results of MASTER partitioning" )
    }

    HArray<IndexType> localSizes;      // sizes of all processors required

    {
        WriteOnlyAccess<IndexType> wLocalSizes( localSizes, comm.getSize() );
        Partitioning::gatherAll( wLocalSizes.get(), nLocal, comm );
    }

    // scatter back the results from newOnwers to local parts newLocalOwners

    ReadAccess<IndexType> rSizes( localSizes );
    WriteOnlyAccess<PartitionId> wNewLocalOwners( newLocalOwners, nLocal );

    if ( comm.getRank() == MASTER )
    {
        ReadAccess<PartitionId> rNewOwners( newOwners );
        comm.scatterV( wNewLocalOwners.get(), nLocal, MASTER, rNewOwners.get(), rSizes.get() );
    }
    else
    {
        PartitionId* allOwners = 0;
        comm.scatterV( wNewLocalOwners.get(), nLocal, MASTER, allOwners, rSizes.get() );
    }
}

/* ---------------------------------------------------------------------------------*/

void MetisPartitioning::squarePartitioningMaster(
    HArray<PartitionId>& newLocalOwners,
    const lama::_Matrix& matrix,
    const HArray<float>& processorWeights ) const
{
    // Note: this method is only called by one processor

    SCAI_REGION( "Partitioning.METIS.square" )

    SCAI_LOG_INFO( logger, "Square partitioning with Metis, matrix = " << matrix << ", #parts = " << processorWeights.size() )

    // Note that here we have: global size == local size of distribution

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions( options );
    // partitioning method: Multilevel recursive bisectioning
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
    options[METIS_OPTION_NUMBERING] = 0;   // C-style, should also be default
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;  // METIS_OBJTYPE_CUT (default), or METIS_OBJTYPE_VOL

    auto tpwghts = hmemo::hostReadAccess( processorWeights ).buildVector<real_t>();

    idx_t nparts     = processorWeights.size();   // desired number of processors

    idx_t nVertices = matrix.getNumRows();

    idx_t minConstraint = 0;   // output result, e.g. edge-cut

    HArray<idx_t> partition;
    WriteOnlyAccess<idx_t> wPartition( partition, nVertices ); 

    SCAI_LOG_INFO( logger, "build CSR graph from this matrix: " << matrix )

    bool isSingle = true;
    CSRGraph<idx_t> graph( matrix, isSingle );

    const HArray<idx_t>& weights = graph.weights();

    idx_t ncon       = weights.size() / nVertices;
  
    SCAI_LOG_INFO( logger, "#vertices = " << nVertices << ", #constraints = " << ncon )

    SCAI_ASSERT_EQ_ERROR( ncon * nVertices, static_cast<idx_t>( weights.size() ), 
                          "Number of constraints = " << ncon << " does not fit" )

    ReadAccess<idx_t> rIA( graph.ia() );
    ReadAccess<idx_t> rJA( graph.ja() );
    ReadAccess<idx_t> rWeights( weights );

    real_t ubvec[]   = { 1.01 };  // ncon entries, imbalance allowed

    idx_t rc = METIS_PartGraphRecursive( &nVertices, 
                                         &ncon, 
                                         const_cast<idx_t*>( rIA.get() ),
                                         const_cast<idx_t*>( rJA.get() ),
                                         const_cast<idx_t*>( rWeights.get() ),
                                         NULL,   // no size of the vertices
                                         NULL,   // no weights for the edges
                                         &nparts,
                                         &tpwghts[0],
                                         ubvec, 
                                         options, 
                                         &minConstraint, 
                                         wPartition.get() );

    SCAI_ASSERT_EQ_ERROR( rc, METIS_OK, "METIS_PartGraphRecursive did not return normally" )

    SCAI_LOG_INFO( logger, "METIS return code = " << rc << ", objval = " << minConstraint )

    wPartition.release();

    if ( newLocalOwners.getValueType() == partition.getValueType() )
    {
        SCAI_LOG_DEBUG( logger, "move partitionining results to new local owners" )
        HArray<IndexType>& partitionT = reinterpret_cast<HArray<IndexType>&>( partition );
        newLocalOwners = std::move( partitionT );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "assign partitionining results to new local owners" )
        // HArrayUtils::assign( newLocalOwners, partition ) might not be instantiated for idx_t
        auto wNewLocalOwners = hostWriteOnlyAccess( newLocalOwners, partition.size() );
        auto rPartition      = hostReadAccess( partition );
        std::copy( rPartition.begin(), rPartition.end(), wNewLocalOwners.begin() );
    }
}

/* ---------------------------------------------------------------------------------*/

void MetisPartitioning::rectangularPartitioning( 
    HArray<PartitionId>& rowMapping,
    HArray<PartitionId>& colMapping,
    const lama::_Matrix& matrix,
    const HArray<float>& processorWeights ) const
{
    SCAI_REGION( "Partitioning.METIS.rectangular" )

    idx_t np = processorWeights.size();

    SCAI_LOG_INFO( logger, "METIS: rectangular partitioning for " << np << " processors, matrix = " << matrix )

    CSRGraph2<idx_t> graph( matrix );

    SCAI_LOG_INFO( logger, "CSRGraph2 built" )

    IndexType numRows    = matrix.getNumRows();
    IndexType numColumns = matrix.getNumColumns();

    idx_t nNodes   = numRows + numColumns;
    idx_t nWeights = graph.vertexWeights().size();
    idx_t ncon     = nWeights / nNodes;

    SCAI_ASSERT_EQ_ERROR( ncon * nNodes, nWeights, "#weights is not multiple of #nodes = " << nNodes )

    SCAI_LOG_INFO( logger, "METIS uses " << ncon << " constraints on the " << nNodes << " vertices" )

    auto metisPWeights = hmemo::hostReadAccess( processorWeights ).buildVector<real_t>();

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions( options );
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;

    HArray<idx_t> mapping( nNodes, invalidIndex );

    idx_t minConstraint = 0;
    idx_t numProcessors = metisPWeights.size();

    real_t ubvec[] = { 1.01, 1.01, 1.01 };

    {
         WriteAccess<idx_t> wMapping( mapping );
         ReadAccess<idx_t> rIA( graph.ia() );
         ReadAccess<idx_t> rJA( graph.ja() );
         ReadAccess<idx_t> rVertexWeights( graph.vertexWeights() );

         // now call METIS partitioning routine

         idx_t rc = METIS_PartGraphRecursive( &nNodes, 
                                              &ncon, 
                                              const_cast<idx_t*>( rIA.get() ), 
                                              const_cast<idx_t*>( rJA.get() ), 
                                              const_cast<idx_t*>( rVertexWeights.get() ), 
                                              NULL, 
                                              NULL, 
                                              &np,
                                              NULL, // &metisPWeights[0]
                                              ubvec, 
                                              options, 
                                              &minConstraint, 
                                              wMapping.get() );

        SCAI_ASSERT_EQ_ERROR( rc, METIS_OK, "METIS_PartGraphRecursive did not return normally" )

        SCAI_LOG_INFO( logger, "METIS return code = " << rc << ", objval = " << minConstraint )
    }

    // split mapping to new ownwers for rows and new owners for columns

    {
        ReadAccess<idx_t> rMapping( mapping );

        WriteOnlyAccess<PartitionId> wRowMapping( rowMapping, numRows );
        WriteOnlyAccess<PartitionId> wColMapping( colMapping, numColumns );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            idx_t targetProcessor = rMapping[i];
            SCAI_ASSERT_VALID_INDEX_DEBUG( rMapping[i], numProcessors, "Illegal mapping for row " << i )
            wRowMapping[i] = targetProcessor;
        }

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            idx_t targetProcessor = rMapping[numRows + j];
            SCAI_ASSERT_VALID_INDEX_DEBUG( targetProcessor, numProcessors, "Illegal mapping for col " << j )
            wColMapping[j] = targetProcessor;
        }
    }
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in Partitioning factory )    *
 * ---------------------------------------------------------------------------------*/

const char MetisPartitioning::theCreateValue[] = "METIS";

std::string MetisPartitioning::createValue()
{
    return theCreateValue;
}

PartitioningPtr MetisPartitioning::create()
{
    return PartitioningPtr( new MetisPartitioning() );
}

} /* end namespace partitioning */

} /* end namespace scai */
