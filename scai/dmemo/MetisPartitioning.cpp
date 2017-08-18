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
#include <scai/dmemo/MetisPartitioning.hpp>

// local library
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/tracing.hpp>

extern "C"
{
#include <metis.h>
}

namespace scai
{

namespace dmemo
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

DistributionPtr MetisPartitioning::partitionIt( const CommunicatorPtr comm, const Distributed& matrix, float weight ) const
{
    std::vector<float> weights;
    gatherWeights( weights, weight, *comm );
    normWeights( weights );

    return computeIt( comm, matrix, weights );
}

DistributionPtr MetisPartitioning::computeIt( const CommunicatorPtr comm, const Distributed& matrix, std::vector<float>& weights ) const
{
    // TODO check only for replicated matrix
    IndexType size = comm->getSize();
    IndexType myRank = comm->getRank();
    IndexType totalRows = matrix.getDistribution().getGlobalSize();

    if ( myRank == MASTER )
    {
        // MASTER must have all values of the matrix to build the graph
        SCAI_ASSERT_EQ_ERROR( totalRows, matrix.getDistribution().getLocalSize(),
                              *comm << ": must have all values of the distributed object" );
    }

    std::vector<IndexType> numRowsPerOwner;
    std::vector<IndexType> distIA;
    std::vector<IndexType> rows;

    if ( myRank == MASTER )
    {
        // test weights (tpwgts)
        // note: no weight can be zero for metis call
        // so reduce number of partitions and map partition index to the right process
        std::vector<real_t> tpwgts( size );
        std::vector<IndexType> mapping( size );
        IndexType count = 0;
        checkAndMapWeights( tpwgts, mapping, count, weights, size );
        // set size only for MASTER
        numRowsPerOwner.resize( size );
        distIA.resize( size + 1 );
        rows.resize( totalRows );

        if ( count > 1 )
        {
            IndexType minConstraint = 0;
            std::vector<IndexType> partition( totalRows );
            callPartitioning( partition, minConstraint, count, tpwgts, comm, matrix );
            IndexType offsetCounter[size];

            // init
            for ( IndexType i = 0; i < size; i++ )
            {
                numRowsPerOwner[i] = 0;
                offsetCounter[i] = 0;
            }

            //count rows per owner
            for ( IndexType i = 0; i < totalRows; i++ )
            {
                ++numRowsPerOwner[mapping[partition[i]]];
            }

            // build "ia" array (running sums) for rows per owner
            distIA[0] = 0;

            for ( IndexType i = 1; i < size + 1; i++ )
            {
                distIA[i] = distIA[i - 1] + numRowsPerOwner[i - 1];
            }

            // sort rows after owner
            for ( IndexType i = 0; i < totalRows; i++ )
            {
                IndexType index = mapping[partition[i]];
                rows[distIA[index] + offsetCounter[index]] = i;
                ++offsetCounter[index];
            }
        }
        else
        {
            SCAI_LOG_WARN( logger,
                           "MetisPartitioning called with 1 processor/1 weight, which is the same as NoPartitioning." )

            for ( IndexType i = 0; i < size; i++ )
            {
                numRowsPerOwner[i] = 0;
            }

            numRowsPerOwner[mapping[0]] = totalRows;

            for ( IndexType i = 0; i < totalRows; i++ )
            {
                rows[i] = i;
            }
        }
    }

    hmemo::HArray<IndexType> myGlobalIndexes;

    // scatter local rows

    int numMyRows = 0;

    comm->scatter( &numMyRows, 1, MASTER, &numRowsPerOwner[0] );

    {
        hmemo::WriteOnlyAccess<IndexType> wLocal2Global( myGlobalIndexes, numMyRows );
        comm->scatterV( wLocal2Global.get(), numMyRows, MASTER, &rows[0], &numRowsPerOwner[0] );
    }

    // verify that indexes are sorted otherwise global to local via binary search won't work

    SCAI_ASSERT_ERROR( utilskernel::HArrayUtils::isSorted( myGlobalIndexes, common::binary::LT ), 
                       *comm << ": local row indexes unsorted." )

    IndexType globalSize = matrix.getDistribution().getGlobalSize();

    return DistributionPtr( new GeneralDistribution( globalSize, myGlobalIndexes, comm ) );
}

MetisPartitioning::~MetisPartitioning()
{
    SCAI_LOG_INFO( logger, "~MetisPartitioning" )
}

template<typename WeightType>
void MetisPartitioning::callPartitioning(
    std::vector<IndexType>& partition,
    IndexType& minConstraint,
    IndexType& parts,
    std::vector<WeightType>& tpwgts,
    const CommunicatorPtr,
    const Distributed& matrix ) const
{
    SCAI_REGION( "METIScall" )
    // Note that here we have: global size == local size of distribution
    IndexType totalRows = matrix.getDistribution().getGlobalSize();
    std::vector<IndexType> csrXadj( totalRows + 1 ); // #rows + 1
    // std::vector<IndexType> csrAdjncy( matrix.getLocalNumValues() - totalRows ); // #values - #diagonalelements(rows)
    std::vector<IndexType> csrAdjncy( matrix.getCSRGraphSize() );
    std::vector<IndexType> csrVwgt( totalRows );
    matrix.buildCSRGraph( &csrXadj[0], &csrAdjncy[0], &csrVwgt[0], NULL );
    IndexType ncon = 1;
    IndexType options[METIS_NOPTIONS];
    METIS_SetDefaultOptions( options );
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
    // use only with PartGraphKway
    // options[ METIS_OPTION_OBJTYPE ] = METIS_OBJTYPE_VOL;//METIS_OBJTYPE_CUT
    // options[ METIS_OPTION_DBGLVL ] = METIS_DBG_TIME;
    // recursive bisection
    METIS_PartGraphRecursive( &totalRows, &ncon, &csrXadj[0], &csrAdjncy[0], &csrVwgt[0], NULL, NULL, &parts,
                              &tpwgts[0], NULL, options, &minConstraint, &partition[0] );
    // multilevel k-way partitioning (used in ParMetis)
    // METIS_PartGraphKway( &totalRows, &ncon, &csrXadj[0], &csrAdjncy[0], &csrVwgt[0], NULL,
    // NULL, &parts, &tpwgts[0], NULL, options, &minConstraint, &partition[0] );
}

template<typename WeightType>
void MetisPartitioning::checkAndMapWeights(
    std::vector<WeightType>& tpwgts,
    std::vector<IndexType>& mapping,
    IndexType& count,
    std::vector<float>& weights,
    IndexType size ) const
{
    for ( IndexType i = 0; i < size; ++i )
    {
        if ( weights[i] != 0 )
        {
            mapping[count] = i;
            tpwgts[count] = weights[i];
            ++count;
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

} /* end namespace dmemo */

} /* end namespace scai */
