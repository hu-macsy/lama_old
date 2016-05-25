/**
 * @file MetisDistribution.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief MetisDistribution.cpp
 * @author Lauretta Schubert
 * @date 01.07.2013
 */

// hpp
#include <scai/dmemo/MetisDistribution.hpp>

// local library
#include <scai/dmemo/NoDistribution.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

extern "C"
{
#include <metis.h>
}

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( MetisDistribution::logger, "Distribution.MetisDistribution" )

#define MASTER 0

MetisDistribution::MetisDistribution( const CommunicatorPtr comm, const Distributed& matrix, std::vector<float>& weights )

    : GeneralDistribution( matrix.getDistribution().getGlobalSize(), comm )
{
    computeIt( comm, matrix, weights );
}

MetisDistribution::MetisDistribution( const CommunicatorPtr comm, const Distributed& matrix, float weight )

    : GeneralDistribution( matrix.getDistribution().getGlobalSize(), comm )
{
    SCAI_LOG_INFO( logger, "construct Metis distribution, weight at " << *comm << ": " << weight )

    // collect weights from all processors

    PartitionId numPartitions = comm->getSize();

    std::vector<float> weights( numPartitions );

    SCAI_LOG_INFO( logger, "#weights = " << weights.size() )

    comm->gather( &weights[0], 1, MASTER, &weight );
    comm->bcast( &weights[0], numPartitions, MASTER );

    // norm weights so that sum gives exactly 1.0

    normWeights( weights );

    SCAI_LOG_INFO( logger, "#weights = " << weights.size() )

    for ( size_t i = 0; i < weights.size(); ++i )
    {
        SCAI_LOG_INFO( logger, "weight[" << i << "] = " << weights[i] )
    }

    computeIt( comm, matrix, weights );
}

void MetisDistribution::normWeights( std::vector<float>& weights )
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

void MetisDistribution::computeIt( const CommunicatorPtr comm, const Distributed& matrix, std::vector<float>& weights )
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
                           "MetisDistribution called with 1 processor/1 weight, which is the same as NoDistribution." )

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

    // scatter local rows
    int numMyRows = 0;
    comm->scatter( &numMyRows, 1, MASTER, &numRowsPerOwner[0] );
    mLocal2Global.resize( numMyRows );
    comm->scatterV( &mLocal2Global[0], numMyRows, MASTER, &rows[0], &numRowsPerOwner[0] );

    std::vector<IndexType>::const_iterator end = mLocal2Global.end();
    std::vector<IndexType>::const_iterator begin = mLocal2Global.begin();

    for ( std::vector<IndexType>::const_iterator it = begin; it != end; ++it )
    {
        IndexType i = static_cast<IndexType>( std::distance( begin, it ) );
        SCAI_ASSERT( 0 <= *it && *it < mGlobalSize,
                     *it << " is illegal index for general distribution of size " << mGlobalSize )
        mGlobal2Local[ *it] = i;
    }
}

MetisDistribution::~MetisDistribution()
{
    SCAI_LOG_INFO( logger, "~MetisDistribution" )
}

void MetisDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "MetisDistribution( size = " << mLocal2Global.size() << " of " << mGlobalSize << ", comm = "
           << *mCommunicator << " )";
}

template<typename WeightType>
void MetisDistribution::callPartitioning(
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

/*  Help routine needed for parallel graphs

    void getVTX( IndexType* vtxdist, CommunicatorPtr com, const Distribted& matrix )
    {
        if ( vtxdist == NULL )
        {
            return;
        }

        const PartitionId MASTER = 0;

        int numLocalRows = matrix->getDistribution().getLocalSize();

        IndexType parts = comm->getSize();

        // Is this valid ?
        // SCAI_ASSERT_ERROR( getDistribution().getNumPartitions() == parts,
        //              "mismatch number of partitions and communicator size" );

        comm->gather( vtxdist, 1, MASTER, &numLocalRows );
        comm->bcast( vtxdist, parts, MASTER );

        // build running sum

        IndexType runningSum = 0;

        for ( IndexType i = 0; i < parts; ++i )
        {
            IndexType tmp = runningSum;
            runningSum += vtxdist[i];
            vtxdist[i] = tmp;
        }
        vtxdist[parts] = runningSum;
    }

*/

template<typename WeightType>
void MetisDistribution::checkAndMapWeights(
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
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

const char MetisDistribution::theCreateValue[] = "METIS";

std::string MetisDistribution::createValue()
{
    return theCreateValue;
}

Distribution* MetisDistribution::create( const DistributionArguments arg )
{
    if ( arg.matrix != NULL )
    {
        return new MetisDistribution( arg.communicator, *arg.matrix, arg.weight );
    }
    else
    {
        DistributionPtr dist( new NoDistribution( arg.globalSize ) );

        Distributed dummy( dist );

        return new MetisDistribution( arg.communicator, dummy, arg.weight );
    }
}

} /* end namespace dmemo */

} /* end namespace scai */
