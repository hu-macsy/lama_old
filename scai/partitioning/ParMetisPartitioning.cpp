/**
 * @file ParMetisPartitioning.cpp
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
 * @brief ParMetisPartitioning.cpp
 * @author ThomasBrandes
 * @date 17.08.2017
 */

// hpp
#include <scai/partitioning/ParMetisPartitioning.hpp>
#include <scai/partitioning/CSRGraph.hpp>

// local library
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/mpi/MPICommunicator.hpp>

// internal scai libraries
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/tracing.hpp>

#include <parmetis.h>

namespace scai
{

using namespace dmemo;
using namespace hmemo;
using namespace lama;

namespace partitioning
{

SCAI_LOG_DEF_LOGGER( ParMetisPartitioning::logger, "Partitioning.ParMetisPartitioning" )

#define MASTER 0

ParMetisPartitioning::ParMetisPartitioning()
{
}

void ParMetisPartitioning::writeAt( std::ostream& stream ) const
{
    stream << "ParMetisPartitioning";
}

ParMetisPartitioning::~ParMetisPartitioning()
{
    SCAI_LOG_INFO( logger, "~ParMetisPartitioning" )
}

/* ---------------------------------------------------------------------------------*/

void ParMetisPartitioning::rectangularPartitioning(
    hmemo::HArray<PartitionId>& /* rowMapping */,
    hmemo::HArray<PartitionId>& /* colMapping */,
    const lama::_Matrix& /* matrix     */,
    const hmemo::HArray<float>& /* processorWeights */ ) const
{
    COMMON_THROWEXCEPTION( "Not available yet" )
}

/* ---------------------------------------------------------------------------------*/

void ParMetisPartitioning::squarePartitioning( 
    HArray<PartitionId>& newLocalOwners,
    const lama::_Matrix& matrix,
    const HArray<float>& processorWeights ) const
{
    SCAI_REGION( "ParMETIScall" )

    auto commPtr = Communicator::getCommunicatorPtr();

    const Distribution& dist = matrix.getRowDistribution();
    const Communicator& comm = *commPtr;
  
    SCAI_LOG_INFO( logger, comm << ": square partitioning with ParMetis" )

    // Note that here we have: global size == local size of distribution

    // use only with PartGraphKway
    // options[ METIS_OPTION_OBJTYPE ] = METIS_OBJTYPE_VOL;//METIS_OBJTYPE_CUT
    // options[ METIS_OPTION_DBGLVL ] = METIS_DBG_TIME;
    // recursive bisection

    SCAI_ASSERT_EQ_ERROR( comm.getType(), Communicator::MPI, "ParMetis only works on MPI distributed data" )

    const MPICommunicator& mpiComm = static_cast<const MPICommunicator&>( comm );

    MPI_Comm lamaComm = mpiComm.getMPIComm();

    idx_t wgtflag    = 0;
    idx_t numflag    = 0;  // C-style numbering, indexes start with 0, is default in LAMA
    idx_t ncon       = 1;  // number of weights per edge
    idx_t nparts     = processorWeights.size();   // desired number of processors
    real_t ubvec[]   = { 1.05 };  // ncon entries, imbalance allowed
    idx_t  options[] = { 0 };   // take default options
    idx_t  edgeCut   = 0;         // is output argument
   
    SCAI_LOG_DEBUG( logger, comm << ": processor weights, same values on all processors" )

    HArray<real_t> tpwghts;
    utilskernel::HArrayUtils::assign( tpwghts, processorWeights );
    WriteAccess<real_t> rTPWeights( tpwghts );       // take write access, otherwise const cast required

    IndexType nlocal = dist.getLocalSize();

    HArray<idx_t> partition; 
    WriteOnlyAccess<idx_t> wPartition( partition, nlocal );   // take write access, otherwise const cast required

    CSRGraph<idx_t> graph( matrix );

    SCAI_LOG_INFO( logger, comm << ": CSR graph ready" )

    ReadAccess<idx_t> rDistOffsets( graph.distOffsets() );
    ReadAccess<idx_t> rIA( graph.ia() );
    ReadAccess<idx_t> rJA( graph.ja() );
    ReadAccess<idx_t> rWeights( graph.weights() );

    idx_t res = ParMETIS_V3_PartKway( const_cast<idx_t*>( rDistOffsets.get() ), 
                                      const_cast<idx_t*>( rIA.get() ),
                                      const_cast<idx_t*>( rJA.get() ),
                                      const_cast<idx_t*>( rWeights.get() ), 
                                      NULL, 
                                      &wgtflag,
                                      &numflag, &ncon, &nparts,
                                      rTPWeights.get(), 
                                      ubvec, 
                                      options, 
                                      &edgeCut, 
                                      wPartition.get(),
                                      &lamaComm );

    SCAI_LOG_INFO( logger, comm << ": parmetis result = " << res << ", edgecut = " << edgeCut )

    wPartition.release();

    if ( newLocalOwners.getValueType() == partition.getValueType() )
    {
        SCAI_LOG_INFO( logger, "move partitionining results to new local owners" )
        HArray<IndexType>& partitionT = reinterpret_cast<HArray<IndexType>&>( partition );
        newLocalOwners = std::move( partitionT );
    }
    else
    {
        SCAI_LOG_INFO( logger, "assign partitionining results to new local owners" )
        utilskernel::HArrayUtils::assign( newLocalOwners, partition );
    }
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in Partitioning factory )    *
 * ---------------------------------------------------------------------------------*/

const char ParMetisPartitioning::theCreateValue[] = "PARMETIS";

std::string ParMetisPartitioning::createValue()
{
    return theCreateValue;
}

PartitioningPtr ParMetisPartitioning::create()
{
    return PartitioningPtr( new ParMetisPartitioning() );
}

} /* end namespace partitioning */

} /* end namespace scai */
