/**
 * @file CSRGraph.cpp
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
 * @brief Implementation of methods for CSR graph
 * @author Thomas Brandes
 * @date 02.01.2018
 */

#include <scai/partitioning/CSRGraph.hpp>

#include <scai/lama/matrix/SparseMatrix.hpp>
#include <scai/lama/io/FileIO.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>

namespace scai
{

using namespace dmemo;
using namespace hmemo;

using lama::CSRSparseMatrix;

namespace partitioning
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename IdxType>, CSRGraph<IdxType>::logger, "CSRGraph" )

/* ---------------------------------------------------------------------- */

template<typename IdxType>
void CSRGraph<IdxType>::joinCSR( 
    HArray<IdxType>& outIA,
    HArray<IdxType>& outJA,
    const HArray<IndexType>& localIA,
    const HArray<IndexType>& localJA,
    const HArray<IndexType>& haloIA,
    const HArray<IndexType>& haloJA )
{
    SCAI_LOG_DEBUG( logger, "join CSR: local with ia = " << localIA << ", ja = " << localJA 
                             << ", halo with ia = " << haloIA << ", ja = " << haloJA )

    SCAI_ASSERT_EQ_ERROR( localIA.size(), haloIA.size(), "serious mismatch for IA arrays of local and halo" )

    IndexType numRows = localIA.size() - 1;

    WriteOnlyAccess<IdxType> ia( outIA, numRows + 1 );

    ReadAccess<IndexType> ia1( localIA );
    ReadAccess<IndexType> ia2( haloIA );
    ReadAccess<IndexType> ja1( localJA );
    ReadAccess<IndexType> ja2( haloJA );

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IdxType nonZeros = 0; // count for each row in parallel

        for ( IndexType jj = ia1[i]; jj < ia1[i+1]; ++jj )
        {
            // CSR data of local might contain invalidIndex 

            if ( ja1[jj] != invalidIndex )
            {
                nonZeros++;
            }
        }

        nonZeros += ia2[i + 1] - ia2[i];

        ia[i] = nonZeros;
    }

    // cannot use kernel sizes2offsets as IdxType might not be IndexType

    IndexType numValues = 0;

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IdxType tmp = numValues;
        numValues += ia[i];
        ia [i] = tmp;
    }

    ia[numRows] = numValues;

    SCAI_LOG_DEBUG( logger, "now build ja, has " << numValues << " entries " )

    WriteOnlyAccess<IdxType> ja( outJA, numValues );

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offset = ia[i];

        for ( IndexType jj = ia1[i]; jj < ia1[i+1]; ++jj )
        {
            if ( ja1[jj] != invalidIndex )
            {
                ja[offset++] = static_cast<IdxType>( ja1[jj] );
            }
        }

        for ( IndexType jj = ia2[i]; jj < ia2[i+1]; ++jj )
        {
            ja[offset++] = static_cast<IdxType>( ja2[jj] );
        }
    }

    SCAI_LOG_DEBUG( logger, "join CSR graph data is ready" )
}

/* ---------------------------------------------------------------------- */

template<typename IdxType>
template<typename ValueType> 
void CSRGraph<IdxType>::buildByCSRSparseMatrix( const CSRSparseMatrix<ValueType>& matrix, bool isSingle )
{
    SCAI_LOG_INFO( logger, "Build CSR graph by this CSR matrix: " << matrix << ", isSingle = " << isSingle )

    SCAI_ASSERT_EQ_ERROR( matrix.getRowDistribution(), matrix.getColDistribution(), 
                          "Not a square matrix with same row/col distribution: " << matrix )

    const Distribution& dist = matrix.getRowDistribution();
    const Communicator& comm = dist.getCommunicator();

    const lama::CSRStorage<ValueType>& localStorage = matrix.getLocalStorage();
    const lama::CSRStorage<ValueType>& haloStorage = matrix.getHaloStorage();
    const Halo& haloSchedule = matrix.getHalo();

    // local2global for local IA, halo2global for halo JA

    HArray<IndexType>& localIA = const_cast<HArray<IndexType>&>( localStorage.getIA() );
    HArray<IndexType>& localJA = const_cast<HArray<IndexType>&>( localStorage.getJA() );

    HArray<IndexType>& haloIA = const_cast<HArray<IndexType>&>( haloStorage.getIA() );
    HArray<IndexType>& haloJA = const_cast<HArray<IndexType>&>( haloStorage.getJA() ); 

    IndexType numLocalRows = dist.getLocalSize();

    SCAI_ASSERT_EQ_DEBUG( numLocalRows + 1, localIA.size(), "illegal offset array for local storage" )
    SCAI_ASSERT_EQ_DEBUG( numLocalRows + 1, haloIA.size(), "illegal offset array for halo storage" )

    // local2global ( localJA ), just done by adding offset, set invalidIndex for diagonal elements

    SCAI_LOG_DEBUG( logger, "update of localJA: mark diagonal elements" )

    IndexType myRank   = comm.getRank();

    IndexType myOffset = 0; 

    if ( !isSingle )
    {
        SCAI_ASSERT_EQ_ERROR( mDistOffsets.size(), comm.getSize() + 1, "distributed offsets must be computed before" )
        auto rOffsets = hostReadAccess( mDistOffsets );
        myOffset = rOffsets[myRank];
    }

    {
        ReadAccess<IndexType> ia( localIA );
        WriteAccess<IndexType> ja( localJA );
 
        for ( IndexType i = 0; i < numLocalRows; ++i )
        {
            for ( IndexType jj = ia[i]; jj < ia[i+1]; jj++ )
            {
                if ( ja[jj] == i )
                {
                     ja[jj] = invalidIndex;   // diagonal element is skiped
                }
                else
                { 
                     ja[jj] += myOffset;
                }
            }
        }
    }

    SCAI_LOG_DEBUG( logger, "update of haloJA: translate halo indexes to global indexes" )

    // Make halo exchange 

    HArray<IndexType> providesGlobalIndexes( haloSchedule.getProvidesIndexes() );

    // add myOffset so the column indexes have the correct numbering
  
    utilskernel::HArrayUtils::binaryOpScalar( providesGlobalIndexes, providesGlobalIndexes, myOffset, common::BinaryOp::ADD, false );

    HArray<IndexType> requiredGlobalIndexes;

    // translate the provides to global indexes by adding myOffset
     
    if ( !haloSchedule.isEmpty() )
    {
        SCAI_LOG_DEBUG( logger, "exchange halo, here global column indexes: send is " << haloSchedule.getProvidesPlan()
                                << ", recv is " << haloSchedule.getRequiredPlan() )

        comm.exchangeByPlan( requiredGlobalIndexes, haloSchedule.getRequiredPlan(), 
                             providesGlobalIndexes, haloSchedule.getProvidesPlan() );
    }

    SCAI_LOG_DEBUG( logger, comm << ": gather translates now" )

    utilskernel::HArrayUtils::gather( haloJA, requiredGlobalIndexes, haloJA, common::BinaryOp::COPY );

    joinCSR( mIA, mJA, localIA, localJA, haloIA, haloJA );

    SCAI_LOG_DEBUG( logger, comm << ": local CSR graph data now complete" )
}

/* ---------------------------------------------------------------------- */

template<typename IdxType>
CSRGraph<IdxType>::CSRGraph( const lama::_Matrix& matrix, bool isSingle )
{
    SCAI_LOG_INFO( logger, "CSRGraph constructor ( isSingle = " << isSingle << " ), matrix = " << matrix )

    // make sure that we have a square matrix

    SCAI_ASSERT_EQ_ERROR( matrix.getNumRows(), matrix.getNumColumns(), "not a square matrix" );
    
    // this is the relevant distribution of how the matrix data is distributed

    const DistributionPtr& dist = matrix.getRowDistributionPtr();

    // determine the distribution offsets 

    if ( !isSingle )
    {
        getDistributionOffsets( mDistOffsets, *dist );
    }

    // we build the CSR graph data from a CSR sparse matrix
    // ToDo: skip this copy if matrix is already sparse matrix with same row/col dist

    CSRSparseMatrix<DefaultReal> csrMatrix;
    SCAI_LOG_INFO( logger, "convert now to csrMatrix" )
    csrMatrix.assign( matrix );
    SCAI_LOG_INFO( logger, "csrMatrix = " << csrMatrix )
    csrMatrix.redistribute( dist, dist );
    SCAI_LOG_INFO( logger, "now redistributed: csrMatrix = " << csrMatrix )

    SCAI_LOG_INFO( logger, "CSRGraph constructor, uses now this csr matrix: " << csrMatrix )

    buildByCSRSparseMatrix( csrMatrix, isSingle );

    SCAI_LOG_INFO( logger, "CSR graph (local) data built, now compue weights" )

    // now set the weights for the vertices

    IndexType nLocal = dist->getLocalSize();

    {
        WriteOnlyAccess<IdxType> wWeight( mVertexWeight, nLocal ); 
        ReadAccess<IdxType> rIA( mIA );
        for ( IndexType i = 0; i < nLocal; ++i )
        {
            wWeight[i] = 1 + ( rIA[i+1] - rIA[i] );
        }
    }

    SCAI_LOG_INFO( logger, "build weights array: " << mVertexWeight )
}

/* ---------------------------------------------------------------------- */

template<typename IdxType>
void CSRGraph<IdxType>::getDistributionOffsets( hmemo::HArray<IdxType>& offsets, const dmemo::Distribution& dist )
{
    const dmemo::Communicator& comm = dist.getCommunicator();

    const PartitionId numPartitions = comm.getSize();
    const PartitionId MASTER = 0;

    const IdxType mySize = dist.getLocalSize();

    {
        hmemo::WriteOnlyAccess<IdxType> wOffsets( offsets, numPartitions + 1 );
        comm.gather( wOffsets.get(), 1, MASTER, &mySize );
        comm.bcast( wOffsets.get(), numPartitions, MASTER );

        IdxType totalSize = 0;
        for ( PartitionId p = 0; p < numPartitions; ++p )
        {
            IdxType tmp = totalSize;
            totalSize += wOffsets[p];
            wOffsets[p] = tmp;
        }

        wOffsets[numPartitions] = totalSize;
    }

    SCAI_LOG_INFO( logger, "computed dist offsets = " << offsets << " for dist = " << dist 
                    << ", comm = " << comm << ", #part = " << numPartitions )
}

/* ---------------------------------------------------------------------- */

template<typename IdxType>
CSRGraph<IdxType>::~CSRGraph()
{
}

/* ====================================================================== */

// Metis and Parmetis come with MPI_INT or MPI_LONG_LONG_INT 

SCAI_COMMON_INST_CLASS( CSRGraph, int )

} /* end namespace partitioning */

} /* end namespace scai */
