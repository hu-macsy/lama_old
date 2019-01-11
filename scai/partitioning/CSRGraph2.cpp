/**
 * @file CSRGraph2.cpp
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
 * @brief Implementation of methods for CSR2 graph
 * @author Thomas Brandes
 * @date 02.01.2018
 */

#include <scai/partitioning/CSRGraph2.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using namespace lama;
using namespace dmemo;
using namespace hmemo;

namespace partitioning
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename IdxType>, CSRGraph2<IdxType>::logger, "CSRGraph2" )

template<typename IdxType>
CSRGraph2<IdxType>::CSRGraph2( const lama::_Matrix& matrix )
{
    // make sure that we have a square matrix

    SCAI_LOG_INFO( logger, "build CSR bipartite graph for matrix = " << matrix )

    DistributionPtr rowDist( new NoDistribution( matrix.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( matrix.getNumColumns() ) );

    typedef DefaultReal ValueType;

    CSRSparseMatrix<ValueType> csrMatrix;
    csrMatrix.assign( matrix );
    csrMatrix.redistribute( rowDist, colDist );
    CSRStorage<ValueType>& storage = csrMatrix.getLocalStorage();
    CSRStorage<ValueType> storageT;
    storageT.assignTranspose( storage );

    SCAI_LOG_INFO( logger, "CSRGraph2: use CSR storage " << storage << " and the transposed " << storageT )

    IndexType numRows    = storage.getNumRows();
    IndexType numColumns = storage.getNumColumns();
    IndexType numValues  = storage.getNumValues();

    // now build the bipartite CSR graph of it

    IndexType nNodes = numRows + numColumns;
    IndexType nEdges = 2 * numValues;

    WriteOnlyAccess<IdxType> wIA( mIA, nNodes + 1 );
    WriteOnlyAccess<IdxType> wJA( mJA, nEdges );
    WriteOnlyAccess<IdxType> wVertexWeights( mVertexWeights, mNumConstraints *  nNodes );

    ReadAccess<IndexType> csrIA( storage.getIA() );
    ReadAccess<IndexType> csrJA( storage.getJA() );
    ReadAccess<IndexType> csrTIA( storageT.getIA() );
    ReadAccess<IndexType> csrTJA( storageT.getJA() );

    SCAI_LOG_INFO( logger, "CSRGraph2: build #nodes = " << nNodes << ", #edges = " << nEdges )

    IndexType offset = 0;

    for ( IndexType i = 0; i < nNodes; ++i )
    {
        if ( i < numRows )
        {
            wIA[i] = offset;

            wVertexWeights[mNumConstraints * i] = csrIA[i+1] - csrIA[i] + 5;

            if ( mNumConstraints == 2 ) wVertexWeights[mNumConstraints * i + 1] = 1;

            if ( mNumConstraints == 3 ) 
            {
                wVertexWeights[mNumConstraints * i + 1] = 100;
                wVertexWeights[mNumConstraints * i + 2] = 1;
            }

            for ( IndexType jj = csrIA[i]; jj < csrIA[i+1]; ++jj )
            {
                // edge from row i to column csrJA[jj]

                SCAI_ASSERT_VALID_INDEX_DEBUG( offset, nEdges, "out of range, i = " << i )
                wJA[offset++] = numRows + csrJA[jj];
            }
        }
        else
        {
            IndexType j = i - numRows;

            wIA[i] = offset;
 
            wVertexWeights[mNumConstraints * i] = csrTIA[ j + 1 ] - csrTIA[ j ];

            if ( mNumConstraints == 2 ) wVertexWeights[mNumConstraints * i + 1] = 1;

            if ( mNumConstraints == 3 ) 
            {
                wVertexWeights[mNumConstraints * i + 1] = 1;
                wVertexWeights[mNumConstraints * i + 2] = 100;
            }

            for ( IndexType ii = csrTIA[j]; ii < csrTIA[j+1]; ++ii )
            {
                // edge from col j to row csrTJA[ii]

                SCAI_ASSERT_VALID_INDEX_DEBUG( offset, nEdges, "out of range, row j = " << j )
                wJA[offset++] = csrTJA[ii];
            }
        }
    }

    wIA[ nNodes ] = offset;

    SCAI_ASSERT_EQ_ERROR( offset, nEdges, "serious mismatch" )

    SCAI_LOG_INFO( logger, "Build csr bipartite graph done." )
}

template<typename IdxType>
CSRGraph2<IdxType>::~CSRGraph2()
{
}

// Metis and Parmetis come with MPI_INT or MPI_LONG_LONG_INT 

SCAI_COMMON_INST_CLASS( CSRGraph2, int )

} /* end namespace partitioning */

} /* end namespace scai */
