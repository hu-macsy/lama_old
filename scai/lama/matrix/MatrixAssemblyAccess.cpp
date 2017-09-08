/**
 * @file MatrixAssemblyAccess.cpp
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
 * @brief Access to a matrix to add matrix elements
 * @author Thomas Brandes
 * @date 07.09.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/MatrixAssemblyAccess.hpp>

#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using namespace hmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, MatrixAssemblyAccess<ValueType>::logger,
                              "MatrixAssemblyAccess" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
MatrixAssemblyAccess<ValueType>::MatrixAssemblyAccess( Matrix& matrix ) : mMatrix( matrix )
{
    SCAI_ASSERT_EQ_ERROR( matrix.getMatrixKind(), Matrix::SPARSE, "Assembly only for sparse matrix supported" )
    SCAI_ASSERT_EQ_ERROR( matrix.getNumValues(), 0, "Assembly only for zero sparse matrices supported" )
 
    mIsReleased = false;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssemblyAccess<ValueType>::global2local( 
    hmemo::HArray<IndexType>& ia,
    const dmemo::Distribution& dist )
{
    IndexType nLocalRows = dist.getLocalSize();
    IndexType nnz = ia.size();

    {
        WriteAccess<IndexType> wIA( ia );
        for ( IndexType i = 0; i < nnz; ++i )
        {
            wIA[i] = dist.global2local( wIA[i] );
            SCAI_ASSERT_VALID_INDEX_DEBUG( wIA[i], nLocalRows, "illegal row index" )
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssemblyAccess<ValueType>::exchangeCOO( 
    HArray<IndexType>& outIA,
    HArray<IndexType>& outJA,
    HArray<ValueType>& outValues,
    const HArray<IndexType> inIA,
    const HArray<IndexType> inJA,
    const HArray<ValueType> inValues,
    const dmemo::Distribution& dist )
{
    using namespace utilskernel;

    HArray<PartitionId> owners;

    dist.computeOwners( owners, inIA );

    SCAI_LOG_DEBUG( logger, "owners = " << owners )

    const dmemo::Communicator& comm = dist.getCommunicator();
    PartitionId np = comm.getSize();

    HArray<IndexType> perm;
    HArray<IndexType> offsets;

    HArrayUtils::bucketSort( offsets, perm, owners, np );

    SCAI_LOG_DEBUG( logger, "sorted, perm = " << perm << ", offsets = " << offsets )

    HArray<IndexType> sendIA;
    HArray<IndexType> sendJA;
    HArray<ValueType> sendValues;

    HArrayUtils::gather( sendIA, inIA, perm, common::binary::COPY );
    HArrayUtils::gather( sendJA, inJA, perm, common::binary::COPY );
    HArrayUtils::gather( sendValues, inValues, perm, common::binary::COPY );

    HArrayUtils::unscan( offsets );  // now we have size

    SCAI_LOG_DEBUG( logger, "sizes = " << offsets )

    dmemo::CommunicationPlan sendPlan;
    dmemo::CommunicationPlan recvPlan;

    {
        ReadAccess<IndexType> rSizes( offsets );
        sendPlan.allocate( rSizes.get(), np );
    }

    recvPlan.allocateTranspose( sendPlan, comm );

    SCAI_LOG_DEBUG( logger, "recv plan: " << recvPlan )

    comm.exchangeByPlan( outIA, recvPlan, sendIA, sendPlan );
    comm.exchangeByPlan( outJA, recvPlan, sendJA, sendPlan );
    comm.exchangeByPlan( outValues, recvPlan, sendValues, sendPlan );
}

template<typename ValueType>
void MatrixAssemblyAccess<ValueType>::release()
{
    SCAI_ASSERT_EQ_DEBUG( mIA.size(), mJA.size(), "serious mismatch" )
    SCAI_ASSERT_EQ_DEBUG( mIA.size(), mValues.size(), "serious mismatch" );

    // Attention: even if mIA.size() == 0, this processor must participate in communication

    if ( mIsReleased )
    {
        return;
    }

    // vector data only read, so we can use HArray references

    HArrayRef<IndexType> ia( mIA );
    HArrayRef<IndexType> ja( mJA );
    HArrayRef<ValueType> values( mValues );

    // These COO array will keep only the values owned by this processor

    HArray<IndexType> ownedIA;
    HArray<IndexType> ownedJA;
    HArray<ValueType> ownedValues;

    const dmemo::Distribution& rowDist = mMatrix.getRowDistribution();

    exchangeCOO( ownedIA, ownedJA, ownedValues, ia, ja, values, rowDist );

    global2local( ownedIA, rowDist );

    // now we add the owned COO data to the local storage

    COOStorage<ValueType> cooLocal;

    cooLocal.allocate( rowDist.getLocalSize(), mMatrix.getNumColumns() );
    cooLocal.swap( ownedIA, ownedJA, ownedValues );

    mMatrix.assign( cooLocal, mMatrix.getRowDistributionPtr(), mMatrix.getColDistributionPtr() );

    // reset the data vectors as they are emptied now

    mIA.clear();
    mJA.clear();
    mValues.clear();

    mIsReleased = true;
}

SCAI_COMMON_INST_CLASS( MatrixAssemblyAccess, SCAI_NUMERIC_TYPES_HOST )

}

}
