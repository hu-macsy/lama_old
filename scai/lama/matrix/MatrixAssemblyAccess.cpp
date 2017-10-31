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
MatrixAssemblyAccess<ValueType>::MatrixAssemblyAccess( _Matrix& matrix, const common::binary::BinaryOp op ) : 

    mMatrix( matrix ),
    mIsReleased( false ),
    mOp( op )
{
    SCAI_ASSERT_EQ_ERROR( matrix.getMatrixKind(), MatrixKind::SPARSE, "Assembly only for sparse matrix supported" )

    // SCAI_ASSERT_EQ_ERROR( matrix.getNumValues(), 0, "Assembly only for zero sparse matrices supported" )
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

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssemblyAccess<ValueType>::shiftAssembledData(
    CSRStorage<ValueType>& localStorage,
    const HArray<IndexType>& myIA,
    const HArray<IndexType>& myJA,
    const HArray<ValueType>& myValues )
{   
    // This method shifts all assembled data and each processor applies a routine on it

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    PartitionId np = comm->getSize();

    SCAI_LOG_INFO( logger, "shift assembled data circular through all processors, "
                           << ", me = " << *comm << " have " << myIA.size() << " entries" )

    COOStorage<ValueType> newCOO( localStorage.getNumRows(), localStorage.getNumColumns(), myIA, myJA, myValues );
    CSRStorage<ValueType> newCSR( newCOO );

    localStorage.binaryOpCSR( localStorage, newCSR, mOp );

    if ( np == 1 )
    {
        return;
    }

    const int COMM_DIRECTION = 1;  // circular shifting from left to right

    // determine the maximal size of assembled data for good allocation of buffers

    IndexType maxSize = comm->max( myIA.size() );

    ContextPtr contextPtr = Context::getHostPtr();

    HArray<IndexType> sendIA;
    HArray<IndexType> sendJA;
    HArray<ValueType> sendValues;
    HArray<IndexType> recvIA;
    HArray<IndexType> recvJA;
    HArray<ValueType> recvValues;

    sendIA.reserve( contextPtr, maxSize );
    sendJA.reserve( contextPtr, maxSize );
    sendValues.reserve( contextPtr, maxSize );
    recvIA.reserve( contextPtr, maxSize );
    recvJA.reserve( contextPtr, maxSize );
    recvValues.reserve( contextPtr, maxSize );

    // np - 1 shift steps are neeed

    for ( PartitionId p = 0; p < np - 1; ++p )
    {
        if ( p == 0 )
        {
            comm->shiftArray( recvIA, myIA, COMM_DIRECTION );
            comm->shiftArray( recvJA, myJA, COMM_DIRECTION );
            comm->shiftArray( recvValues, myValues, COMM_DIRECTION );
        }
        else
        {
            comm->shiftArray( recvIA, sendIA, COMM_DIRECTION );
            comm->shiftArray( recvJA, sendJA, COMM_DIRECTION );
            comm->shiftArray( recvValues, sendValues, COMM_DIRECTION );
        }

        newCOO.swap( recvIA, recvJA, recvValues );

        newCSR = newCOO;
        localStorage.binaryOpCSR( localStorage, newCSR, mOp );

        newCOO.swap( recvIA, recvJA, recvValues );

        // prepare for next step, the received values from left will be sent to right

        sendIA.swap( recvIA );
        sendJA.swap( recvJA );
        sendValues.swap( recvValues );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssemblyAccess<ValueType>::addLocalCOO( CSRStorage<ValueType>& localStorage )
{
    SCAI_LOG_INFO( logger, "add " << mLocalIA.size() << " locally assembled entries" )

    // Build CSR storage from the local COO data

    HArrayRef<IndexType> l_ia( mLocalIA );
    HArrayRef<IndexType> l_ja( mLocalJA );
    HArrayRef<ValueType> l_values( mLocalValues );

    COOStorage<ValueType> newCOO( localStorage.getNumRows(), localStorage.getNumColumns(), l_ia, l_ja, l_values );
    CSRStorage<ValueType> newCSR( newCOO );

    localStorage.binaryOpCSR( localStorage, newCSR, mOp );

    mLocalIA.clear();
    mLocalJA.clear();
    mLocalValues.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssemblyAccess<ValueType>::addCOO( CSRStorage<ValueType>& localStorage )
{
    SCAI_LOG_INFO( logger, "add " << mIA.size() << " assembled entries" )

    // vector data only read, so we can use HArray references

    HArrayRef<IndexType> ia( mIA );
    HArrayRef<IndexType> ja( mJA );
    HArrayRef<ValueType> values( mValues );

    const dmemo::Distribution& rowDist = mMatrix.getRowDistribution();

    if ( rowDist.isReplicated() )
    {
        shiftAssembledData( localStorage, ia, ja, values );
    }
    else
    {
        // These COO array will keep only the values owned by this processor

        HArray<IndexType> ownedIA;
        HArray<IndexType> ownedJA;
        HArray<ValueType> ownedValues;

        exchangeCOO( ownedIA, ownedJA, ownedValues, ia, ja, values, rowDist );

        rowDist.global2local( ownedIA );

        // now we add the owned COO data to the local storage

        COOStorage<ValueType> cooLocal;
        cooLocal.allocate( rowDist.getLocalSize(), mMatrix.getNumColumns() );
        cooLocal.swap( ownedIA, ownedJA, ownedValues );

        CSRStorage<ValueType> csrLocal( cooLocal );  // resorts also the entries corresponding to the rows

        localStorage.binaryOpCSR( localStorage, csrLocal, mOp );
    }

    // reset the data vectors as they are emptied now

    mIA.clear();
    mJA.clear();
    mValues.clear();
}

/* -------------------------------------------------------------------------- */

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

    dmemo::DistributionPtr saveColDist = mMatrix.getColDistributionPtr();

    // adding matrix data is only possible on replicated local storage as halo must be rebuilt

    if ( !saveColDist->isReplicated() )
    {
        dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( mMatrix.getNumColumns() ) );
        mMatrix.redistribute( mMatrix.getRowDistributionPtr(), repColDist );
    }

    CSRStorage<ValueType> matrixCSR( mMatrix.getLocalStorage() );

    addCOO( matrixCSR );

    addLocalCOO( matrixCSR );

    SCAI_LOG_DEBUG( logger, "merged CSR = " << matrixCSR )

    mMatrix.assign( matrixCSR, mMatrix.getRowDistributionPtr(), saveColDist );

    mIsReleased = true;
}

SCAI_COMMON_INST_CLASS( MatrixAssemblyAccess, SCAI_NUMERIC_TYPES_HOST )

}

}
