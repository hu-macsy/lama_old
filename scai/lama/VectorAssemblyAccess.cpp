/**
 * @file VectorAssemblyAccess.cpp
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

#include <scai/lama/VectorAssemblyAccess.hpp>

#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using namespace hmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, VectorAssemblyAccess<ValueType>::logger,
                              "VectorAssemblyAccess" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
VectorAssemblyAccess<ValueType>::VectorAssemblyAccess( _Vector& vector, const common::binary::BinaryOp op ) : 

    mVector( vector ),
    mIsReleased( false ),
    mOp( op )
{
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssemblyAccess<ValueType>::exchangeCOO( 
    HArray<IndexType>& outIA,
    HArray<ValueType>& outValues,
    const HArray<IndexType> inIA,
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
    HArray<ValueType> sendValues;

    HArrayUtils::gather( sendIA, inIA, perm, common::binary::COPY );
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
    comm.exchangeByPlan( outValues, recvPlan, sendValues, sendPlan );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssemblyAccess<ValueType>::shiftAssembledData(
    const HArray<IndexType>& myIA, 
    const HArray<ValueType>& myValues )
{
    // This method shifts all assembled data and each processor applies a routine on it

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    PartitionId np = comm->getSize();

    mVector.fillSparseData( myIA, myValues, mOp );

    if ( np == 1 )
    {
        return;
    }

    const int COMM_DIRECTION = 1;  // circular shifting from left to right

    // determine the maximal size of assembled data for good allocation of buffers

    IndexType maxSize = comm->max( myIA.size() );

    ContextPtr contextPtr = Context::getHostPtr();

    HArray<IndexType> sendIA;
    HArray<ValueType> sendValues;
    HArray<IndexType> recvIA;
    HArray<ValueType> recvValues;

    sendIA.reserve( contextPtr, maxSize );
    sendValues.reserve( contextPtr, maxSize );
    recvIA.reserve( contextPtr, maxSize );
    recvValues.reserve( contextPtr, maxSize );

    // np - 1 shift steps are neeed

    for ( PartitionId p = 0; p < np - 1; ++p )
    {
        if ( p == 0 )
        {
            comm->shiftArray( recvIA, myIA, COMM_DIRECTION );
            comm->shiftArray( recvValues, myValues, COMM_DIRECTION );
        }
        else
        {
            comm->shiftArray( recvIA, sendIA, COMM_DIRECTION );
            comm->shiftArray( recvValues, sendValues, COMM_DIRECTION );
        }

        mVector.fillSparseData( recvIA, recvValues, mOp );
        
        // prepare for next step, the received values from left will be sent to right

        std::swap( sendValues, recvValues );
        std::swap( sendIA, recvIA );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssemblyAccess<ValueType>::release()
{
    SCAI_ASSERT_EQ_DEBUG( mIA.size(), mValues.size(), "serious mismatch" );

    // Attention: even if mIA.size() == 0, this processor must participate in communication

    if ( mIsReleased )
    {
        return;
    }

    // vector data only read, so we can use HArray references

    HArrayRef<IndexType> ia( mIA );
    HArrayRef<ValueType> values( mValues );

    // These COO array will keep only the values owned by this processor


    const dmemo::Distribution& dist = mVector.getDistribution();

    if ( dist.isReplicated() )
    {
        shiftAssembledData( ia, values );
    }
    else
    {
        HArray<IndexType> ownedIA;
        HArray<ValueType> ownedValues;

        exchangeCOO( ownedIA, ownedValues, ia, values, dist );

        dist.global2local( ownedIA );   // translate global indexes to local indexes

        // now we add the owned COO data to the local vector data

        mVector.fillSparseData( ownedIA, ownedValues, mOp );
    }

    // reset the data vectors as they are emptied now

    mIA.clear();
    mValues.clear();

    HArrayRef<IndexType> l_ia( mLocalIA );
    HArrayRef<ValueType> l_values( mLocalValues );
    mVector.fillSparseData( l_ia, l_values, mOp );

    mLocalIA.clear();
    mLocalValues.clear();

    mIsReleased = true;
}

SCAI_COMMON_INST_CLASS( VectorAssemblyAccess, SCAI_ARRAY_TYPES_HOST )

}

}
