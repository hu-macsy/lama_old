/**
 * @file MatrixAssembly.cpp
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
 * @brief  to a matrix to add matrix elements
 * @author Thomas Brandes
 * @date 07.09.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/MatrixAssembly.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/sparsekernel/COOUtils.hpp>

#include <scai/common/macros/instantiate.hpp>

#include <algorithm>

namespace scai
{

using namespace hmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, MatrixAssembly<ValueType>::logger,
                              "MatrixAssembly" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
MatrixAssembly<ValueType>::MatrixAssembly( dmemo::CommunicatorPtr comm ) 
{
    mComm = comm;
}

template<typename ValueType>
MatrixAssembly<ValueType>::~MatrixAssembly()
{
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssembly<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "MatrixAssembly<" << common::TypeTraits<ValueType>::id() 
           << ">( " << *mComm << ": " << mIA.size() << " entries )";
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssembly<ValueType>::reserve( const IndexType n )
{
    mIA.reserve( n );
    mJA.reserve( n );
    mValues.reserve( n );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType MatrixAssembly<ValueType>::getNumRows() const
{
    auto it = std::max_element( std::begin(mIA), std::end(mIA) ); 

    IndexType numRows = it == end( mIA ) ? 0 : *it;

    return mComm->max( numRows ) + 1;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType MatrixAssembly<ValueType>::getNumColumns() const
{
    auto it = std::max_element( std::begin(mJA), std::end(mJA) ); 

    IndexType numCols = it == end( mJA ) ? 0 : *it;

    return mComm->max( numCols ) + 1;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType MatrixAssembly<ValueType>::getNumValues() const
{
    return mComm->sum( mIA.size() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssembly<ValueType>::exchangeCOO( 
    HArray<IndexType>& outIA,
    HArray<IndexType>& outJA,
    HArray<ValueType>& outValues,
    const HArray<IndexType> inIA,
    const HArray<IndexType> inJA,
    const HArray<ValueType> inValues,
    const dmemo::Distribution& dist ) const
{
    using namespace utilskernel;

    HArray<PartitionId> owners;

    dist.computeOwners( owners, inIA );

    const dmemo::Communicator& comm = dist.getCommunicator();

    SCAI_LOG_DEBUG( logger, comm << ": owners = " << owners )

    PartitionId np = comm.getSize();

    HArray<IndexType> perm;
    HArray<IndexType> offsets;

    HArrayUtils::bucketSort( offsets, perm, owners, np );

    SCAI_LOG_DEBUG( logger, comm << ": sorted, perm = " << perm << ", offsets = " << offsets )

    HArray<IndexType> sendIA;
    HArray<IndexType> sendJA;
    HArray<ValueType> sendValues;

    HArrayUtils::gather( sendIA, inIA, perm, common::BinaryOp::COPY );
    HArrayUtils::gather( sendJA, inJA, perm, common::BinaryOp::COPY );
    HArrayUtils::gather( sendValues, inValues, perm, common::BinaryOp::COPY );

    auto sendPlan = dmemo::CommunicationPlan::buildByOffsets( hostReadAccess( offsets ).get(), np );
    auto recvPlan = sendPlan.transpose( comm );

    SCAI_LOG_DEBUG( logger, comm << ": send plan: " << sendPlan )
    SCAI_LOG_DEBUG( logger, comm << ": recv plan: " << recvPlan )

    comm.exchangeByPlan( outIA, recvPlan, sendIA, sendPlan );
    comm.exchangeByPlan( outJA, recvPlan, sendJA, sendPlan );
    comm.exchangeByPlan( outValues, recvPlan, sendValues, sendPlan );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssembly<ValueType>::checkLegalIndexes( const IndexType numRows, const IndexType numColumns ) const
{
    using namespace utilskernel;

    HArrayRef<IndexType> ia( mIA );
    HArrayRef<IndexType> ja( mJA );
    HArrayRef<ValueType> values( mValues );

    bool okay1 = HArrayUtils::validIndexes( ia, numRows );
    bool okay2 = HArrayUtils::validIndexes( ja, numColumns );

    if ( !okay1 )
    {
        SCAI_LOG_ERROR( logger, *this << " has illegal row index, #rows =  " << numRows );
    }

    if ( !okay2 )
    {
        SCAI_LOG_ERROR( logger, *this << " has illegal col index, #cols = " << numColumns );
    }

    bool okay = mComm->all( okay1 && okay2 );

    if ( !okay )
    {
        COMMON_THROWEXCEPTION( "MatrixAssembly has entries out of range "
                                << numRows << " x " << numColumns );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType> MatrixAssembly<ValueType>::buildLocalCOO( 
    const dmemo::Distribution& dist, 
    const IndexType numColumns,
    common::BinaryOp op ) const
{
    SCAI_ASSERT_EQ_ERROR( dist.getCommunicator(), *mComm, "dist has illegal communicator" )

    checkLegalIndexes( dist.getGlobalSize(), numColumns );

    HArrayRef<IndexType> ia( mIA );
    HArrayRef<IndexType> ja( mJA );
    HArrayRef<ValueType> values( mValues );

    // These COO arrays will keep the matrix items owned by this processor

    HArray<IndexType> ownedIA;
    HArray<IndexType> ownedJA;
    HArray<ValueType> ownedValues;

    exchangeCOO( ownedIA, ownedJA, ownedValues, ia, ja, values, dist );

    dist.global2localV( ownedIA, ownedIA );   // translates global row indexes to local ones

    sparsekernel::COOUtils::normalize( ownedIA, ownedJA, ownedValues, op, Context::getHostPtr() );

    return COOStorage<ValueType>( dist.getLocalSize(), numColumns,
                                  std::move( ownedIA ), std::move( ownedJA ), std::move( ownedValues ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType> MatrixAssembly<ValueType>::buildGlobalCOO( 
    const IndexType numRows,
    const IndexType numColumns,
    common::BinaryOp op ) const
{
    using utilskernel::HArrayUtils;

    ContextPtr contextPtr = Context::getHostPtr();

    checkLegalIndexes( numRows, numColumns );

    const IndexType numValues = mComm->sum( mIA.size() );
    const IndexType maxSize   = mComm->max( mIA.size() );

    // These COO arrays will keep all entries; reserve numValues entries to avoid reallocations

    HArray<IndexType> allIA;
    HArray<IndexType> allJA;
    HArray<ValueType> allValues;

    allIA.reserve( contextPtr, numValues );
    allJA.reserve( contextPtr, numValues );
    allValues.reserve( contextPtr, numValues );

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

    HArrayUtils::appendArray( allIA, HArrayRef<IndexType>( mIA ) );
    HArrayUtils::appendArray( allJA, HArrayRef<IndexType>( mJA ) );
    HArrayUtils::appendArray( allValues, HArrayRef<ValueType>( mValues ) );

    // np - 1 shift steps are neeed

    const PartitionId np = mComm->getSize();

    const int COMM_DIRECTION = 1;  // circular shifting from left to right

    for ( PartitionId p = 0; p < np - 1; ++p )
    {   
        if ( p == 0 )
        {   
            mComm->shiftArray( recvIA, allIA, COMM_DIRECTION );
            mComm->shiftArray( recvJA, allJA, COMM_DIRECTION );
            mComm->shiftArray( recvValues, allValues, COMM_DIRECTION );
        }
        else
        {
            mComm->shiftArray( recvIA, sendIA, COMM_DIRECTION );
            mComm->shiftArray( recvJA, sendJA, COMM_DIRECTION );
            mComm->shiftArray( recvValues, sendValues, COMM_DIRECTION );
        }
        
        HArrayUtils::appendArray( allIA, recvIA );
        HArrayUtils::appendArray( allJA, recvJA );
        HArrayUtils::appendArray( allValues, recvValues );

        // prepare for next step, the received values from left will be sent to right

        sendIA.swap( recvIA );
        sendJA.swap( recvJA );
        sendValues.swap( recvValues );
    }

    SCAI_ASSERT_EQ_ERROR( allIA.size(), numValues, "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( allJA.size(), numValues, "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( allValues.size(), numValues, "serious mismatch" )

    sparsekernel::COOUtils::normalize( allIA, allJA, allValues, op, Context::getHostPtr() );

    return COOStorage<ValueType>( numRows, numColumns,
                                  std::move( allIA ), std::move( allJA ), std::move( allValues ) );
}

/* ================================================================================ */
/*   Instantion of MatrixAssembly with numeric types                                */
/* ================================================================================ */

SCAI_COMMON_INST_CLASS( MatrixAssembly, SCAI_NUMERIC_TYPES_HOST )

}

}
