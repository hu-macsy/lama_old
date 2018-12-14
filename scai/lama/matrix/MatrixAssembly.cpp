/**
 * @file MatrixAssembly.cpp
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
 * @brief to a matrix to add matrix elements
 * @author Thomas Brandes
 * @date 07.09.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/matrix/MatrixAssembly.hpp>

#include <scai/dmemo/GlobalExchangePlan.hpp>
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

    HArray<PartitionId> owners;

    dist.computeOwners( owners, ia );
 
    // send the COO data to their owners and receive it

    dmemo::globalExchange( ownedIA, ownedJA, ownedValues, ia, ja, values, owners, dist.getCommunicator() );

    dist.global2LocalV( ownedIA, ownedIA );   // translates global row indexes to local ones

    sparsekernel::COOUtils::normalize( ownedIA, ownedJA, ownedValues, op, Context::getHostPtr() );

    return COOStorage<ValueType>( dist.getLocalSize(), numColumns,
                                  std::move( ownedIA ), std::move( ownedJA ), std::move( ownedValues ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType> MatrixAssembly<ValueType>::buildOwnedCOO(
    const dmemo::Distribution& dist,
    const IndexType numColumns,
    common::BinaryOp op ) const
{
    // These COO arrays will keep the matrix items owned by this processor

    HArray<IndexType> ownedIA;
    HArray<IndexType> ownedJA;
    HArray<ValueType> ownedValues;

    if ( dist.isReplicated() )
    {
        // all entries are used, we can save overhead for dist.isLocal and dist.global2Local

        IndexType numValues = mIA.size();

        WriteOnlyAccess<IndexType> wIA( ownedIA, numValues );
        WriteOnlyAccess<IndexType> wJA( ownedJA, numValues );
        WriteOnlyAccess<ValueType> wValues( ownedValues, numValues );

        for ( IndexType k = 0; k < numValues; ++k )
        {
            wIA[k] = mIA[k];
            wJA[k] = mJA[k];
            wValues[k] = mValues[k];
        }
    }
    else
    { 
        // Step 1 : count number of local entries

        IndexType numValues = 0;

        for ( size_t k = 0; k < mIA.size(); ++k )
        {
            if ( dist.isLocal( mIA[k] ) )
            {
                numValues++;
            }
        }
    
        // Step 2 : allocate the COO arrays and copy the local values

        WriteOnlyAccess<IndexType> wIA( ownedIA, numValues );
        WriteOnlyAccess<IndexType> wJA( ownedJA, numValues );
        WriteOnlyAccess<ValueType> wValues( ownedValues, numValues );

        IndexType count = 0;

        for ( size_t k = 0; k < mIA.size(); ++k )
        {
            const IndexType local = dist.global2Local( mIA[k] );

            if ( local != invalidIndex )
            {
                wIA[count] = local;       // IMPORTANT: global->local index required
                wJA[count] = mJA[k];
                wValues[count] = mValues[k];
                count++;
            }
        }

        SCAI_ASSERT_EQ_ERROR( count, numValues, "serious mismatch" )
    }

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

/* -------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType> MatrixAssembly<ValueType>::buildCOO( 
    const dmemo::Distribution& dist, 
    const IndexType numColumns,
    common::BinaryOp op ) const
{
    if ( mComm->getType() == dmemo::Communicator::NO )
    {
        // replicated assembly, so just select the owned entries

        SCAI_LOG_INFO( logger, "build (distributed) COO from replicated assembly, select owned entries" )

        return buildOwnedCOO( dist, numColumns, op );
    }
    else if ( dist.getCommunicator() == *mComm )
    {
        SCAI_LOG_INFO( logger, "build COO from assembly, same processor set" )

        // distributed assembly -> distributed matrix, sends non-local entries to owners

        return buildLocalCOO( dist, numColumns, op );
    }
    else if ( dist.getCommunicator().getType() == dmemo::Communicator::NO )
    {
        // distributed assembly -> replicated matrix 

        SCAI_LOG_INFO( logger, "build replicated COO from distributed assembly, circular shift it around" )

        return buildGlobalCOO( dist.getGlobalSize(), numColumns, op );
    }
    else
    {
        COMMON_THROWEXCEPTION( "unhandled, mismatch of assembly processor set and new processor set" )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixAssembly<ValueType>::truncate( 
    const IndexType numRows,
    const IndexType numColumns )
{
    IndexType offset = 0;

    for ( size_t k = 0; k < mIA.size(); ++k )
    {
        if ( mIA[k] >= numRows || mJA[k] >= numColumns )
        {
            continue;   // skip this element
        }

        mIA[offset] = mIA[k];
        mJA[offset] = mJA[k];
        mValues[offset] = mValues[k];

        ++offset;
    }

    mIA.resize( offset );
    mJA.resize( offset );
    mValues.resize( offset );
}

/* ================================================================================ */
/*   Instantion of MatrixAssembly with numeric types                                */
/* ================================================================================ */

SCAI_COMMON_INST_CLASS( MatrixAssembly, SCAI_NUMERIC_TYPES_HOST )

}

}
