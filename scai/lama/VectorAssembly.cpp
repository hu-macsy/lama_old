/**
 * @file VectorAssembly.cpp
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
 * @brief Implementation of methods for vector assembly.
 * @author Thomas Brandes
 * @date 07.09.2017
 */

#include <scai/lama.hpp>

#include <scai/lama/VectorAssembly.hpp>

#include <scai/common/macros/instantiate.hpp>
#include <algorithm>

namespace scai
{

using namespace hmemo;
using utilskernel::HArrayUtils;

namespace lama
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, VectorAssembly<ValueType>::logger,
                              "VectorAssembly" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
VectorAssembly<ValueType>::VectorAssembly( dmemo::CommunicatorPtr comm )
{
    mComm = comm;
}

template<typename ValueType>
VectorAssembly<ValueType>::~VectorAssembly()
{
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssembly<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "VectorAssembly<" << common::TypeTraits<ValueType>::id()
           << ">( " << *mComm << ": " << mIA.size() << " entries )";
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssembly<ValueType>::reserve( const IndexType n )
{
    mIA.reserve( n );
    mValues.reserve( n );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType VectorAssembly<ValueType>::getSize() const
{
    auto it = std::max_element( std::begin(mIA), std::end(mIA) );

    IndexType numRows = it == end( mIA ) ? 0 : *it;

    return mComm->max( numRows ) + 1;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType VectorAssembly<ValueType>::getNumValues() const
{
    return mComm->sum( mIA.size() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssembly<ValueType>::exchangeCOO( 
    HArray<IndexType>& outIA,
    HArray<ValueType>& outValues,
    const HArray<IndexType> inIA,
    const HArray<ValueType> inValues,
    const dmemo::Distribution& dist ) const
{
    HArray<PartitionId> owners;

    dist.computeOwners( owners, inIA );

    SCAI_LOG_DEBUG( logger, "owners = " << owners )

    const dmemo::Communicator& comm = dist.getReduceCommunicator();
    PartitionId np = comm.getSize();

    HArray<IndexType> perm;
    HArray<IndexType> offsets;

    HArrayUtils::bucketSort( offsets, perm, owners, np );

    SCAI_LOG_DEBUG( logger, "sorted, perm = " << perm << ", offsets = " << offsets )

    HArray<IndexType> sendIA;
    HArray<ValueType> sendValues;

    HArrayUtils::gather( sendIA, inIA, perm, common::BinaryOp::COPY );
    HArrayUtils::gather( sendValues, inValues, perm, common::BinaryOp::COPY );

    auto sendPlan = dmemo::CommunicationPlan::buildByOffsets( hostReadAccess( offsets ).get(), np );
    auto recvPlan = sendPlan.transpose( comm );

    SCAI_LOG_DEBUG( logger, "recv plan: " << recvPlan )

    comm.exchangeByPlan( outIA, recvPlan, sendIA, sendPlan );
    comm.exchangeByPlan( outValues, recvPlan, sendValues, sendPlan );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssembly<ValueType>::checkLegalIndexes( const IndexType n ) const
{
    using namespace utilskernel;

    HArrayRef<IndexType> ia( mIA );

    bool okay = HArrayUtils::validIndexes( ia, n );

    if ( !okay )
    {
        SCAI_LOG_ERROR( logger, *mComm << ": illegal index in assembly, n = " << n )
    }

    okay = mComm->all( okay );

    if ( !okay )
    {
        COMMON_THROWEXCEPTION( "VectorAssembly has entries out of range " << n );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssembly<ValueType>::buildLocalData(
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<ValueType>& values,
    const dmemo::Distribution& dist ) const
{
    SCAI_ASSERT_EQ_ERROR( dist.getCommunicator(), *mComm, "VectorAssembly collected for comm = " << *mComm 
                          << ", cannot be used to distribute onto " << dist.getReduceCommunicator() )

    checkLegalIndexes( dist.getGlobalSize() );

    if ( dist.isReplicated() )
    {
        const HArrayRef<IndexType> localIA( mIA );
        const HArrayRef<ValueType> localValues( mValues );

        ia = localIA;
        values = localValues;
    }
    else
    {
        const HArrayRef<IndexType> localIA( mIA );
        const HArrayRef<ValueType> localValues( mValues );

        // These arrays will keep the matrix items owned by this processor

        ia.clear();
        values.clear();

        exchangeCOO( ia, values, localIA, localValues, dist );
        dist.global2localV( ia, ia );   // translates global indexes to local ones
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssembly<ValueType>::buildGlobalData(
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<ValueType>& values,
    const IndexType n ) const
{
    ContextPtr contextPtr = Context::getHostPtr();

    checkLegalIndexes( n );

    const IndexType numValues = mComm->sum( mIA.size() );
    const IndexType maxSize   = mComm->max( mIA.size() );

    // These arrays will keep all entries; reserve numValues entries to avoid reallocations

    ia.reserve( contextPtr, numValues );
    values.reserve( contextPtr, numValues );

    HArray<IndexType> sendIA;
    HArray<ValueType> sendValues;

    HArray<IndexType> recvIA;
    HArray<ValueType> recvValues;

    sendIA.reserve( contextPtr, maxSize );
    sendValues.reserve( contextPtr, maxSize );

    recvIA.reserve( contextPtr, maxSize );
    recvValues.reserve( contextPtr, maxSize );

    HArrayUtils::appendArray( ia, HArrayRef<IndexType>( mIA ) );
    HArrayUtils::appendArray( values, HArrayRef<ValueType>( mValues ) );

    // np - 1 shift steps are neeed

    const PartitionId np = mComm->getSize();

    const int COMM_DIRECTION = 1;  // circular shifting from left to right

    for ( PartitionId p = 0; p < np - 1; ++p )
    {
        if ( p == 0 )
        {
            mComm->shiftArray( recvIA, ia, COMM_DIRECTION );
            mComm->shiftArray( recvValues, values, COMM_DIRECTION );
        }
        else
        {
            mComm->shiftArray( recvIA, sendIA, COMM_DIRECTION );
            mComm->shiftArray( recvValues, sendValues, COMM_DIRECTION );
        }

        HArrayUtils::appendArray( ia, recvIA );
        HArrayUtils::appendArray( values, recvValues );

        // prepare for next step, the received values from left will be sent to right

        sendIA.swap( recvIA );
        sendValues.swap( recvValues );
    }

    SCAI_ASSERT_EQ_ERROR( ia.size(), numValues, "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( values.size(), numValues, "serious mismatch" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void VectorAssembly<ValueType>::truncate( const IndexType size )
{
    IndexType offset = 0;

    for ( size_t k = 0; k < mIA.size(); ++k )
    {
        if ( mIA[k] >= size )
        {
            continue;   // skip this element
        }

        mIA[offset] = mIA[k];
        mValues[offset] = mValues[k];

        ++offset;
    }

    SCAI_LOG_INFO( logger, "truncate size = " << size << ", now " << offset << " entries, were " << mIA.size() << " before" )

    mIA.resize( offset );
    mValues.resize( offset );
}

/* -------------------------------------------------------------------------- */


SCAI_COMMON_INST_CLASS( VectorAssembly, SCAI_ARRAY_TYPES_HOST )

}

}
