/**
 * @file CommunicationPlan.cpp
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
 * @brief Implementation of methods for class CommunicationPlan
 * @author Thomas Brandes
 * @date 10.03.2011
 */

// hpp
#include <scai/dmemo/CommunicationPlan.hpp>

// local library
#include <scai/dmemo/Communicator.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace dmemo
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( CommunicationPlan::logger, "CommunicationPlan" )

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan()
    : mAllocated( false ), mCompressed( true ), mQuantity( 0 )
{
    SCAI_LOG_INFO( logger, "Communication plan constructed, not allocated" )
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan(
    const PartitionId nPartitions,
    const PartitionId owners[],
    const IndexType nOwners,
    const bool compressFlag )
    : mAllocated( false ), mQuantity( 0 )
{
    allocateByOwners( nPartitions, owners, nOwners, compressFlag );
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::multiplyConst( const IndexType n )
{
    for ( size_t i = 0; i < mEntries.size(); ++i )
    {
        Entry& entry = mEntries[i];
        entry.quantity *= n;          // each quantitiy is multiplied by n
        entry.offset *= n;            // offsets are also multiplied by n
    }

    mQuantity *= n;                   // update of total quantity
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::multiplyRagged( const IndexType quantities[] )
{
    SCAI_LOG_INFO( logger, "extend plan with individual multiplicators: " << *this )

    SCAI_ASSERT_ERROR( allocated(), "non allocated plan cannot be extended" )

    IndexType newOffset = 0;    // new running sum 
    IndexType oldOffset = 0;    // old running sum 

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        Entry& entry = mEntries[i];

        SCAI_ASSERT_EQ_ERROR( entry.offset, oldOffset, "serious mismatch" )

        IndexType newQuantity = 0;

        for ( IndexType k = 0; k < entry.quantity; k++ )
        {
            newQuantity += quantities[entry.offset + k];
        }

        oldOffset += entry.quantity;

        entry.quantity = newQuantity;
        entry.offset = newOffset;
        newOffset += newQuantity;
    }

    mQuantity = newOffset;

    SCAI_LOG_INFO( logger, "extended quantity plan: " << *this )
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::~CommunicationPlan()
{
    SCAI_LOG_DEBUG( logger, "~CommunicationPlan" )
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::clear()
{
    mAllocated = false;
    mQuantity = 0;
    mCompressed = true; // does not matter, but no entry has quantity 0
    mEntries.clear();
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::purge()
{
    mAllocated = false;
    mQuantity = 0;
    mCompressed = true; // does not matter, but no entry has quantity 0
    std::vector<Entry>().swap( mEntries ); // clear mEntries reallocating
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::allocateBySizes( const IndexType quantities[], const PartitionId nPartitions, bool compressFlag )
{
    mCompressed = false;

    SCAI_LOG_INFO( logger, "allocate plan for " << nPartitions << " partitions from quantities" )
    mEntries.resize( nPartitions );
    mQuantity = 0; // counts total quantity

    for ( PartitionId i = 0; i < nPartitions; ++i )
    {
        Entry& entry = mEntries[i];
        entry.quantity = quantities[i];
        entry.offset = mQuantity;
        entry.partitionId = i;
        mQuantity += entry.quantity;
        SCAI_LOG_TRACE( logger, "Entries[" << i << "].quantity = " << mEntries[i].quantity )
        SCAI_LOG_TRACE( logger, " mQuantity = " << mQuantity )
    }

    mAllocated = true;

    if ( compressFlag )
    {
        compress();
    }
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::allocateByOffsets( const IndexType offsets[], const PartitionId nPartitions, bool compressFlag )
{
    mCompressed = false;

    SCAI_LOG_INFO( logger, "allocate plan for " << nPartitions << " partitions by offsets" )

    mEntries.resize( nPartitions );

    SCAI_ASSERT_EQ_ERROR( offsets[0], 0, "Illegal offsets array" )

    for ( PartitionId i = 0; i < nPartitions; ++i )
    {
        Entry& entry = mEntries[i];
        SCAI_ASSERT_LE_ERROR( offsets[i], offsets[i + 1], "illegal offsets for i = " << i )
        entry.quantity = offsets[i + 1] - offsets[i];
        entry.offset = offsets[i];
        entry.partitionId = i;
        SCAI_LOG_TRACE( logger, "Entries[" << i << "].quantity = " << mEntries[i].quantity )
        SCAI_LOG_TRACE( logger, " mQuantity = " << mQuantity )
    }

    mQuantity = offsets[nPartitions];   // total quantity

    mAllocated = true;

    if ( compressFlag )
    {
        compress();
    }
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::allocateByOwners(
    const PartitionId nPartitions,
    const PartitionId owners[],
    const IndexType nOwners,
    bool compressFlag )
{
    mCompressed = false;

    mEntries.resize( nPartitions );
    SCAI_LOG_INFO( logger, "allocate plan for " << nPartitions << " partitions from owners" )

    for ( PartitionId p = 0; p < nPartitions; ++p )
    {
        mEntries[p].quantity = 0;
        mEntries[p].partitionId = p;
    }

    for ( IndexType i = 0; i < nOwners; ++i )
    {
        const PartitionId& p = owners[i];
        SCAI_ASSERT_VALID_INDEX( p, nPartitions, "Illegal owner value at owners[ " << i << "]" )
        ++mEntries[p].quantity;
        SCAI_LOG_TRACE( logger, " entry for p = " << p << ", total = " << mEntries[p].quantity )
    }

    mQuantity = 0; // counts total quantity

    for ( PartitionId p = 0; p < nPartitions; ++p )
    {
        mEntries[p].offset = mQuantity;
        mQuantity += mEntries[p].quantity;
    }

    // this assertion should be valid by the above algorithm
    SCAI_ASSERT_EQ_DEBUG( mQuantity, nOwners, "wrong sum up" )
    mAllocated = true;

    if ( compressFlag )
    {
        compress();
    }
}

/* ------------------------------------------------------------------------- */

IndexType CommunicationPlan::maxQuantity() const
{
    SCAI_ASSERT( allocated(), "Plan not allocated." )
    // Get the maximal quantity of another proc
    IndexType quantity = 0;

    for ( PartitionId id = 0; id < size(); id++ )
    {
        quantity = std::max( quantity, mEntries[id].quantity );
    }

    return quantity;
}

/* ------------------------------------------------------------------------- */

bool CommunicationPlan::allocated() const
{
    return mAllocated;
}

/* ------------------------------------------------------------------------- */

bool CommunicationPlan::compressed() const
{
    return mCompressed;
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::compress()

{
    // remove all entries with zero quantity
    PartitionId count = 0;

    for ( PartitionId pid = 0; pid < size(); pid++ )
    {
        SCAI_LOG_TRACE( logger, "Entries[" << pid << "].quantity = " << mEntries[pid].quantity )

        if ( mEntries[pid].quantity == 0 )
        {
            continue;
        }

        if ( count != pid )
        {
            mEntries[count] = mEntries[pid];
        }

        count++;
    }

    SCAI_LOG_INFO( logger,
                   "CommunicationPlan compressed from " << size() << " to " << count << " entries, total quantity = " << mQuantity )
    mEntries.resize( count );
    mCompressed = true;
}

/* ----------------------------------------------------------------------- */

void CommunicationPlan::getInfo( IndexType& quantity, IndexType& offset, PartitionId p ) const
{
    // set initial values in case we do not find an entry for p

    quantity = 0;
    offset   = 0;

    for ( PartitionId pid = 0; pid < size(); pid++ )
    {
        if ( mEntries[pid].partitionId == p )
        {
            quantity = mEntries[pid].quantity;
            offset = mEntries[pid].offset;
            break;
        }
    }
}

/* ----------------------------------------------------------------------- */

void CommunicationPlan::extractPlan( const CommunicationPlan& oldPlan, const PartitionId p )
{
    mEntries.clear();
    mQuantity = 0;

    for ( PartitionId pid = 0; pid < oldPlan.size(); pid++ )
    {
        const Entry& entry = oldPlan.mEntries[pid];

        if ( entry.partitionId == p  && entry.quantity > 0 )
        {
            Entry newEntry( entry );

            newEntry.offset = 0;

            mEntries.push_back( newEntry );

            mQuantity += entry.quantity;
        }
    }

    mAllocated = true;
    mCompressed = true;
}

/* ----------------------------------------------------------------------- */

void CommunicationPlan::singleEntry( const PartitionId p, const IndexType quantity )
{
    mEntries.clear();
    mQuantity = 0;

    if ( quantity > 0 )
    {
        Entry entry;

        entry.partitionId = p;
        entry.quantity    = quantity;
        entry.offset      = 0;

        mQuantity += quantity;

        mEntries.push_back( entry );
    }

    mAllocated = true;
    mCompressed = true;
}

/* ----------------------------------------------------------------------- */

void CommunicationPlan::writeAt( std::ostream& stream ) const
{
    // stream output of a communication plan
    stream << "CommunicationPlan(size=" << mEntries.size() << ",quantity=" << mQuantity;

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        stream << ",->" << mEntries[i].partitionId << ":" << mEntries[i].quantity;
    }

    stream << ")";
}

/* ----------------------------------------------------------------------- */

CommunicationPlan CommunicationPlan::buildBySizes( const IndexType sizes[], const PartitionId nPartitions, bool compressFlag )
{
    CommunicationPlan plan;
    plan.allocateBySizes( sizes, nPartitions, compressFlag );
    return plan;
}

CommunicationPlan CommunicationPlan::buildByOffsets( const IndexType offsets[], const PartitionId nPartitions, bool compressFlag )
{
    CommunicationPlan plan;
    plan.allocateByOffsets( offsets, nPartitions, compressFlag );
    return plan;
}

CommunicationPlan CommunicationPlan::transpose( const Communicator& comm ) const
{
    SCAI_ASSERT_ERROR( allocated(), "plan not allocated: " << *this )

    const PartitionId np = comm.getSize();

    // plan might be compressed, so build values again

    std::vector<IndexType> sendSizes( np, 0 );

    for ( PartitionId i = 0; i < size(); ++i )
    {
        const Entry& entry = mEntries[i];
        sendSizes[entry.partitionId] = entry.quantity;
    }

    std::vector<IndexType> recvSizes( np, 0 );

    // send each processor the number of indexes I require
    // and receive the number of indexes that I have to provide

    comm.all2all( &recvSizes[0], &sendSizes[0] );

    // now we can allocate by quantities

    return buildBySizes( &recvSizes[0], np );
}


} /* end namespace dmemo */

} /* end namespace scai */
