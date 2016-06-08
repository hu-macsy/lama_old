/**
 * @file CommunicationPlan.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
    : mAllocated( false ), mQuantity( 0 )
{
    SCAI_LOG_INFO( logger, "Communication plan constructed, not allocated" )
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan(
    const PartitionId noPartitions,
    const PartitionId owners[],
    const IndexType nOwners,
    const bool compressFlag )
    : mAllocated( false ), mQuantity( 0 )
{
    allocate( noPartitions, owners, nOwners, compressFlag );
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan(
    const IndexType quantities[],
    const PartitionId noPartitions,
    const bool compressFlag )
    : mAllocated( false ), mQuantity( 0 )
{
    allocate( quantities, noPartitions, compressFlag );
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan( const CommunicationPlan& other, const IndexType quantities[] )

    : mAllocated( false ), mQuantity( 0 )
{
    SCAI_LOG_INFO( logger, "extend plan for quantity of quantites: " << other )

    SCAI_ASSERT_ERROR( other.allocated(), "non allocated plan cannot be extended" )

    // copy the entries of the available plan, partitionIds do not change

    mEntries = other.mEntries;

    // This routine only works for running offsets, but should be the case

    IndexType offset = 0;

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        Entry& entry = mEntries[i];
        SCAI_ASSERT_ERROR( entry.offset == offset, "illegal plan to extend for quantities" )
        offset += entry.quantity;
    }

    offset = 0; // running for new offsets

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        Entry& entry = mEntries[i];

        IndexType newQuantity = 0;

        for ( IndexType k = 0; k < entry.quantity; k++ )
        {
            newQuantity += quantities[entry.offset + k];
        }

        entry.quantity = newQuantity;
        entry.offset = offset;
        offset += newQuantity;
    }

    // set other member variables

    mAllocated = true;
    mCompressed = other.compressed();
    mQuantity = offset;

    SCAI_LOG_INFO( logger, "extended quantity plan: " << *this )
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan( const CommunicationPlan& other, const IndexType n )

    : mAllocated( false ), mQuantity( 0 )
{
    SCAI_LOG_INFO( logger, "Construct multiply plan: " << other << ", factor = " << n )

    SCAI_ASSERT_ERROR( other.allocated(), "non allocated plan cannot be extended" )

    // copy the entries of the available plan, partitionIds do not change

    mEntries = other.mEntries;

    // This routine only works for running offsets, but should be the case

    IndexType offset = 0; // running for new offsets

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        Entry& entry = mEntries[i];
        SCAI_ASSERT_EQ_DEBUG( entry.offset * n, offset, "offset mismatch" )
        entry.quantity *= n; // each quantitiy is multiplied by n
        entry.offset *= n; // offsets are also multiplied by n
        offset += entry.quantity;
    }

    // set other member variables

    mAllocated = true;
    mCompressed = other.compressed();
    mQuantity = offset;

    SCAI_LOG_INFO( logger, "multiplied quantity plan: " << *this )
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

void CommunicationPlan::allocate( const IndexType quantities[], const PartitionId noPartitions, bool compressFlag )
{
    SCAI_LOG_INFO( logger, "allocate plan for " << noPartitions << " partitions from quantities" )

    mEntries.resize( noPartitions );

    mQuantity = 0; // counts total quantity

    for ( PartitionId i = 0; i < noPartitions; ++i )
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

void CommunicationPlan::allocate(
    const PartitionId noPartitions,
    const PartitionId owners[],
    const IndexType nOwners,
    bool compressFlag )
{
    mEntries.resize( noPartitions );

    SCAI_LOG_INFO( logger, "allocate plan for " << noPartitions << " partitions from owners" )

    for ( PartitionId p = 0; p < noPartitions; ++p )
    {
        mEntries[p].quantity = 0;
        mEntries[p].partitionId = p;
    }

    for ( IndexType i = 0; i < nOwners; ++i )
    {
        const PartitionId& p = owners[i];
        SCAI_ASSERT( p >= 0 && p < noPartitions, "Illegal owner value: " << p << " at Position " << i )
        ++mEntries[p].quantity;
        SCAI_LOG_TRACE( logger, " entry for p = " << p << ", total = " << mEntries[p].quantity )
    }

    mQuantity = 0; // counts total quantity

    for ( PartitionId p = 0; p < noPartitions; ++p )
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

void CommunicationPlan::allocateTranspose( const CommunicationPlan& plan, const Communicator& comm )
{
    SCAI_ASSERT( plan.allocated(), "plan to reverse not allocated" )

    const int size = comm.getSize();
    comm.synchronize();

    // plan might be compressed, so build values again

    std::vector<IndexType> sendSizes( size, 0 );

    for ( PartitionId i = 0; i < plan.size(); ++i )
    {
        sendSizes[plan[i].partitionId] = plan[i].quantity;
    }

    std::vector<IndexType> recvSizes( size, 0 );

    // send each processor the number of indexes I require
    // and receive the number of indexes that I have to provide

    comm.all2all( &recvSizes[0], &sendSizes[0] );

    // now we can allocate by quantities

    allocate( &recvSizes[0], recvSizes.size() );
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

} /* end namespace dmemo */

} /* end namespace scai */
