/**
 * @file CommunicationPlan.cpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Implementation of methods for class CommunicationPlan
 * @author Thomas Brandes
 * @date 10.03.2011
 * @since 1.0.0
 */

// hpp
#include <lama/CommunicationPlan.hpp>

// others
#include <lama/Communicator.hpp>

#include <lama/exception/LAMAAssert.hpp>

namespace lama
{

/* ------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( CommunicationPlan::logger, "CommunicationPlan" )

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan()
    : mAllocated( false ), mQuantity( 0 )
{
    LAMA_LOG_INFO( logger, "Communication plan constructed, not allocated" )
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan(
    const PartitionId noPartitions,
    const PartitionId owners[],
    const IndexType   nOwners,
    const bool compressFlag )
    : mAllocated( false ), mQuantity( 0 )
{
    allocate( noPartitions, owners, nOwners, compressFlag );
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan( const IndexType quantities[], const PartitionId noPartitions, const bool compressFlag )
    : mAllocated( false ), mQuantity( 0 )
{
    allocate( quantities, noPartitions, compressFlag );
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan( const CommunicationPlan& other, const IndexType quantities[] )

    : mAllocated( false ), mQuantity( 0 )
{
    LAMA_LOG_INFO( logger, "extend plan for quantity of quantites: " << other )

    LAMA_ASSERT_ERROR( other.allocated(), "non allocated plan cannot be extended" )

    // copy the entries of the available plan, partitionIds do not change

    mEntries = other.mEntries;

    // This routine only works for running offsets, but should be the case

    IndexType offset = 0;

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        Entry& entry = mEntries[i];
        LAMA_ASSERT_ERROR( entry.offset == offset, "illegal plan to extend for quantities" )
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

    LAMA_LOG_INFO( logger, "extended quantity plan: " << *this )
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan( const CommunicationPlan& other, const IndexType n )

    : mAllocated( false ), mQuantity( 0 )
{
    LAMA_LOG_INFO( logger, "Construct multiply plan: " << other << ", factor = " << n )

    LAMA_ASSERT_ERROR( other.allocated(), "non allocated plan cannot be extended" )

    // copy the entries of the available plan, partitionIds do not change

    mEntries = other.mEntries;

    // This routine only works for running offsets, but should be the case

    IndexType offset = 0; // running for new offsets

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        Entry& entry = mEntries[i];
        LAMA_ASSERT_EQUAL_DEBUG( entry.offset * n, offset )
        entry.quantity *= n; // each quantitiy is multiplied by n
        entry.offset *= n; // offsets are also multiplied by n
        offset += entry.quantity;
    }

    // set other member variables

    mAllocated = true;
    mCompressed = other.compressed();
    mQuantity = offset;

    LAMA_LOG_INFO( logger, "multiplied quantity plan: " << *this )
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::~CommunicationPlan()
{
    LAMA_LOG_DEBUG( logger, "~CommunicationPlan" )
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

    std::vector<Entry>().swap( mEntries );   // clear mEntries reallocating 
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::allocate( const IndexType quantities[], const PartitionId noPartitions, bool compressFlag )
{
    LAMA_LOG_INFO( logger, "allocate plan for " << noPartitions << " partitions from quantities" )

    mEntries.resize( noPartitions );

    mQuantity = 0; // counts total quantity

    for ( PartitionId i = 0; i < noPartitions; ++i )
    {
        Entry& entry = mEntries[i];

        entry.quantity    = quantities[i];
        entry.offset      = mQuantity;
        entry.partitionId = i;

        mQuantity += entry.quantity;

        LAMA_LOG_TRACE( logger, "Entries["<<i<<"].quantity = "<<mEntries[i].quantity )
        LAMA_LOG_TRACE( logger, " mQuantity = " << mQuantity )
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

    LAMA_LOG_INFO( logger, "allocate plan for " << noPartitions << " partitions from owners" )

    for ( PartitionId p = 0; p < noPartitions; ++p )
    {
        mEntries[p].quantity = 0;
        mEntries[p].partitionId = p;
    }

    for ( IndexType i = 0; i < nOwners; ++i )
    {
        const PartitionId& p = owners[i];
        LAMA_ASSERT( p >= 0 && p < noPartitions,
                     "Illegal owner value: " << p << " at Position " << i )
        ++mEntries[p].quantity;
        LAMA_LOG_TRACE( logger, " entry for p = " << p << ", total = " << mEntries[p].quantity )
    }

    mQuantity = 0; // counts total quantity

    for ( PartitionId p = 0; p < noPartitions; ++p )
    {
        mEntries[p].offset = mQuantity;
        mQuantity += mEntries[p].quantity;
    }

    // this assertion should be valid by the above algorithm

    LAMA_ASSERT_EQUAL_DEBUG( mQuantity, nOwners )

    mAllocated = true;

    if ( compressFlag )
    {
        compress();
    }
}

/* ------------------------------------------------------------------------- */

IndexType CommunicationPlan::maxQuantity() const
{
    LAMA_ASSERT( allocated(), "Plan not allocated." )

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

{   // remove all entries with zero quantity

    PartitionId count = 0;

    for ( PartitionId pid = 0; pid < size(); pid++ )
    {
        LAMA_LOG_TRACE( logger, "Entries["<<pid<<"].quantity = "<<mEntries[pid].quantity )
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

    LAMA_LOG_INFO( logger,
                   "CommunicationPlan compressed from " << size() << " to " << count << " entries, total quantity = " << mQuantity )

    mEntries.resize( count );

    mCompressed = true;
}

/* ----------------------------------------------------------------------- */

void CommunicationPlan::allocateTranspose( const CommunicationPlan& plan, const Communicator& comm )
{
    LAMA_ASSERT( plan.allocated(), "plan to reverse not allocated" )

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

    comm.all2all( recvSizes.data(), sendSizes.data() );

    // now we can allocate by quantities

    allocate( recvSizes.data(), recvSizes.size() );
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

}
