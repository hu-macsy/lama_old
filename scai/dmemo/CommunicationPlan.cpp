/**
 * @file CommunicationPlan.cpp
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
 * @brief Implementation of methods for class CommunicationPlan
 * @author Thomas Brandes
 * @date 10.03.2011
 */

// hpp
#include <scai/dmemo/CommunicationPlan.hpp>

// local library
#include <scai/dmemo/Communicator.hpp>
#include <scai/hmemo/HArray.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace dmemo
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( CommunicationPlan::logger, "CommunicationPlan" )

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan() :

    mQuantity( 0 )

{
    SCAI_LOG_INFO( logger, "Communication plan constructed, not allocated" )
}

/* ------------------------------------------------------------------------- */

CommunicationPlan::CommunicationPlan( const IndexType quantities[], const PartitionId size )
{
    defineByQuantities( quantities, size );
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::defineBySingleEntry( const IndexType quantity, const PartitionId partner )
{
    mQuantity = 0;

    if ( quantity > 0 )
    {
        mEntries.resize( 1 );

        Entry& entry = mEntries[0];

        entry.partitionId = partner;
        entry.quantity    = quantity;
        entry.offset      = 0;

        mQuantity += quantity;
    }
    else
    {
        mEntries.clear();
    }
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::swap( CommunicationPlan& other )
{
    std::swap( mEntries, other.mEntries );
    std::swap( mQuantity, other.mQuantity );
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

void CommunicationPlan::multiplyRaggedBySizes( const IndexType sizes[] )
{
    SCAI_LOG_INFO( logger, "extend plan with individual multiplicators: " << *this )

    IndexType newOffset = 0;    // new running sum 
    IndexType oldOffset = 0;    // old running sum 

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        Entry& entry = mEntries[i];

        SCAI_ASSERT_EQ_DEBUG( entry.offset, oldOffset, "serious mismatch, illegal communication plan" )

        IndexType newQuantity = 0;

        for ( IndexType k = 0; k < entry.quantity; k++ )
        {
            newQuantity += sizes[entry.offset + k];
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

void CommunicationPlan::multiplyRaggedByOffsets( const IndexType offsets[] )
{
    SCAI_LOG_INFO( logger, "extend plan with individual multiplicators: " << *this )

    IndexType newOffset = 0;    // new running sum 
    IndexType oldOffset = 0;    // old running sum 

    for ( size_t i = 0; i < mEntries.size(); i++ )
    {
        Entry& entry = mEntries[i];

        SCAI_ASSERT_EQ_DEBUG( entry.offset, oldOffset, "serious mismatch, illegal communication plan" )

        IndexType newQuantity = offsets[entry.offset + entry.quantity] - offsets[entry.offset];

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
    mEntries.clear();
    mQuantity = 0;
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::purge()
{
    std::vector<Entry>().swap( mEntries ); // clear mEntries reallocating
    mQuantity = 0;
}

/* ------------------------------------------------------------------------- */

void CommunicationPlan::defineByQuantities( const IndexType quantities[], const PartitionId size )
{
    SCAI_LOG_INFO( logger, "define communication plan for " << size << " processors from array with quantities" )

    IndexType countEntries = 0;

    for ( PartitionId i = 0; i < size; ++i )
    {
        if ( quantities[i] > 0 )
        {
            countEntries++;
        }
    }

    mEntries.resize( countEntries );

    mQuantity = 0; // counts total quantity

    countEntries = 0;  // reset 

    for ( PartitionId i = 0; i < size; ++i )
    {
        if ( quantities[i] > 0 )
        {
            Entry& entry = mEntries[countEntries];
            entry.quantity = quantities[i];
            entry.offset = mQuantity;
            entry.partitionId = i;
            mQuantity += entry.quantity;
            countEntries++;
        }  
    }
}

/* ------------------------------------------------------------------------- */

IndexType CommunicationPlan::maxQuantity() const
{
    // Get the maximal quantity of another proc

    IndexType quantity = 0;

    for ( PartitionId id = 0; id < size(); id++ )
    {
        quantity = std::max( quantity, mEntries[id].quantity );
    }

    return quantity;
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

CommunicationPlan CommunicationPlan::buildByQuantities( const IndexType quantities[], const PartitionId size )
{
    CommunicationPlan plan;
    plan.defineByQuantities( quantities, size );
    return plan;
}

/* ----------------------------------------------------------------------- */

CommunicationPlan CommunicationPlan::constructRaggedBySizes( const hmemo::HArray<IndexType>& sizes ) const
{
    SCAI_ASSERT_EQ_ERROR( totalQuantity(), sizes.size(), "serious mismatch" )

    CommunicationPlan planV( *this );
    planV.multiplyRaggedBySizes( hmemo::hostReadAccess( sizes ).get() );
    return planV;
}

CommunicationPlan CommunicationPlan::constructRaggedByOffsets( const hmemo::HArray<IndexType>& offsets ) const
{
    SCAI_ASSERT_EQ_ERROR( totalQuantity() + 1, offsets.size(), "serious mismatch" )

    CommunicationPlan planV( *this );
    planV.multiplyRaggedByOffsets( hmemo::hostReadAccess( offsets ).get() );
    return planV;
}

CommunicationPlan CommunicationPlan::constructN( const IndexType N ) const
{
    CommunicationPlan planN( *this );
    planN.multiplyConst( N );
    return planN;
}

} /* end namespace dmemo */

} /* end namespace scai */
