/**
 * @file CommunicationPlan.hpp
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
 * @brief Definition of a class that holds for each partition a quantity value.
 * @author Thomas Brandes, Jiri Kraus
 * @date 10.03.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

// std
#include <vector>

namespace scai
{

namespace dmemo
{

// forward declaration

class Communicator;

/**
 * A communication plan describes a schedule of data to send
 * or data to receive from other processors.
 *
 * A communication plan contains a certain number of entries where
 * one entry contains the partner (source or target) of the communication
 * and the number of elements to send or to receive.
 *
 * Note: the number of entries can be smaller than the number of available
 * partitions/processors. A compressed plan will not contain any entry
 * where the number of elements is zero.
 *
 */
class COMMON_DLL_IMPORTEXPORT CommunicationPlan: public scai::common::Printable
{
public:

    /**
     * @brief Record that is used for the entries of the communication plan.
     */
    struct Entry
    {
        PartitionId partitionId; //!< partner for communication
        IndexType quantity; //!< number of values to communicate
        IndexType offset; //!< running sum for quantities of all entries
    };

    /** @brief Default constructor creates an empty plan
     */
    CommunicationPlan();

    /** Default copy constructor can be used. */

    CommunicationPlan( const CommunicationPlan& ) = default;

    /** Default move constructor can be used. */

    CommunicationPlan( CommunicationPlan&& ) = default;

    /** Default copy assignment operator can be used. */

    CommunicationPlan& operator=( const CommunicationPlan& ) = default;

    /** Default move assignment operator can be used. */

    CommunicationPlan& operator=( CommunicationPlan&& ) = default;

    /** Destructor */

    ~CommunicationPlan();

    void swap( CommunicationPlan& other );

    /** Reset to a zero communication plan. */

    void clear();

    /** Clear an free allocated data. */

    void purge();

    /** Construct a communication plan by the number of entries either required or provided for each processor.
     *
     *  @param[in] quantities   array of sizes for each processor
     *  @param[in] size         number of processors
     *
     *  \code
     *  IndexType quantities[] = { 3, 1, 0, 1 };
     *  PartitionId size = 4;
     *  CommunicationPlan plan( sizes, 4 );
     *  \endcode
     */
    CommunicationPlan( const IndexType quantities[], const PartitionId size );

    /** Construct a communication plan by iterating over an object with quantities 
     *
     *  @param object must have an iterable range with for quantities
     */
    template<typename T>
    inline explicit CommunicationPlan( const T& object );

    /** Define this communication plan by iterating over an object with quantities 
     *
     *  @param object must have an iterable range with for quantities
     */
    template<typename T>
    void defineByQuantities( const T& object );

    /** @brief Build a communication plan from this one where each entry is multiplied by same factor.
     *
     * @param[in] n       is the multiplicator (n >= 1)
     *
     * The new communication plan can be used to send / receive whole arrays of size n instead
     * of single values.
     */
    void multiplyConst( const IndexType n );

    /** 
     *  @brief Build in-place a communication plan where where each communicated element is replaced
     *         with an array whose size is individually for each element (ragged array, i.e. array of arrays)
     * 
     *  @param[in] sizes specifies for each entry to be communicated the size of its array
     */
    void multiplyRaggedBySizes( const IndexType sizes[] );

    /**
     *  Update this communication plan with variant sizes for each entry, using the offset array.
     */
    void multiplyRaggedByOffsets( const IndexType offsets[] );

    /** Define a communication plan by sizes for each partition
     *
     *  @param quantities array of non-negative values, size is number of partitions
     *  @param nPartitions number of entries to take from quantities
     *
     *  \code
     *  CommunicationPlan plan;
     *  IndexType quantities[] = { 0, 3, 4 };
     *  PartitionId nPartitions = sizeof( quantities ) / sizeof ( IndexType );
     *  plan.defineByQuantities( quantities, nPartitions );
     *  \endcode
     */
    void defineByQuantities( const IndexType quantities[], const PartitionId nPartitions );

    /** Define a communication plan by one quantity for one other processor.
     *
     *  @param[in] quantity  number of entries for exchange with processor rank
     *  @param[in] partner   rank of the processor with which data is exchanged.
     *
     *  \code
     *    PartitionId partner = 1;
     *    IndexType N = 125;
     *    CommunicationPlan plan;
     *    plan.defineBySingleEntry( N, partner );
     *  \endcode
     */
    void defineBySingleEntry( const IndexType quantity, const PartitionId partner );

    /** Query the number of entries.
     *
     * @returns number of entries stored in the communication plan.
     *
     * Attention: keep in mind that entries with a zero quantity are deleted
     * and that this method returns a value that might be smaller than the original
     * number of partitions.
     */
    inline PartitionId size() const;

    /** Get the sum of all quantities, is total amount of data to send or receive. */

    inline IndexType totalQuantity() const;

    /**
     * Get the maximal quantity of any partition.
     *
     * This method is useful for the allocation of send or receive buffers that
     * will be used for all communications.
     */

    IndexType maxQuantity() const;

    /** Index operator
     *
     * @param[in] index   is an index value, 0 <= id < size()
     *
     * Attention: keep in mind that size() might be smaller than
     * the number of partitions involved.
     */
    inline const Entry& operator[]( PartitionId index ) const;

    /** @brief Overrides the routine of base class Printable. */

    void writeAt( std::ostream& stream ) const;

    /** @brief Query info for a certain partition.
     *
     *  @param[out] quantity is number of entries to communicate with p
     *  @param[out] offset is the corresponding offset
     *  @param[in] p is the partition for which info is needed
     */
    void getInfo( IndexType& quantity, IndexType& offset, const PartitionId p ) const;

    /** @brief Remove an entry in the communication plan and get its values
     *
     *  @param[in] p is the processor for which entry is removed
     *  @param[out] quantity is number of entries to communicate with p
     *  @param[out] offset is the corresponding offset
     */
    void removeEntry( IndexType& quantity, IndexType& offset, const PartitionId p );

    /** @brief build a communnication plan from existing one for one processor only */

    void extractPlan( const CommunicationPlan& oldPlan, const PartitionId p );

    CommunicationPlan constructRaggedBySizes( const hmemo::HArray<IndexType>& sizes ) const;

    CommunicationPlan constructRaggedByOffsets( const hmemo::HArray<IndexType>& offsets ) const;

    CommunicationPlan constructN( const IndexType N ) const;

    /** Build a communication plan by sizes for each partition
     *
     *  @param quantities array of non-negative values, has size entries
     *  @param size number of entries for quantities
     *
     *  \code
     *  IndexType quantities[] = { 0, 3, 4 };
     *  PartitionId np = sizeof( quantities ) / sizeof ( IndexType );
     *  auto plan = CommunicationPlan::buildByQuantities( quantities, np );
     *  \endcode
     */
    static CommunicationPlan buildByQuantities( const IndexType quantities[], const PartitionId size );

    /** Build a communication plan with a single entry.
     *
     *  @param quantity  number of entries to communicate with partner p 
     *  @param p         is the only processor for which communication is done
     */
    static CommunicationPlan buildSingle( const IndexType quantity, const PartitionId p );

private:

    std::vector<Entry> mEntries; //!< vector of entries

    IndexType mQuantity; //!< sum of quantities for all entries

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* ----------------------------------------------------------------------*/

template<typename T>
CommunicationPlan::CommunicationPlan( const T& object )
{
    defineByQuantities( object );
}

template<>
inline CommunicationPlan::CommunicationPlan( const hmemo::HArray<IndexType>& object )
{
    defineByQuantities( hostReadAccess( object ) );
}

template<typename T>
void CommunicationPlan::defineByQuantities( const T& object )
{
    IndexType countEntries = 0;

    for ( typename T::const_iterator iter = object.begin(); iter != object.end(); ++iter )
    {
        if ( *iter > 0 )
        {
            countEntries++;
        }
    }

    mEntries.resize( countEntries );

    mQuantity = 0; // counts total quantity

    countEntries = 0;  // reset 

    PartitionId p = 0;

    for ( typename T::const_iterator iter = object.begin(); iter != object.end(); iter++ )
    {
        if ( *iter > 0 )
        {
            Entry& entry = mEntries[countEntries];
            entry.quantity = static_cast<IndexType>( *iter );
            entry.offset = mQuantity;
            entry.partitionId = p;
            mQuantity += entry.quantity;
            countEntries++;
        }
        ++p;
    }
}

PartitionId CommunicationPlan::size() const
{
    return static_cast<PartitionId>( mEntries.size() );
}

const CommunicationPlan::Entry& CommunicationPlan::operator[]( PartitionId id ) const
{
    return mEntries.at( id );
}

IndexType CommunicationPlan::totalQuantity() const
{
    return mQuantity;
}

} /* end namespace dmemo */

} /* end namespace scai */
