/**
 * @file CommunicationPlan.hpp
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

    /** @brief Default constructor creates an empty plan with allocated() returns false
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

    /** Clear an existing object for communication plan. */

    void clear();

    /** Clear an free allocated data. */

    void purge();

    /** Construct a communication plan by array of owners; quantity
     *  is computed by counting values for each partition.
     *
     *  @param[in] nPartitions  number of partitions that exist
     *  @param[in] owners        array of owners, each owner must be a value from 0 to nPartitions - 1
     *  @param[in] nOwners       number of entries in owners
     *  @param[in] compressFlag  will compress the computed plan
     *
     *  Each partition id that appears in owners will be counted to get the number of entries to send or receive.
     *
     *  \code
     *  PartitionId owners[] = { 2, 1, 2, 1, 2, 1, 2 };
     *  PartitionId nOwners = sizeof( owners ) / sizeof ( PartitionId );
     *  CommunicationPlan plan( 3, owners, nOwner );
     *  \endcode
     */
    CommunicationPlan(
        const PartitionId nPartitions,
        const PartitionId owners[],
        const IndexType nOwners,
        bool compressFlag = true );

    /** @brief Build a communication plan from this one where each entry is multiplied by same factor.
     *
     * @param[in] n       is the multiplicator (n >= 1)
     *
     * The new communication plan can be used to send / receive whole arrays of size n instead
     * of single values.
     */
    void multiplyConst( const IndexType n );

    /** Update this communication plan 
     *
     *  @param[in] quantities must have mQuantity entries 
     *
     *  ToDo: example
     */
    void multiplyRagged( const IndexType quantities[] );

    /** Allocate a communication plan by sizes for each partition
     *
     *  @param quantities array of non-negative values, size is number of partitions
     *  @param nPartitions number of entries to take from quantities
     *  @param compressFlag if true zero-entries are removed.
     *
     *  \code
     *  CommunicationPlan plan;
     *  IndexType quantities[] = { 0, 3, 4 };
     *  PartitionId nPartitions = sizeof( quantities ) / sizeof ( IndexType );
     *  plan.allocateBySizes( quantities, nPartitions );
     *  \endcode
     */
    void allocateBySizes( const IndexType quantities[], const PartitionId nPartitions, bool compressFlag = true );

    /** Allocate a communication plan by offsets instead of quantities for each partition.
     */
    void allocateByOffsets( const IndexType offsets[], const PartitionId nPartitions, bool compressFlag = true );

    /** @brief Allocate communication plan by an array of owners.
     *
     *  @param[in] nPartitions  number of partitions that exist
     *  @param[in] owners        array of owners, each owner must be a value from 0 to nPartitions - 1
     *  @param[in] nOwners       number of entries in owners
     *  @param[in] compressFlag  will compress the computed plan
     *
     *  Each partition id that appears in owners will be counted to get the number of entries to send or receive.
     *
     *  \code
     *  CommunicationPlan plan;
     *  PartitionId owners[] = { 2, 1, 2, 1, 2, 1, 2 };
     *  PartitionId nOwners = sizeof( owners ) / sizeof ( PartitionId );
     *  plan.allocateByOwners( 3, owners, nOwner );
     *  \endcode
     */
    void allocateByOwners( const PartitionId nPartitions, const PartitionId owners[], IndexType nOwners, bool compressFlag = true );

    /** Build a new communication plan that is the inverse of this plan 
     *
     *  processor[p].plan->entry[q].quantity  = processor[q].entry[p].quantity
     */
    CommunicationPlan transpose( const Communicator& comm ) const;

    /** @brief Query whether this communication plan has already been allocated */

    bool allocated() const;

    /** Compress a communication plan to remove all entries where quantity is zero. */

    void compress();

    /** Predicate to ask whether a communication plan is already compressed. */

    bool compressed() const;

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
    inline const Entry& operator[]( const PartitionId index ) const;

    /** @brief Overrides the routine of base class Printable. */

    void writeAt( std::ostream& stream ) const;

    /** @brief Query info for a certain partition.
     *
     *  @param[out] quantity is number of entries to communicate with p
     *  @param[out] offset is the corresponding offset
     *  @param[in] p is the partition for which info is needed
     */
    void getInfo( IndexType& quantity, IndexType& offset, PartitionId p ) const;

    /** @brief build a communnication plan from existing one for one processor only */

    void extractPlan( const CommunicationPlan& oldPlan, const PartitionId p );

    /** @brief build a communication plan with a single entry only */

    void singleEntry( const PartitionId p, const IndexType quantity );

    /** Build a communication plan by sizes for each partition
     *
     *  @param sizes array of non-negative values, size is nPartitions
     *  @param nPartitions number of entries for quantities
     *  @param compressFlag if true zero-entries are removed.
     *
     *  \code
     *  IndexType quantities[] = { 0, 3, 4 };
     *  PartitionId nPartitions = sizeof( quantities ) / sizeof ( IndexType );
     *  auto plan = CommunicationPlan::buildBySizes( quantities, nPartitions );
     *  \endcode
     */
    static CommunicationPlan buildBySizes( const IndexType sizes[], const PartitionId nPartitions, bool compressFlag = true );

    /** Build a new communication plan by offset array. */
 
    static CommunicationPlan buildByOffsets( const IndexType offsets[], const PartitionId nPartitions, bool compressFlag = true );

private:

    bool mAllocated; //!< true, if plan has been allocated

    bool mCompressed; //!< true, if no entry has quantity 0

    std::vector<Entry> mEntries; //!< vector of entries

    IndexType mQuantity; //!< sum of quantities for all entries

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* ----------------------------------------------------------------------*/

PartitionId CommunicationPlan::size() const
{
    return static_cast<PartitionId>( mEntries.size() );
}

const CommunicationPlan::Entry& CommunicationPlan::operator[]( const PartitionId id ) const
{
    return mEntries.at( id );
}

IndexType CommunicationPlan::totalQuantity() const
{
    return mQuantity;
}

} /* end namespace dmemo */

} /* end namespace scai */
