/**
 * @file CommunicationPlan.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Definition of a class that holds for each partition a quantity value.
 * @author Thomas Brandes, Jiri Kraus
 * @date 10.03.2011
 * $Id$
 */
#ifndef LAMA_COMMUNICATION_PLAN_HPP_
#define LAMA_COMMUNICATION_PLAN_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Printable.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/LAMAArray.hpp>

// logging
#include <logging/logging.hpp>

#include <vector>

namespace lama
{

class Communicator;
// forward declaration

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
class LAMA_DLL_IMPORTEXPORT CommunicationPlan: public Printable
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

    CommunicationPlan();

    /** Clear an existing object for communication plan. */

    void clear();

    /** Construct a communication plan by quantity for each partition
     *
     *  @param quantities vector of non-negative values, size is number of partitions
     *  @param compressFlag if true zero entries are removed
     *
     *  \code
     *  std:vector<IndexType> quantities;
     *  quantities.push_back(0);
     *  quantities.push_back(3);
     *  quantities.push_back(4);
     *  CommunicationPlan plan( quantities );
     *  \endcode
     */
    CommunicationPlan( const LAMAArray<IndexType>& quantities, bool compressFlag = true );

    /** Construct a communication plan by array of owners; quantity
     *  is computed by counting values for each partition.
     */
    CommunicationPlan(
        const PartitionId noPartitions,
        const std::vector<PartitionId>& owners,
        bool compressFlag = true );

    /** Construct a communication plan by an existing one */

    CommunicationPlan( const CommunicationPlan& other, const IndexType quantities[] );

    /** @brief Construct a communication plan by an existing one where each entry is multiplied by same factor.
     *
     * @param[in] other   is the original communication plan
     * @param[in] n       is the multiplicator (n >= 1)
     *
     * The new communication plan can be used to send / receive whole arrays of size n instead
     * of single values.
     */
    CommunicationPlan( const CommunicationPlan& other, const IndexType n );

    /** Destructor. */

    virtual ~CommunicationPlan();

    /** Allocate a communication plan by quantity for each partition
     *
     *  @param quantities vector of non-negative values, size is number of partitions
     *  @param compressFlag if true zero-entries are removed.
     *
     *  \code
     *  CommunicationPlan plan;
     *  std:vector<IndexType> quantities;
     *  quantities.push_back( 0 );
     *  quantities.push_back( 3 );
     *  quantities.push_back( 4 );
     *  plan.allocate( quantities );
     *  \endcode
     */
    void allocate( const LAMAArray<IndexType>& quantities, bool compressFlag = true );

    void allocate( const PartitionId noPartitions, const std::vector<PartitionId>& owners, bool compressFlag = true );

    /** Allocate a communication plan as the transposed plan of another communication plan.
     *
     *  processor[p].plan->entry[q].quantity  = processor[q].entry[p].quantity
     */
    void allocateTranspose( const CommunicationPlan& plan, const Communicator& comm );

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

    void writeAt( std::ostream& stream ) const;

private:

    bool mAllocated; //!< true, if plan has been allocated

    bool mCompressed; //!< true, if no entry has quantity 0

    std::vector<Entry> mEntries; //!< vector of entries

    IndexType mQuantity; //!< sum of quantities for all entries

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
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

}

#endif // LAMA_COMMUNICATION_PLAN_HPP_
