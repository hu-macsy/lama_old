/**
 * @file Distribution.hpp
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
 * @brief Definition of abstract base class for a one-dimensional distribution.
 * @author Jiri Kraus
 * @date 22.02.2011
 * $Id$
 */
#ifndef LAMA_DISTRIBUTION_HPP_
#define LAMA_DISTRIBUTION_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>
#include <lama/Printable.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/Communicator.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

namespace lama
{

typedef boost::shared_ptr<const class Distribution> DistributionPtr;

/** Abstract base class for a one-dimensional distribution.
 *
 * A distribution specifies a mapping from a global range to the
 * partitions of a Communicator.
 *
 * Default and copy constructor are not available for this class (noncopyable).
 */
class LAMA_DLL_IMPORTEXPORT Distribution: public Printable, private NonCopyable
{
public:

    /** Constructor for a distribution.
     *
     * @param[in] globalSize is the number of elements to distribute
     * @param[in] communicator specifies the partitions onto which the elements
     *            are distributed
     */
    Distribution( const IndexType globalSize, const CommunicatorPtr communicator );

    /** Same as Distribution( globalSize, NoCommunicator() )
     *
     * @param[in] globalSize is the number of elements to distribute
     */
    Distribution( const IndexType globalSize );

    /** Destructor of distribution. */

    virtual ~Distribution();

    /** Getter routine for the communicator of the distribution. */

    const Communicator& getCommunicator() const;

    /** Getter routine for the communicator as shared pointer. */

    CommunicatorPtr getCommunicatorPtr() const;

    /** Query for the number of partitions onto which the distribution is done. */

    PartitionId getNumPartitions() const;

    /** Query whether the distribution is a replication. */

    bool isReplicated() const
    {
        return getNumPartitions() == 1;
    }

    /** Query if the given global index is local for the calling rank
     *                          (e.g. process for an MPI Communicator)
     *
     * @param[in]    index the global index to query for locality to the calling rank
     * @return       if the passed global index is local to the calling rank
     */
    virtual bool isLocal( const IndexType index ) const = 0;

    /** Getter for the global number of elements that are distributed. */

    inline IndexType getGlobalSize() const;

    /** This method returns the number of local elements on the calling
     *  processor.
     */
    virtual IndexType getLocalSize() const = 0;

    /** Abstract method that translates a local index back to a global index.
     *
     * @param[in]   localIndex is the local index, 0 <= localIndex < getLocalSize()
     * @return      globalIndex with 0 <= globalIndex < getGlobalSize
     *
     * This method must be implemented by all base classes. It should throw
     * an exception if the argument is not in the valid range.
     */
    virtual IndexType local2global( const IndexType localIndex ) const = 0;

    /** Abstract method that translates a global index into a local index.
     *
     * @param[in] globalIndex with 0 <= globalIndex < getGlobalSize
     * @return    localIndex with 0 <= localIndex < getLocalSize() if local, nIndex otherwise
     *
     * This method must be implemented by all base classes. It should throw
     * an exception if the argument is not in the valid range.
     */
    virtual IndexType global2local( const IndexType globalIndex ) const = 0;

    /** Compute ownership for required indexes.
     *
     * The default solution is to communicate required indexes around
     * all partitions and each partition marks indexes with its id
     * if it is local. If ownership can be computed without communication,
     * this routine might be implemented more efficiently.
     *
     * @param[in] requiredIndexes   TODO[doxy] Complete Description.
     * @param[in] owners            TODO[doxy] Complete Description.
     */
    virtual void computeOwners( const std::vector<IndexType>& requiredIndexes, std::vector<PartitionId>& owners ) const;

    /**
     * TODO[doxy] Complete Description.
     *
     * @param[in] other   TODO[doxy] Complete Description.
     */
    virtual bool isEqual( const Distribution& other ) const = 0;

    virtual void writeAt( std::ostream& stream ) const;

    /** Check for equality of two distributions.
     *
     *  As verification of same distribution can be rather expensive,
     *  the operator might return false.
     *
     *  Attention: The operator must be conservative and only return
     *             true if the distributions are really equal.
     */
    bool operator==( const Distribution& other ) const;

    /** Check for inequality. Even if two distributions are not equal
     *  it might be the case that the mapping of elements to partititons
     *  is the same.
     */
    bool operator!=( const Distribution& other ) const;

    /** Replication of distributed data, one entry for each element of the global range
     *
     * @tparam     T1           TODO[doxy] Complete Description.
     * @tparam     T2           TODO[doxy] Complete Description.
     * @param[out] allValues    array over the global range will contain all values
     * @param[in]  localValues  array over the local range has only values for each partition
     */
    template<typename T1,typename T2>
    void replicate( T1* allValues, const T2* localValues ) const;

    /** Replication of distributed data, one line of n entries for each element of the global range
     *
     * @tparam     T1           TODO[doxy] Complete Description.
     * @tparam     T2           TODO[doxy] Complete Description.
     * @param[out] allValues    array over the global range will contain all values
     * @param[in]  localValues  array over the local range has only values for each partition
     * @param[in]  n            is number of entries in each line of data
     */
    template<typename T1,typename T2>
    void replicateN( T1* allValues, const T2* localValues, const IndexType n ) const;

    /** Replication of distributed data, several entries for each element of the global range
     *
     * @tparam     T             TODO[doxy] Complete Description.
     * @param[out] allValues     array with all values from all partitions
     * @param[in]  localValues   elements available on this partition
     * @param[in]  allOffsets    contains unique offset for each element of the global range
     *
     * Note: the size of allValues must be globalSize + 1, allOffsets[globalSize] is total number of values
     *
     * The offset array is used like in CSR sparse matrix storage for offsets and number of values.
     */
    template<typename T>
    void replicateRagged( T* allValues, const T* localValues, const IndexType* allOffsets ) const;

    /**
     * Master process prints out the distribution vector to file named "name.part".
     * Every row contains a single number: the index of the process, where the row is local.
     *
     * @param[in] name   output file gets named "name.part".
     */
    virtual void printDistributionVector( std::string name ) const = 0;

protected:

    IndexType mGlobalSize;
    CommunicatorPtr mCommunicator;

private:

    Distribution(); // disable default constructor

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

IndexType Distribution::getGlobalSize() const
{
    return mGlobalSize;
}

}

#endif // LAMA_DISTRIBUTION_HPP_
