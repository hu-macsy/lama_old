/**
 * @file Distribution.hpp
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
 * @brief Definition of abstract base class for a one-dimensional distribution.
 * @author Thomas Brandes, Jiri Kraus
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>
#include <scai/common/Factory1.hpp>

// local library
#include <scai/dmemo/Communicator.hpp>
#include <scai/hmemo/HArray.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

// std
#include <map>
#include <utility>
#include <memory>

namespace scai
{

namespace dmemo
{

typedef std::shared_ptr<const class Distribution> DistributionPtr;

class Distributed;

/** Structure that keeps all kind of arguments used to create a distribution.
 *
 *  An element of this struct is used as an argument of the create method of the factory.
 */

struct DistributionArguments
{
    DistributionArguments( CommunicatorPtr comm, IndexType size, const Distributed* m, float w )
    {
        communicator = comm;
        globalSize   = size;
        matrix       = m;
        weight       = w;
    }

    CommunicatorPtr communicator;
    IndexType globalSize;
    const Distributed* matrix;
    float weight;
};

/** Abstract base class for a one-dimensional distribution.
 *
 * A distribution specifies a mapping from a global range to the
 * partitions of a Communicator.
 *
 * Default and copy constructor are not available for this class (noncopyable).
 */
class COMMON_DLL_IMPORTEXPORT Distribution:

    public common::Factory1<std::string, DistributionArguments, Distribution*>,
    public common::Printable
{

public:

    /**
     * @brief Distribution factory to get a distribution of a certain kind and a certain type
     *
     * @param[in] kind specfies the name of the distribution, e.g. BLOCK, CYCLIC, GEN_BLOCK, GENERAL, METIS
     * @param[in] comm is the communicator used for the distribution
     * @param[in] globalSize is the number of elements to distribute
     * @param[in] weight is an individual weight for each partition of the communicator
     *
     * @returns pointer to a new distribution of the specified kind, NULL if kind is not supported
     *
     *  /code
     *  // Using a MetisDistribution requires its availabilty
     *  Distribution* dist = MetisDistribution( comm, size, weight )
     *  // code using the factory does not require the availability
     *  Distribution* dist = Distribution::getDistributionPtr( "METIS", comm, size, weight )
     *  if ( dist == NULL )
     *  {
     *      dist = Distribution::getDistributionPtr( "GEN_BLOCK", comm, size, weight )
     *  }
     *  /endcode
     *
     *  Note: Internally, this routine requires that all derived classes implement a corresponding
     *        create method that will be registered during static initialization.
     */
    static Distribution* getDistributionPtr(
        const std::string& kind,
        const CommunicatorPtr comm,
        const IndexType globalSize,
        const float weight = 1.0 );

    /**
     * @brief Distribution factory to get a distribution of a certain kind and a certain type
     *
     * @param[in] kind specfies the name of the distribution, e.g. BLOCK, CYCLIC, GEN_BLOCK, GENERAL, METIS
     * @param[in] comm is the communicator used for the distribution
     * @param[in] matrix is the matrix for which a good row distribution is determined
     * @param[in] weight is an individual weight for each partition of the communicator
     *
     * @returns pointer to a new distribution of the specified kind, NULL if kind is not supported
     *
     * Note: the current distribution of matrix does not matter.
     */
    static Distribution* getDistributionPtr(
        const std::string& kind,
        const CommunicatorPtr comm,
        const Distributed& matrix,
        const float weight = 1.0 );

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

    /** Distributions should never be copied. */

    Distribution( const Distribution& other ) = delete;

    /** Destructor of distribution. */

    virtual ~Distribution();

    /** Distributions should never be assigned */

    Distribution& operator=( const Distribution& other ) = delete;

    /** Getter routine for the communicator of the distribution. */

    const Communicator& getCommunicator() const;

    /** Getter routine for the communicator as shared pointer. */

    CommunicatorPtr getCommunicatorPtr() const;

    /** Query for the number of processors/partitions onto which the distribution is done. */

    PartitionId getNumPartitions() const;

    /** Query whether the distribution is a replication, i.e. each processor is owner of all data
     *
     *  same as getNumPartitions() == 1, always true for NoDistribution, always true when running
     *  an application on a single processor.
     */
    inline bool isReplicated() const;

    /** Each derived class has to return a specific kind string to specify its kind. */

    virtual const char* getKind() const = 0;

    /** Query if the given global index is local for the calling rank
     *                          (e.g. process for an MPI Communicator)
     *
     * @param[in]    globalIndex the global index to query for locality to the calling rank
     * @return       true if this processor/partition  is owner
     */
    virtual bool isLocal( const IndexType globalIndex ) const = 0;

    /** Getter for the global number of elements that are distributed. */

    inline IndexType getGlobalSize() const;

    /** This method returns the number of local elements on the calling
     *  processor.
     */
    virtual IndexType getLocalSize() const = 0;

    /**
     * @brief This method returns the maximal number of local elements on any processor.
     *
     * This method can be used to allocate data with the right size for any communication
     * that uses circular shifting for partition data.
     *
     * Default implementation: getCommunicator().max( getLocalSize() )
     */
    virtual IndexType getMaxLocalSize() const;

    /** Abstract method that translates a local index back to a global index.
     *
     * @param[in]   localIndex is the local index, 0 <= localIndex < getLocalSize()
     * @return      globalIndex with 0 <= globalIndex < getGlobalSize
     *
     * This method must be implemented by all base classes. It should throw
     * an exception if the argument is not in the valid range.
     */
    virtual IndexType local2global( const IndexType localIndex ) const = 0;

    /** This method translates a global index into a local index.
     *
     * @param[in] globalIndex with 0 <= globalIndex < getGlobalSize
     * @return    localIndex with 0 <= localIndex < getLocalSize() if local, invalidIndex otherwise
     *
     * This method must be implemented by all base classes. It should throw
     * an exception if the argument is not in the valid range.
     */
    virtual IndexType global2local( const IndexType globalIndex ) const = 0;

    /** This method translates a whole array of global indexes to local indexes.
     *
     *  @param[in]  globalIndexes is the array with the indexes that are translated
     *  @param[out] localIndexes is the array with local indexes
     *
     *  \code
     *  for ( IndexType i = 0; i < globalIndexes.size(); ++i )
     *  {
     *      localIndexes[i] = global2local( globalIndexes[i] );
     *  }
     *  \endcode
     *
     *  Note: alias of localIndexes and globalIndexes is supported
     */
    virtual void global2localV( hmemo::HArray<IndexType>& localIndexes, const hmemo::HArray<IndexType>& globalIndexes ) const;

    /** Get the owners for a set of (global) indexes
     *
     * The default solution is to communicate required indexes around
     * all partitions and each partition marks indexes with its id
     * if it is local. If ownership can be computed without communication,
     * this routine might be implemented more efficiently.
     *
     * @param[in] indexes is an array with global indexes, 0 <= indexes[i] < getGlobalSize()
     * @param[out] owners are the corresponing processors that own the indexes
     */
    virtual void computeOwners( hmemo::HArray<PartitionId>& owners, const hmemo::HArray<IndexType>& indexes ) const;

    /**
     * Get the owner of a global index, all processors call with same value.
     *
     *  @param[in] globalIndex index for which owner is required
     *  @return id of the owner partition.
     */
    virtual PartitionId findOwner( const IndexType globalIndex ) const;

    /** Get the owners of all global indexes.
     *
     *  @param[out] owners owners[i] is owner of element i, 0 <= i < globalSize
     *  @param[in]  root   is the processor where the values are needed
     */
    virtual void allOwners( hmemo::HArray<PartitionId>& owners, const PartitionId root ) const;

    /** This method returns the owned indexes by this processor.
     *
     *  @param[out] myGlobalIndexes array with localSize 'global' indexes that are owned by this processor
     */
    virtual void getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const;

    /** The following function verifies if the distribution is nothing else than a block
     *  or general block distribution.
     *
     *  @returns the local size of the block distribution if it is one, invalidIndex if it is not
     *
     *  Note: The call of this function might involve communication. It returns invalidIndex on all processors if it is invalidIndex on one.
     *  Note: If it is a block distribution, the distribution of a distributed vector/matrix can be easily reconstructed without a mapping file.
     *
     *  getBlockDistributionSize() != invalidIndex iff isSorted( owners( {0, ..., globalSize-1 }, ascending = true )
     */
    virtual IndexType getBlockDistributionSize() const = 0;

    /**
     * Determine whether the distribution has AnyAddressing.
     *
     * See the documentation for enableAnyAddressing for more information.
     */
    virtual bool hasAnyAddressing() const = 0;

    /**
     * @brief This method sets up local data structures in such a way that afterwards on each
     *        partition/processor it is possible to get any owner/local index for any global index
     *        without communication.
     *
     * Most distributions can do any addressing by closed formulas ( block, cyclic, grid distributions )
     * but general distributions have no information about the ownership of global indexes that are not
     * local. If a general distribution is used in dense matrices for the column distribution it is essential
     * to know for all indexes where the data is mapped to.
     */
    virtual void enableAnyAddressing() const = 0;

    /** Get the local size of any partititon. */

    virtual IndexType getAnyLocalSize( const PartitionId partition ) const = 0;

    /** Get the owner for any global index, same as findOwner but here no communication is guaranteed */

    virtual IndexType getAnyOwner( const IndexType globalIndex ) const = 0;

    /** Get the local index for any global index, owner is not really required but might be helpful if already available. */

    virtual IndexType getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const = 0;

    /** Get the global index for any local index on any partition */

    virtual IndexType getAnyGlobalIndex( const IndexType locaIndex, const PartitionId owner ) const = 0;

    /** This method returns the permutation of global indexes that sorts them by the different owners
     *  (same as bucket sort of array with all owners).
     *
     *  @param[out] offsets local sizes of all partitions as offset array, size is number of partitions + 1
     *  @param[out] local2global contains all global indexes sorted by the owners, is permutation
     *
     *  With the output arrays, the call of getAnyGlobalIndex can be done directly on heterorgeneous arrays
     *
     *  \code
     *     getAnyGlobalIndex( localIndex, owner ) == local2global[ offsets[owner] + localIndex]
     *  \endcode
     *
     *  /code
     *     PartitionId p = ... // some partition
     *     // traverse all global indexes belonging to partition p
     *     for ( IndexType k = offsets[p]; k < offsets[p+1]; ++k )
     *     {
     *         IndexType localIndex = k - offsets[p];
     *         IndexType globalIndex = local2global[k];
     *         ....
     *     }
     *  /endcode
     *
     *  Note: If this distribution is a block distribution, the permutation is the identitiy.
     */
    virtual void getAnyLocal2Global( hmemo::HArray<IndexType>& offsets, hmemo::HArray<IndexType>& local2global ) const;

    /** This method returns the inverse permutation as called by getAnyLocal2Global.
     *
     *  \code
     *     getAnyLocalIndex( globalIndex, owner ) == global2local[globalIndex] - offsets[owner]
     *  \endcode
     */
    virtual void getAnyGlobal2Local( hmemo::HArray<IndexType>& offsets, hmemo::HArray<IndexType>& global2local ) const;

    /**
     * Virtual method to check two distributions for equality.
     *
     * @param[in] other distribution used for comparison with this distribution
     * @returns true if distributions are proven to be equal, false if equality is not guaranteed
     *
     * This method must be implemented by each derived class.
     */
    virtual bool isEqual( const Distribution& other ) const = 0;

    /** Override virtual method Printable::writeAt */

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

    /** Replication of distributed data, e.g. a distributed vector, where all local
     *  values will be replicated on all processors corresponding to their order
     *
     * @tparam     T1           Value type of output data
     * @tparam     T2           Value type of input data
     * @param[out] allValues    array over the global range will contain all values
     * @param[in]  localValues  array over the local range has only values for each partition
     *
     * \code
     *   Partitions:          0        1      2
     *   local indexes:    0, 3, 5    1, 4   2, 6
     *   local Values:     a, b, c    d, e   f, g
     *
     *                     0  1  2  3  4  5  6
     *   allvalues  :      a, d, f, b, e, c, g
     * \endcode
     */
    template<typename T1, typename T2>
    void replicate( T1* allValues, const T2* localValues ) const;

    /** Replication of distributed data, e.g. a distributed dense array where all
     *  local values will be replicated. In contrary to replicate there are n elements for
     *  each distributed element.
     *
     * @tparam     T1           Value type of output data
     * @tparam     T2           Value type of input data
     * @param[out] allValues    array over the global range will contain all values
     * @param[in]  localValues  array over the local range has only values for each partition
     * @param[in]  n            is number of entries in each line of data
     */
    template<typename T1, typename T2>
    void replicateN( T1* allValues, const T2* localValues, const IndexType n ) const;

    /** Replication of distributed data, several entries for each element of the global range
     *
     * @tparam     ValueType     Value type of input and output data
     * @param[out] allValues     array with all values from all partitions
     * @param[in]  localValues   elements available on this partition
     * @param[in]  allOffsets    contains unique offset for each element of the global range
     *
     * Note: the size of allValues must be globalSize + 1, allOffsets[globalSize] is total number of values
     *
     * The offset array is used like in CSR sparse matrix storage for offsets and number of values.
     */
    template<typename ValueType>
    void replicateRagged( ValueType allValues[], const ValueType localValues[], const IndexType allOffsets[] ) const;

protected:

    IndexType mGlobalSize;
    CommunicatorPtr mCommunicator;

    /** Type definition of a function to create a distribution.
     *
     *  @param[in] commPtr is the communicator
     *  @param[in] globalSize number of elements to distribute
     *  @param[in] weight  individual weight for each partition
     *
     *  The weight can be used to determine different load on the partitions depending on
     *  their compute power. E.g. a partition with weight = 4 will get two times the load
     *  of a partition with weight = 2 and four times the load of a partition with weigth = 1.
     */
    typedef Distribution* ( *CreateFn1 )( CommunicatorPtr commPtr, IndexType globalSize, float weight );

    /** Type definition of a function to create a distribution with a connectivity matrix.
     *
     *  @param[in] commPtr is the communicator
     *  @param[in] matrix  is a sparse matrix that provides the connectivity
     *  @param[in] weight  individual weight for each partition
     *
     *  The weight can be used to determine different load on the partitions depending on
     *  their compute power. E.g. a partition with weight = 4 will get two times the load
     *  of a partition with weight = 2 and four times the load of a partition with weigth = 1.
     */
    typedef Distribution* ( *CreateFn2 )( CommunicatorPtr commPtr, const Distributed& matrix, float weight );

    /** Common routine for proving two distributions to be equal
     *
     *  @param[in,out]  isSame will be set to true or false if one general criterion applies
     *  @param[in]      other is the distribution compared with this one
     *  @returns        true if equality or inequality could be proven
     *
     *  Note: if this routine returns false ( nothing proven ), the two distributions have at least
     *        the same global size and the same communicator
     */

    bool proveEquality( bool& isSame, const Distribution& other ) const;

private:

    Distribution(); // disable default constructor

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* ======================================================================== */
/*             Inline methods                                               */
/* ======================================================================== */

IndexType Distribution::getGlobalSize() const
{
    return mGlobalSize;
}

inline bool Distribution::isReplicated() const
{
    return getNumPartitions() == 1;
}

} /* end namespace dmemo */

} /* end namespace scai */
