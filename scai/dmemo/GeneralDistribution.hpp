/**
 * @file GeneralDistribution.hpp
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
 * @brief GeneralDistribution.hpp
 * @author brandes
 * @date 25.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/hmemo/HArray.hpp>
#include <scai/dmemo/Distribution.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>

// std
#include <unordered_map>

namespace scai
{

namespace dmemo
{

/** A general distribution allows to map a global range of values
 to the partitions of a communicator completely arbitrarily.

 Each partition has the information which values it holds.
 A partition has no information where to find values
 not owned by itself.
 */

class COMMON_DLL_IMPORTEXPORT GeneralDistribution: public Distribution
{
public:

    /** Constructor of a general distribution where each processor knows it indexes.
     *
     *  \param globalSize is the size of the distributed range
     *  \param myGlobalIndexes contains all indexes of range owned by this processor
     *  \param communicator specifies the set on processors for the distribution
     *
     *  Important: each global index from 0 to globalSize-1 must appear exactly once in
     *  the vector myGlobalIndexes on one partition.
     *
     *  Note: the (large) arrays for enable any addressing are not built by default.
     */
    GeneralDistribution(
        const IndexType globalSize,
        hmemo::HArray<IndexType> myGlobalIndexes,
        const CommunicatorPtr communicator = Communicator::getCommunicatorPtr() );

    virtual ~GeneralDistribution();

    /** Implementation of pure method Distribution::isLocal */

    virtual bool isLocal( const IndexType globalIndex ) const;

    /** Implementation of pure method Distribution::getLocalSize */

    virtual IndexType getLocalSize() const;

    /** Implementation of pure method Distribution::local2Global */

    virtual IndexType local2Global( const IndexType localIndex ) const;

    /** Implementation of pure method Distribution::global2Local */

    virtual IndexType global2Local( const IndexType globalIndex ) const;

    /** Override Distribution::global2LocalV with a more efficient version */

    virtual void global2LocalV( hmemo::HArray<IndexType>& localIndexes, const hmemo::HArray<IndexType>& globalIndexes ) const;

    /** Implementation of pure function Distribution::getBlockDistributionSize.
     *
     *  Each processor must have a contiguous part of indexes and their order
     *  must match the rank order of processors.
     */
    virtual IndexType getBlockDistributionSize() const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** Override Distribution::computeOwners with more efficient version.
     *
     *  This routine uses a block distributed vector that holds the owners, i.e. each processor
     *  knows the owners of a corresponding contiguous section.
     *  For the queried indexes each processor queries the owners from the processors that holds
     *  the corresponding information.
     */

    virtual void computeOwners( hmemo::HArray<PartitionId>& owners, const hmemo::HArray<IndexType>& indexes ) const;

    /** Override the default implementation of Distribution::allOwners */

    virtual void allOwners( hmemo::HArray<PartitionId>& owners, const PartitionId root ) const;

    /** Override Distribution::getOwnedIndexes */

    virtual void getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const;

    /** Implementation of pure method Distribution::hasAnyAddressing */

    virtual bool hasAnyAddressing() const;

    /** Implementation of pure method Distribution::enableAnyAddressing */

    virtual void enableAnyAddressing() const;

    /** Implementation of pure method Distribution::getAnyLocalSize */

    virtual IndexType getAnyLocalSize( const PartitionId partition ) const;

    /** Implementation of pure method Distribution::getAnyOwner */

    virtual PartitionId getAnyOwner( const IndexType globalIndex ) const;

    /** Implementation of pure method Distribution::getAnyLocalIndex */

    virtual IndexType getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const;

    /** Implementation of pure method Distribution::getAnyGlobalIndex */

    virtual IndexType getAnyGlobalIndex( const IndexType localIndex, const PartitionId owner ) const;

    /** This method returns the array that contains for this processors all owned indexes. */

    inline const hmemo::HArray<IndexType>& getMyIndexes() const;

    /** This method returns the array of owners that this processor knowns about */

    inline const hmemo::HArray<IndexType>& getMyBlockDistributedOwners() const;

    /** Implementation of pure method Distribution::getKind */

    virtual inline const char* getKind() const;

protected:

    static const char theCreateValue[];

    hmemo::HArray<IndexType> mLocal2Global;   //!< for each local index its global index, entries are sorted

    /** 
     * Hash map to find for a global index the local index (if owned by this processor)
     *
     * Note: this solution is more efficient than a binary search in mLocal2Global.
     */
    std::unordered_map<IndexType, IndexType> mGlobal2Local;

    /** Block distributed array of owners (used for compute owners and for sanity check)
     *
     *  In the global view the array mBlockDistributedOwners[i] contains the owner of index i.
     */
    hmemo::HArray<PartitionId> mBlockDistributedOwners;

    // the following arrays will only be available if enableAnyAddressing has been called
    // Note: if set the array mGlobal2Local is no more needed

    mutable hmemo::HArray<PartitionId> mAllOwners;         // will have globalSize entries on each processor
    mutable hmemo::HArray<IndexType> mAllLocalOffsets;     // local size on each partition
    mutable hmemo::HArray<IndexType> mAllLocal2Global;     // sorts elements into buckets
    mutable hmemo::HArray<IndexType> mAllGlobal2Local;     // sorts elements into buckets

    // Example
    // index       0    1    2    3   4    5    6    7   8   9   10   11   12
    // mOwners:    0    1    2    0   2    0    1    0   0   1    1    2    2
    // Offsets:    0                       5                 9                    13
    // perm   :    0    3    5    7   8    1    6    9  10   2    4   11   12     local2Global
    // perm'  :    0    5    9    1  10    2    6    3   4   7    8   11   12     global2Local
    //
    // Note: perm is identity iff we have a block distribution

private:

    GeneralDistribution();

    GeneralDistribution& operator=( const GeneralDistribution& other ) = delete;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Help routine that computes the map from global indexes to local indexes */

    void fillIndexMap();

    /** Determine block-distributed ownership, i.e. mBlockDistributedOwners */
 
    void setBlockDistributedOwners();

    static void computeBlockDistributedOwners( hmemo::HArray<IndexType>& blockDistributedOwners,
                                               const IndexType globalSize,
                                               const hmemo::HArray<IndexType>& ownedIndexes,
                                               CommunicatorPtr comm );
};

/** Constructor of a general distribution here as a function for convenience
 *
 *  @param[in] globalSize is the size of the distributed range
 *  @param[in] myGlobalIndexes contains all indexes of range owned by this processor
 *  @param[in] comm specifies the set of processors for the distribution
 *
 *  The method must be called by all processors of comm at the same time. 
 */
std::shared_ptr<GeneralDistribution> generalDistribution( 
    const IndexType globalSize,
    hmemo::HArray<IndexType> myGlobalIndexes,
    const CommunicatorPtr comm = Communicator::getCommunicatorPtr() );

/** This function creates a general distribution by an array containing the owner for each element
 *
 *  @param[in] owners array with ower for each element, 0 <= owners[i] < communicator->size()
 *  @param[in] root is the processor that has the valid copy of owners, must be same value on all processors
 *  @param[in] comm specifies the processor set used for the distribution
 *
 *  The method must be called by all processors of comm at the same time. 
 */
std::shared_ptr<GeneralDistribution> generalDistributionByOwners( 
    const hmemo::HArray<PartitionId>& owners, 
    const PartitionId root, 
    CommunicatorPtr comm );

/** This function creates a new general distributuion by an existing one and a mapping of the
 *  local elements to new owners. 
 *
 *  @param[in] dist is the actual distribution
 *  @param[in] newOwners contains the new owner for each of the local indexes.
 */
std::shared_ptr<GeneralDistribution> generalDistributionNew( 
    const Distribution& dist,
    const hmemo::HArray<PartitionId>& newOwners );

/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods                                         */
/* ------------------------------------------------------------------------- */

const hmemo::HArray<IndexType>& GeneralDistribution::getMyIndexes() const
{
    return mLocal2Global;
}

const hmemo::HArray<PartitionId>& GeneralDistribution::getMyBlockDistributedOwners() const
{
    return mBlockDistributedOwners;
}

const char* GeneralDistribution::getKind() const
{
    return theCreateValue;
}

} /* end namespace dmemo */

} /* end namespace scai */
