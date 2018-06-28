/**
 * @file CyclicDistribution.hpp
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
 * @brief Distribution class to support block-cyclic distributions.
 * @author Lauretta Schubert
 * @date 20.05.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Distribution.hpp>

namespace scai
{

namespace dmemo
{

/** For a block cyclic distribution the range is partitioned into chunks of
 * a certain block size and these chunks are assigned in a round-robin manner
 * to the available partitions.
 */
class COMMON_DLL_IMPORTEXPORT CyclicDistribution:

    public Distribution,
    private Distribution::Register<CyclicDistribution>

{
public:

    /**
     * Construct a new object for a block-cyclic distribution.
     *
     * @param[in] globalSize number of elements to distribution
     * @param[in] chunkSize is the size of the chunk block.
     * @param[in] communicator specifies the set of processors onto which the elements are distributed
     */
    CyclicDistribution(
        const IndexType globalSize,
        const IndexType chunkSize,
        const CommunicatorPtr communicator = Communicator::getCommunicatorPtr() );

    virtual ~CyclicDistribution();

    /**
     * @brief Implementation of pure method Distribution::isLocal
     */
    virtual bool isLocal( const IndexType index ) const;

    /**
     * @brief Computes the owner of the passed globalIndex.
     *
     * @param[in] globalIndex   the global index to compute the owner for.
     * @return                  the owner of the passed globalIndex.
     */
    virtual PartitionId getOwner( const IndexType globalIndex ) const;

    /**
     * @brief Query the number of owned indexes on this processor.
     */
    virtual IndexType getLocalSize() const;

    /** Override default implementation Distribution::getMaxLocalSize() */

    virtual IndexType getMaxLocalSize() const;

    /**
     * @brief get number of elements nb in chunk as defined by Cyclic( nb )
     */
    inline IndexType chunkSize() const;

    /**
     * @brief Returns the number of local chunks
     *
     * @returns number of chunks, also included the remaining chunk that might
     * have less elements than chunkSize
     */
    IndexType getNumLocalChunks() const;

    /**
     * @brief Returns the number of chunks for a given partition
     *
     * @returns number of chunks for any partition, includes the remaining chunk
     * that might have less elements than chunkSize
     */
    IndexType getNumChunks( const PartitionId partition ) const;

    /**
     * @brief Returns the global number of chunks
     *
     * @returns number of chunks, also included the remaining chunk that might
     * have less elements than chunkSize
     */
    IndexType getNumTotalChunks() const;

    /**
     * @brief Returns the size of the passed partition
     *
     * @returns size of the passed partition
     */
    IndexType getPartitionSize( const PartitionId partition ) const;

    virtual IndexType local2global( const IndexType localIndex ) const;

    virtual IndexType global2local( const IndexType globalIndex ) const;

    /** Implementation of pure function Distribution::getBlockDistributionSize.
     *
     *  A cyclic distribution Cyclic( globalSize, chunkSize ) is a block distribution
     *  iff chunkSize * numPartitions <= globalSize, i.e. each processor has maximal one chunk
     */
    virtual IndexType getBlockDistributionSize() const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    virtual void computeOwners( hmemo::HArray<PartitionId>& owners, const hmemo::HArray<IndexType>& indexes ) const;

    /** Override Distribution::getOwnedIndexes with more efficient version. */

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

    /** Static method required for create to use in Distribution::Register */

    static Distribution* create( const DistributionArguments args );

    /** Static method required for Distribution::Register */

    static std::string createValue();

    virtual const char* getKind() const
    {
        return getId();
    }

    static const char* getId()
    {
        return "CYCLIC";
    }

private:

    CyclicDistribution(); // disable default constructor as it has no global size

    IndexType allGlobal2local( const IndexType globalIndex ) const;

    /** Help routine to get the number of local chunks and info about rmemaining elements */

    void getChunkInfo( IndexType& localChunks, IndexType& extra, const PartitionId rank ) const;

    IndexType mChunkSize; // chunk size of each partition

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

IndexType CyclicDistribution::chunkSize() const
{
    return mChunkSize;
}

} /* end namespace dmemo */

} /* end namespace scai */
