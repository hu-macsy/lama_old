/**
 * @file CyclicDistribution.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Distribution class to support block-cyclic distributions.
 * @author Lauretta Schubert
 * @date 20.05.2011
 * @since 1.0.0
 */
#ifndef LAMA_CYCLIC_DISTRIBUTION_HPP_
#define LAMA_CYCLIC_DISTRIBUTION_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <lama/distribution/Distribution.hpp>

namespace lama
{

/** For a block cyclic distribution the range is partitioned into chunks of
 * a certain block size and these chunks are assigned in a round-robin manner
 * to the available partitions.
 */
class COMMON_DLL_IMPORTEXPORT CyclicDistribution: public Distribution
{
public:

    /**
     * Construct a new object for a block-cyclic distribution.
     *
     * @param[in] globalSize number of elements to distribution
     * @param[in] chunkSize is the size of the chunk block.
     * @param[in] communicator TODO[doxy] Complete Description.
     */
    CyclicDistribution( const IndexType globalSize, const IndexType chunkSize, const CommunicatorPtr communicator );

    virtual ~CyclicDistribution();

    /**
     * TODO[doxy] Complete Description.
     *
     * @param[in] index TODO[doxy] Complete Description.
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
     * @brief TODO[doxy] Complete Description.
     */
    virtual IndexType getLocalSize() const;

    /**
     * @brief TODO[doxy] Complete Description.
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

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    virtual void computeOwners( const std::vector<IndexType>& requiredIndexes, std::vector<PartitionId>& owners ) const;

    /**
     * @brief TODO[doxy] Complete Description.
     *
     * @param[in] problem TODO[doxy] Complete Description.
     */
    void printDistributionVector( std::string problem ) const;

private:

    CyclicDistribution(); // disable default constructor as it has no global size

    IndexType allGlobal2local( const IndexType globalIndex ) const;

    /** Help routine to get the number of local chunks and info about rmemaining elements */

    void getChunkInfo( IndexType& localChunks, IndexType& extra, const PartitionId rank ) const;

    IndexType mChunkSize; // chunk size of each partition

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

IndexType CyclicDistribution::chunkSize() const
{
    return mChunkSize;
}

}

#endif // LAMA_CYCLIC_DISTRIBUTION_HPP_
