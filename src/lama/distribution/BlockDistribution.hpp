/**
 * @file BlockDistribution.hpp
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
 * @brief BlockDistribution.hpp
 * @author thomas Brandes
 * @date 18.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_BLOCKDISTRIBUTION_HPP_
#define LAMA_BLOCKDISTRIBUTION_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/distribution/Distribution.hpp>

namespace lama
{

/** Block distribution in contiguous chunks of same size for each partition.
 *
 *  BlockDistribution is noncopyable as Distribution is noncopyable
 *
 */
class LAMA_DLL_IMPORTEXPORT BlockDistribution: public Distribution
{
public:

    /** Construct a block distribution for a number of elements on to the partitions of the passed communicator.
     *
     *  @param[in] globalSize   number of elements to distribute
     *  @param[in] communicator used for the partitions onto which elements are distributed.
     */
    BlockDistribution( const IndexType globalSize, const CommunicatorPtr communicator );

    virtual ~BlockDistribution();

    /** Static method that allows to compute the range of a block distribution for
     *  arbitrary rank and communicators.
     */

    static void getRange(
        IndexType& lb,
        IndexType& ub,
        const IndexType n,
        const PartitionId rank,
        const PartitionId size );

    virtual bool isLocal( const IndexType index ) const;

    /**
     * @brief Computes the owner of the passed globalIndex.
     *
     * @param[in] globalIndex   the global index to compute the owner for.
     * @return                  the owner of the passed globalIndex.
     */
    virtual PartitionId getOwner( const IndexType globalIndex ) const;

    virtual IndexType getLocalSize() const;

    virtual IndexType local2global( const IndexType localIndex ) const;

    virtual IndexType global2local( const IndexType globalIndex ) const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** Method to compute directly the owners of global indexes without
     *  any communication.
     */
    virtual void computeOwners( const std::vector<IndexType>& requiredIndexes, std::vector<PartitionId>& owners ) const;

    /** Static method to construct a new block distribution. */

    static DistributionPtr create( const IndexType globalSize, const CommunicatorPtr communicator );

    void printDistributionVector( std::string problem ) const;

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private:

    BlockDistribution(); // disable default constructor as it has no size

    IndexType mBlockSize;//!< block size of each partition
    IndexType lb, ub;//!< local range of full size in global values
};

}

#endif // LAMA_BLOCKDISTRIBUTION_HPP_
