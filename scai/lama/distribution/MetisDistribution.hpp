/**
 * @file MetisDistribution.hpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief MetisDistribution.hpp
 * @author Lauretta Schubert
 * @date 01.07.2013
 * @since 1.1.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/distribution/GeneralDistribution.hpp>

// local library
#include <scai/lama/matrix/SparseMatrix.hpp>

// std
#include <vector>

namespace scai
{

namespace lama
{

/** Metis creates a distribution with load balance in the pieces and less communication volume
 *  for a dedicated sparse matrix.
 *
 *  MetisDistribution is noncopyable as Distribution is noncopyable
 *
 */

class COMMON_DLL_IMPORTEXPORT MetisDistribution: public GeneralDistribution
{
public:

    /** Construct a block distribution for a number of elements on to the partitions of the passed communicator.
     *
     *  @param[in] comm  used for the partitions onto which elements are distributed.
     *  @param[in] matrix  the matrix the distribution is dedicated for
     *  @param[in] weights  weights for the computational load to the processors
     */
    MetisDistribution( const CommunicatorPtr comm, const Matrix& matrix, std::vector<float>& weights );

    /** Same as above but with individual weight of each processor. */

    MetisDistribution( const CommunicatorPtr comm, const Matrix& matrix, float weight );

    virtual ~MetisDistribution();

    virtual void writeAt( std::ostream& stream ) const;

    static MetisDistribution* create(
        const CommunicatorPtr communicator,
        const IndexType globalSize,
        const float weight = 1.0 );

    static MetisDistribution* create(
        const CommunicatorPtr communicator,
        const Matrix& matrix,
        const float weight = 1.0 );


    virtual std::string getDistributionKind() const
    {
        return "METIS";
    }

private:

    MetisDistribution();

    void computeIt( const CommunicatorPtr comm, const Matrix& matrix, std::vector<float>& weights );

    template<typename weightType>
    void callPartitioning(
        std::vector<IndexType>& partition,
        IndexType& minConstraint,
        IndexType& parts,
        std::vector<weightType>& tpwgts,
        const CommunicatorPtr comm,
        const Matrix& matrix ) const;

    template<typename weightType>
    void checkAndMapWeights(
        std::vector<weightType>& tpwgts,
        std::vector<IndexType>& mapping,
        IndexType& count,
        std::vector<float>& weights,
        IndexType size ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    std    ::vector<float> mWeights;

    void normWeights( std::vector<float>& weights );

    static bool registerCreator(); //!< used in static initialization for registration

    static bool initialized;//!< static initialization used for registration of create in Distribution factory
};

} /* end namespace lama */

} /* end namespace scai */
