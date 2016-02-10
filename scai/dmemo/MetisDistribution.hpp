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
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/Distributed.hpp>

// std
#include <vector>

namespace scai
{

namespace dmemo
{

/** Metis creates a distribution with load balance in the pieces and less communication volume
 *  for a dedicated sparse matrix.
 *
 *  MetisDistribution is noncopyable as Distribution is noncopyable
 *
 */

class COMMON_DLL_IMPORTEXPORT MetisDistribution: 

    public GeneralDistribution,
    private Distribution::Register<MetisDistribution>

{
public:

    /** Construct a new general distribution for a number of elements on to the partitions of the passed communicator.
     *
     *  @param[in] comm  used for the partitions onto which elements are distributed.
     *  @param[in] matrix is an object whose size and connectivity is used for the new distribution
     *  @param[in] weights  weights for the computational load to the processors
     */
    MetisDistribution( const CommunicatorPtr comm, const Distributed& matrix, std::vector<float>& weights );

    /** Same as above but with individual weight of each processor. */

    MetisDistribution( const CommunicatorPtr comm, const Distributed& matrix, float weight );

    virtual ~MetisDistribution();

    virtual void writeAt( std::ostream& stream ) const;

    /** Static method required for create to use in Distribution::Register */

    static Distribution* create( const DistributionArguments args );

    /** Static method required for Distribution::Register */

    static std::string createValue();

private:

    MetisDistribution();

    void computeIt( const CommunicatorPtr comm, const Distributed& matrix, std::vector<float>& weights );

    template<typename weightType>
    void callPartitioning(
        std::vector<IndexType>& partition,
        IndexType& minConstraint,
        IndexType& parts,
        std::vector<weightType>& tpwgts,
        const CommunicatorPtr comm,
        const Distributed& matrix ) const;

    template<typename weightType>
    void checkAndMapWeights(
        std::vector<weightType>& tpwgts,
        std::vector<IndexType>& mapping,
        IndexType& count,
        std::vector<float>& weights,
        IndexType size ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    std::vector<float> mWeights;  //!< The weights of all partitions, mWeights.size() == mComm.size()

    /** Norm the weights that its sum is exactly 1. */

    static void normWeights( std::vector<float>& weights );
};

} /* end namespace dmemo */

} /* end namespace scai */
