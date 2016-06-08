/**
 * @file MetisDistribution.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @endlicense
 *
 * @brief MetisDistribution.hpp
 * @author Lauretta Schubert
 * @date 01.07.2013
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

    virtual const char* getKind() const
    {
        return theCreateValue;
    }

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

    static const char theCreateValue[];
};

} /* end namespace dmemo */

} /* end namespace scai */
