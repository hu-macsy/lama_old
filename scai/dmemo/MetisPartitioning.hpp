/**
 * @file MetisPartitioning.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief MetisPartitioning.hpp
 * @author Thomas Brandes
 * @date 18.08.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Partitioning.hpp>

// std
#include <vector>

namespace scai
{

namespace dmemo
{

/** Metis creates a distribution with load balance in the pieces and less communication volume
 *  for a dedicated sparse matrix.
 *
 *  MetisPartitioning is noncopyable as Partitioning is noncopyable
 *
 */

class COMMON_DLL_IMPORTEXPORT MetisPartitioning:

    public Partitioning,
    private Partitioning::Register<MetisPartitioning>

{
public:

    /** Constructor of an object that can partitition 'serial' graph data
     *
     */
    MetisPartitioning();

    virtual ~MetisPartitioning();

    /** Implementation of pure method Partitioning::partitionIt */

    virtual DistributionPtr partitionIt( const CommunicatorPtr comm, const Distributed& matrix, float weight ) const;

    /** Override Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    /** Static method required for create to use in Partitioning::Register */

    static PartitioningPtr create();

    /** Static method required for Partitioning::Register */

    static std::string createValue();

private:

    DistributionPtr computeIt( const CommunicatorPtr comm, const Distributed& matrix, std::vector<float>& weights ) const;

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

    static const char theCreateValue[];
};

} /* end namespace dmemo */

} /* end namespace scai */
