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
#include <scai/partitioning/Partitioning.hpp>

// std
#include <vector>

namespace scai
{

namespace partitioning
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

    /** Implementation of pure method Partitioning::rectangularPartitioning */

    virtual void rectangularPartitioning( hmemo::HArray<PartitionId>& rowMapping, 
                                          hmemo::HArray<PartitionId>& colMapping, 
                                          const lama::_Matrix& matrix, 
                                          const hmemo::HArray<float>& processorWeights ) const;

    using Partitioning::squarePartitioning;

    /** Implementation of pure method Partitioning::squarePartitioning 
     *
     *  Be careful: if the input matrix is distributed the CSR graph will be build on
     *              the MASTER node that calls a routine of the serial software package METIS;
     *              the computed owners will be distributed back to the array newLocalOwners
     *              corresponding to the distribution of the matrix.
     */
    virtual void squarePartitioning( hmemo::HArray<PartitionId>& newLocalOwners,
                                     const lama::_Matrix& matrix,
                                     const hmemo::HArray<float>& processorWeights ) const;

    /** Override Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    /** Static method required for create to use in Partitioning::Register */

    static PartitioningPtr create();

    /** Static method required for Partitioning::Register */

    static std::string createValue();

private:

    dmemo::DistributionPtr computeIt( const dmemo::CommunicatorPtr comm, const lama::_Matrix& matrix, std::vector<float>& weights ) const;

    template<typename weightType>
    void callPartitioning(
        std::vector<IndexType>& partition,
        IndexType& minConstraint,
        IndexType& parts,
        std::vector<weightType>& tpwgts,
        const dmemo::CommunicatorPtr comm,
        const lama::_Matrix& matrix ) const;

    template<typename weightType>
    void checkAndMapWeights(
        std::vector<weightType>& tpwgts,
        std::vector<IndexType>& mapping,
        IndexType& count,
        std::vector<float>& weights,
        IndexType size ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static const char theCreateValue[];

    /** squarePartitioning called for whole matrix on the master node. */

    virtual void squarePartitioningMaster( hmemo::HArray<PartitionId>& newLocalOwners,
                                           const lama::_Matrix& matrix,
                                           const hmemo::HArray<float>& processorWeights ) const;
};

} /* end namespace partitioning */

} /* end namespace scai */
