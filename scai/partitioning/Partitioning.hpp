/**
 * @file Partitioning.hpp
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
 * @brief Definition of abstract base class for partitioning methods like Metis
 * @author Thomas Brandes
 * @date 18.07.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/NonCopyable.hpp>
#include <scai/common/Printable.hpp>
#include <scai/common/Factory.hpp>
#include <scai/common/Factory.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/lama/matrix/_Matrix.hpp>

#include <vector>

namespace scai
{

/** The namespace partitioning holds a relevant classes that divide graph into smaller components
 *  and their use for finding a good mapping of sparse matrices to multiple processsors.
 */
namespace partitioning
{

typedef common::shared_ptr<const class Partitioning> PartitioningPtr;

/** Abstract base class for partititioners
 *
 *  A (graph) partitioner computes a new distribution for a given distribution that takes its connections
 *  (edges) into account.
 * 
 *  Typically, the graph is being partitioned so that a computational task can be correspondingly divided up
 *  the processors of the communicator. Thus, it is important that each part of the partition has roughly the 
 *  same number of elements, and that "connections" between distinct partitions are limited, as much as possible, 
 *  to reduce the amount of interprocess communication necessary. 
 *
 *  Default and copy constructor are not available for this class (noncopyable).
 */
class COMMON_DLL_IMPORTEXPORT Partitioning:

    public common::Factory<std::string, PartitioningPtr>,
    public common::Printable,
    private common::NonCopyable
{

public:

    /** Constructor for a partitioning.
     *
     */
    Partitioning();

    /** Destructor of partitionng */

    virtual ~Partitioning();

    /** Pure method that must be implemented by each partitioning class */

    virtual dmemo::DistributionPtr partitionIt( const dmemo::CommunicatorPtr comm, const lama::Matrix& matrix, float weight ) const = 0;

    /** Partitioning of rectangular matrix
     *
     *  @param[in]  matrix            might be any distributed matrix
     *  @param[in]  processorWeights  is an array that specifies the desired weight for the partitions
     *  @param[out] rowMapping        array with new mapping for the local rows
     *  @param[out] colMapping        array with new mapping for the local columns
     * 
     *  Notes: 
     *
     *   - the number of partitions is implicitly given by processorWeights.size()
     *     (so this routine might also be called on a serial machine)
     *   - rowMapping.size() will be matrix.getRowDistribution.getLocalSize()
     *   - colMapping.size() will be matrix.getColDistribution.getLocalSize()
     *   - this routine does not any remapping here
     */
    virtual void rectangularPartitioning( hmemo::HArray<PartitionId>& rowMapping,
                                          hmemo::HArray<PartitionId>& colMapping,
                                          const lama::Matrix& matrix,
                                          const hmemo::HArray<float>& processorWeights ) const = 0;

    /** Repartitioning of a rectangular matrix with redistribution.
     *
     *  @param[in,out]  matrix  that will be redistributed
     *  @param[in]      weight  is the weight used by the calling processor
     *
     *  Matrix must be distributed among the calling processors with a given communicator
     *  that is taken for the redistribution.
     *
     *  Default implementation uses rectangular partititioning by determining the owners and
     *  building a general distribution.
     */
    virtual void rectangularRedistribute( lama::Matrix& matrix, const float weight ) const;

    /** Override Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

protected:

    static void gatherWeights( std::vector<float>& weights, const float weight, const dmemo::Communicator& comm );

    /** Norm the weights that its sum is exactly 1. */

    static void normWeights( std::vector<float>& weights );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace partitioning */

} /* end namespace scai */
