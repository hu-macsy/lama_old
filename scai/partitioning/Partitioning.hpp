/**
 * @file Partitioning.hpp
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
#include <scai/utilskernel/HArrayUtils.hpp>

#include <vector>

namespace scai
{

/** The namespace partitioning holds a relevant classes that divide graph into smaller components
 *  and their use for finding a good mapping of sparse matrices to multiple processsors.
 */
namespace partitioning
{

typedef std::shared_ptr<const class Partitioning> PartitioningPtr;

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

    /** Determine a new distribution for a sparse matrix 
     *
     *  @param[in] comm specifies the set of processors onto which the matrix is distributed
     *  @param[in] matrix is the square sparse matrix that should be partitioned
     *  @param[in] weight individual value for each processor 
     */
    virtual dmemo::DistributionPtr partitionIt( 
        const dmemo::CommunicatorPtr comm, 
        const lama::_Matrix& matrix, 
        float weight ) const;

    /** This method determines a new mapping of the local elements of a given distribution
     * 
     *  @param[out]  newOwners contains the new owners of the local rows, newOwners.size() == matrix.getRowDistribution().getLocalSize()
     *  @param[in]   matrix    (sparse) matrix that specifies connectivity 
     *  @param[in]   processorWeights considered weight for each partition
     *
     *  The size of the array processorWeights determines the number of partitions that will be used for the partitioning.
     *  It does not have to be the same number as the size of the communicator used for the distribution of the matrix.
     *
     *  Note: with the new mapping a redistributor can be built that contains the new distribution but also the
     *        communication schedule for redistribution of data structures that have the same mapping
     *
     *  \code
     *     CSRSparseMatrix<double> matrixA( "xyzmatrix.mtx" );
     *     matrixA.redistribute( dist, dist );   // initial distribution
     *     HArray<float> processorWeights( np, 1.0f );
     *     HArray<PartitionId> newOwners;
     *     squarePartitioning( newOwners, matrixA, processorWeights );
     *  \endcode
     *
     *  Note: Internally a CSR graph for the adjacency matrix is built. It is recommended to
     *        use a matrix with same row and column distribution.
     */
    virtual void squarePartitioning( hmemo::HArray<PartitionId>& newOwners,
                                     const lama::_Matrix& matrix,
                                     const hmemo::HArray<float>& processorWeights ) const = 0;

    virtual void squarePartitioningW( hmemo::HArray<PartitionId>& newOwners,
                                      const lama::_Matrix& matrix,
                                      const hmemo::HArray<IndexType>& vertexWeights,
                                      const hmemo::HArray<float>& processorWeights ) const;

    /** This method is a special case of the above one but here the number of the processor
     *  weights are gathered before, i.e. number of new and old partitions remains the same.
     */
    void squarePartitioning( hmemo::HArray<PartitionId>& newLocalOwners,
                             const lama::_Matrix& matrix,
                             const float weight ) const;

    void squarePartitioningW( hmemo::HArray<PartitionId>& newLocalOwners,
                              const lama::_Matrix& matrix,
                              const hmemo::HArray<IndexType>& vertexWeights,
                              const float weight ) const;

    /** Partitioning of rectangular matrix
     *
     *  @param[in]  matrix            might be any distributed matrix
     *  @param[in]  processorWeights  is an array that specifies the desired weight for the partitions
     *  @param[out] rowMapping        array with new mapping for the local rows
     *  @param[out] colMapping        array with new mapping for the local columns
     * 
     *  Notes: 
     *
     *   - the number of (new) partitions is implicitly given by processorWeights.size(),
     *     so this routine might also be called on a serial machine
     *   - rowMapping.size() will be matrix.getRowDistribution.getLocalSize()
     *   - colMapping.size() will be matrix.getColDistribution.getLocalSize()
     *   - this routine does not any remapping here
     */
    virtual void rectangularPartitioning( hmemo::HArray<PartitionId>& rowMapping,
                                          hmemo::HArray<PartitionId>& colMapping,
                                          const lama::_Matrix& matrix,
                                          const hmemo::HArray<float>& processorWeights ) const = 0;


    /** Repartitioning of a rectangular matrix with redistribution.
     *
     *  @param[in,out]  matrix  that will be redistributed
     *  @param[in]      weight  is the weight used by the calling processor
     *
     *  _Matrix must be distributed among the calling processors with a given communicator
     *  that is taken for the redistribution.
     *
     *  Default implementation uses rectangular partititioning by determining the owners and
     *  building a general distribution.
     */
    virtual void rectangularRedistribute( lama::_Matrix& matrix, const float weight ) const;

    /** Override Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    /* The following functions are planned, i.e. rectangular partitiong but either column 
       or row distribution is fixed.

    virtual rowPartitioning( HArray<PartitionId>& rowMapping, 
                             const HArray<PartitionId>& colMapping,
                             const lama::_Matrix& matrix, 
                             const HArray<float>& processorWeights );

    virtual colPartitioning( HArray<PartitionId>& colMapping, 
                             const HArray<PartitionId>& rowMapping,
                             const lama::_Matrix& matrix, 
                             const HArray<float>& processorWeights );
    */

protected:

    static void gatherWeights( std::vector<float>& weights, const float weight, const dmemo::Communicator& comm );

    static void gatherWeights( hmemo::HArray<float>& weights, const float weight, const dmemo::Communicator& comm );

    /** Norm the weights that its sum is exactly 1. */

    static void normWeights( std::vector<float>& weights );

    static void normWeights( hmemo::HArray<float>& weights );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template<typename IdxType>
    static void getDistributionOffsets( hmemo::HArray<IdxType>& offsets, const dmemo::Distribution& dist );

    template<typename ValueType>
    static void gatherAll( ValueType vals[], const ValueType val, const dmemo::Communicator& comm );

    static void normWeights( float weights[], IndexType np );
};

template<typename ValueType>
void Partitioning::gatherAll( ValueType values[], const ValueType value, const dmemo::Communicator& comm )
{
    const PartitionId MASTER = 0;

    PartitionId numPartitions = comm.getSize();

    SCAI_LOG_INFO ( logger, comm << ": gather, my value = " << value )
    comm.gather( values, 1, MASTER, &value );
    SCAI_LOG_INFO ( logger, comm << ": bcast all values" )
    comm.bcast( values, numPartitions, MASTER );
    SCAI_LOG_INFO ( logger, comm << ": bcast done" )
}

template<typename IdxType>
void Partitioning::getDistributionOffsets( hmemo::HArray<IdxType>& offsets, const dmemo::Distribution& dist )
{
    const dmemo::Communicator& comm = dist.getCommunicator();

    const IndexType np = comm.getSize();

    const IdxType mySize = dist.getLocalSize();

    {
        hmemo::WriteOnlyAccess<IdxType> wOffsets( offsets, np + 1 );
        wOffsets.resize( np );
        Partitioning::gatherAll( wOffsets.get(), mySize, dist.getCommunicator() );
    }

    utilskernel::HArrayUtils::scan1( offsets );
}

} /* end namespace partitioning */

} /* end namespace scai */
