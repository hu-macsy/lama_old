/**
 * @file GenBlockDistribution.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief General block distribution where partitions might have different sizes.
 * @author Thomas Brandes
 * @date 18.03.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Distribution.hpp>

// internal scai libraries
#include <scai/common/unique_ptr.hpp>

namespace scai
{

namespace dmemo
{

/** Derived distribution class for general block distributions.
 *
 *  For the general block distribution a local block size is given
 *  for each partition.
 *
 *  GenBlockDistribution is like all Distribution classes non-copyable.
 *  Shared pointers might be used to use the same distribution for several
 *  objects.
 */

class COMMON_DLL_IMPORTEXPORT GenBlockDistribution:

    public Distribution,
    private Distribution::Register<GenBlockDistribution>

{
public:


    /** Construct a general block distribution by a vector of localSizes.
     *
     *  @param[in] globalSize is the number of elements to distribute
     *  @param[in] localSizes contains the sizes for each partition
     *  @param[in] communicator specifies the communicator used for this distribution
     *
     *  Note: All partitions must call this constructor with the 'same' arguments.
     */

    GenBlockDistribution(
        const IndexType globalSize,
        const std::vector<IndexType>& localSizes,
        const CommunicatorPtr communicator );

    /** Construct a general block distribution by an interval.
     *
     *  @param[in] globalSize is the number of elements to distribute
     *  @param[in] firstGlobalIdx is the smallest global index in partition
     *  @param[in] lastGlobalIdx is the largest global index in partition
     *  @param[in] communicator specifies the communicator used for this distribution
     *
     */

    GenBlockDistribution(
        const IndexType globalSize,
        const IndexType firstGlobalIdx,
        const IndexType lastGlobalIdx,
        const CommunicatorPtr communicator );

    /** Construct a general block distribution by individual localSize
     *
     *  @param[in] globalSize is the number of elements to distribute
     *  @param[in] localSize is the number of elements for this partition
     *  @param[in] communicator specifies the communicator used for this distribution
     *
     */

    GenBlockDistribution( const IndexType globalSize, const IndexType localSize, const CommunicatorPtr communicator );

    /** Construct a general block distribution by a weight so that each partition
     *  will have a size corresponding to its weight.
     *
     *  @param globalSize is the number of elements to distribute
     *  @param weight is the weight of this partition (must be positive)
     *  @param communicator specifies the communicator used for this distribution
     *
     *  Note: GenBlockDistribution( globalSize, localSize, comm) is the same
     *        as GenBlockDistribution( globalSize, (float) localSize, comm )
     */

    GenBlockDistribution( const IndexType globalSize, const float weight, const CommunicatorPtr communicator );

    virtual ~GenBlockDistribution();

    /** Get the local range of the calling partition. */

    void getLocalRange( IndexType& lb, IndexType& ub ) const;

    // Implementation of abstract methods, override default virtual methods

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

    /** Override Distribution::computeOwners
     *
     *  Each processor knowns the sizes of each partition and can therefore compute
     *  owners without any communication.
     */
    virtual void computeOwners( const std::vector<IndexType>& requiredIndexes, std::vector<PartitionId>& owners ) const;

    void printDistributionVector( std::string name ) const;

    /** Static method required for create to use in Distribution::Register */

    static Distribution* create( const DistributionArguments args );

    /** Static method required for Distribution::Register */

    static std::string createValue();

    virtual const char* getKind() const
    {
        return theCreateValue;
    }

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    static const char theCreateValue[];

    void setOffsets( const IndexType rank, const IndexType numPartitions, const IndexType localSizes[] );

    void setOffsets( const IndexType rank, const IndexType numPartitions, const IndexType mySize );

    GenBlockDistribution(); // disable default destructor

    common::scoped_array<IndexType> mOffsets;//!< offset for each partition

    IndexType mLB, mUB;//!< local range of full size in global values
};

} /* end namespace dmemo */

} /* end namespace scai */
