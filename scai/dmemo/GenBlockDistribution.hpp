/**
 * @file GenBlockDistribution.hpp
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
 * @brief General block distribution where partitions might have different sizes.
 * @author Thomas Brandes
 * @date 18.03.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Distribution.hpp>

#include <memory>

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

    /** Construct a general block distribution by an offset array with the sizes
     *
     *  @param[in] offsets is the offset array, must have comm->getSize() + 1 entries
     *  @param[in] comm specifies the communicator used for this distribution
     *
     *  Note: All processors must call this constructor with the 'same' values.
     *
     *  The global size is given by the entry offsets[size].
     */
    GenBlockDistribution(
        std::unique_ptr<IndexType[]>&& offsets,
        CommunicatorPtr comm = Communicator::getCommunicatorPtr() );

    virtual ~GenBlockDistribution();

    /** Get the local range of the calling partition.
     *
     *  @param[out] lb, ub is the local range, i.e all elements i with lb <= i < ub
     *
     *  Note: lb == ub stands for zero size, ub < lb can never happen
     *
     *  Be careful: older version returned ub with lb <= i <= ub
     */

    void getLocalRange( IndexType& lb, IndexType& ub ) const;

    /** Implemenation of pure method Distribution::isLocal */

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

    /** Implementation of pure function Distribution::getBlockDistributionSize, here same as getLocalSize */

    virtual IndexType getBlockDistributionSize() const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** Override Distribution::computeOwners
     *
     *  Each processor knowns the sizes of each partition and can therefore compute
     *  owners without any communication.
     */
    virtual void computeOwners( hmemo::HArray<PartitionId>& owners, const hmemo::HArray<IndexType>& indexes ) const;

    /** Override Distribution::getOwnedIndexes with more efficient version. */

    virtual void getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const;

    /** Implementation of pure method Distribution::hasAnyAddressing */
    virtual bool hasAnyAddressing() const;

    /** Implementation of pure method Distribution::enableAnyAddressing */

    virtual void enableAnyAddressing() const;

    /** Implementation of pure method Distribution::getAnyLocalSize */

    virtual IndexType getAnyLocalSize( const PartitionId partition ) const;

    /** Implementation of pure method Distribution::getAnyOwner */

    virtual PartitionId getAnyOwner( const IndexType globalIndex ) const;

    /** Implementation of pure method Distribution::getAnyLocalIndex */

    virtual IndexType getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const;

    /** Implementation of pure method Distribution::getAnyGlobalIndex */

    virtual IndexType getAnyGlobalIndex( const IndexType localIndex, const PartitionId owner ) const;

    /** Static method required for create to use in Distribution::Register */

    static DistributionPtr create( const DistributionArguments args );

    /** Static method required for Distribution::Register */

    static std::string createValue();

    /** Implementation of pure method Distribution::getKind */

    virtual inline const char* getKind() const;

    /** Unique identification for this derived distribution class. */

    static inline const char* getId();

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    GenBlockDistribution(); // disable default destructor

    std::unique_ptr<IndexType[]> mOffsets;  //!< offset for each partition

    // this processor owns mLB, ..., mUB - 1

    IndexType mLB;
    IndexType mUB;   //!< local range of full size in global values
};

const char* GenBlockDistribution::getKind() const
{
    return getId();
}

const char* GenBlockDistribution::getId()
{
    return "GEN_BLOCK";
}

/** Construct a general block distribution by individual local sizes
 *
 *  @param[in] localSize is the number of elements owned by this processor
 *  @param[in] comm specifies the communicator used for this distribution
 *
 */
std::shared_ptr<GenBlockDistribution> genBlockDistribution( const IndexType localSize,
                                                            CommunicatorPtr comm = Communicator::getCommunicatorPtr() );

/** Construct a general block distribution by a vector of localSizes.
 *
 *  @param[in] localSizes contains the sizes for each partition, must be same for all processors
 *  @param[in] communicator specifies the communicator used for this distribution
 *
 *  Note: All partitions must call this constructor with the 'same' arguments.
 *
 *  This constructor does not involve any global communication.
 */
std::shared_ptr<GenBlockDistribution> genBlockDistribution( const std::vector<IndexType>& localSizes,
                                                            CommunicatorPtr comm = Communicator::getCommunicatorPtr() );

/** Generate a general block distribution by a weight so that each partition
 *  will have a size corresponding to its weight.
 *
 *  @param globalSize is the number of elements to distribute, must be same for all processors
 *  @param weight is the weight of this partition (must not be negative)
 *  @param comm specifies the communicator used for this distribution
 */
std::shared_ptr<GenBlockDistribution> genBlockDistributionByWeight( const IndexType globalSize,
                                                                    const float weight,
                                                                    const CommunicatorPtr comm = Communicator::getCommunicatorPtr() );

} /* end namespace dmemo */

} /* end namespace scai */
