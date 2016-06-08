/**
 * @file BlockDistribution.hpp
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
 * @brief BlockDistribution.hpp
 * @author thomas Brandes
 * @date 18.03.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/dmemo/Distribution.hpp>

namespace scai
{

namespace dmemo
{

/** Block distribution in contiguous chunks of same size for each partition.
 *
 *  BlockDistribution is noncopyable as Distribution is noncopyable
 *
 */
class COMMON_DLL_IMPORTEXPORT BlockDistribution:

    public Distribution,
    private Distribution::Register<BlockDistribution>

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

    void printDistributionVector( std::string problem ) const;

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

    BlockDistribution(); // disable default constructor as it has no size

    static const char theCreateValue[];

    IndexType mBlockSize;//!< block size of each partition
    IndexType mLB;//!< lower bound value of local range
    IndexType mUB;//!< upper bound value of local range
};

} /* end namespace dmemo */

} /* end namespace scai */
