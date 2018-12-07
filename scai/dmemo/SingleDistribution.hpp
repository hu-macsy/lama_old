/**
 * @file SingleDistribution.hpp
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
 * @brief SingleDistribution.hpp
 * @author thomas Brandes
 * @date 30.01.2017
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

/** Single distribution is a distribution where all data is mapped to one single processor.
 *
 *  In contrary to a replicated distribution where each processor is owner of all data a single
 *  distribution is given if only one processor has all data.
 *
 */
class COMMON_DLL_IMPORTEXPORT SingleDistribution:

    public Distribution,
    private Distribution::Register<SingleDistribution>

{
public:

    /** Construct a single distribution for a number of elements on to the partitions of the passed communicator.
     *
     *  @param[in] globalSize   number of elements to distribute
     *  @param[in] communicator used for the partitions onto which elements are distributed.
     *  @param[in] owner        partition id that owns all data
     */
    SingleDistribution( const IndexType globalSize, const CommunicatorPtr communicator, const PartitionId owner );

    virtual ~SingleDistribution();

    /** Implementation of pure method Distribution::isLocal */

    virtual bool isLocal( const IndexType index ) const;

    /** Override the default implementation Distribution::findOwner of Distribution with solution without communication */

    virtual PartitionId findOwner( const IndexType globalIndex ) const;

    /** Implementation of pure method Distribution::getLocalSize
     *
     *  Returns global size for the owner partition, all other return 0
     */
    virtual IndexType getLocalSize() const;

    /** Override default implementation Distribution::getMaxLocalSize() */

    virtual IndexType getMaxLocalSize() const;

    /** Implementation of pure method Distribution::getBlockDistributionSize
     *
     *  Note: a single distribution is a special case of a general block distribution.
     */
    virtual IndexType getBlockDistributionSize() const;

    /** Implementation of pure method Distribution::local2global */

    virtual IndexType local2global( const IndexType localIndex ) const;

    /** Implementation of pure method Distribution::global2local */

    virtual IndexType global2local( const IndexType globalIndex ) const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** Override Distribution::computeOwners with more efficient version. */

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

    virtual const char* getKind() const
    {
        return getId();
    }

    static const char* getId()
    {
        return "SINGLE";
    }

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    SingleDistribution(); // disable default constructor as it has no size

    PartitionId mOwner;   //!< Partition id for the owner process.
};

} /* end namespace dmemo */

} /* end namespace scai */
