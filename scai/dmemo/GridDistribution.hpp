/**
 * @file GridDistribution.hpp
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
 * @brief GridDistribution.hpp
 * @author Thomas Brandes
 * @date 30.01.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/common/Grid.hpp>

// base classes
#include <scai/dmemo/Distribution.hpp>

namespace scai
{

namespace dmemo
{

/** Grid distribution stands for a block distribution in multiple dimensions.
 *
 */
class COMMON_DLL_IMPORTEXPORT GridDistribution:

    public Distribution,
    private Distribution::Register<GridDistribution>

{
public:

    /** Construct a grid distribution
     *
     *  @param[in] globalGrid specifies the problem grid to be distributed
     *  @param[in] communicator used for the partitions onto which grid is distributed
     *  @param[in] procGrid specifies the processor grid
     *
     *  Note: procGrid.size() must be less or equal than communicator.size()
     */
    GridDistribution( const common::Grid& globalGrid, const CommunicatorPtr communicator, const common::Grid& procGrid );

    virtual ~GridDistribution();

    virtual bool isLocal( const IndexType index ) const;

    virtual PartitionId findOwner( const IndexType globalIndex ) const;

    virtual IndexType getLocalSize() const;

    /** Override default implementation Distribution::getMaxLocalSize() */

    virtual IndexType getMaxLocalSize() const;

    virtual IndexType local2global( const IndexType localIndex ) const;

    virtual IndexType global2local( const IndexType globalIndex ) const;

    /** Implementation of pure function Distribution::getBlockDistributionSize
     *
     *  A grid distribution is block distributed iff only the first dimension is distributed
     */

    virtual IndexType getBlockDistributionSize() const;

    virtual bool isEqual( const Distribution& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** Override Distribution::computeOwners with more efficient version. */

    virtual void computeOwners( hmemo::HArray<PartitionId>& owners, const hmemo::HArray<IndexType>& indexes ) const;

    /** Override Distribution::getOwnedIndexes with more efficient version. */

    virtual void getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const;

    /** Static method required for create to use in Distribution::Register */

    static Distribution* create( const DistributionArguments args );

    /** Static method required for Distribution::Register */

    static std::string createValue();

    virtual const char* getKind() const
    {
        return getId();
    }

    static const char* getId()
    {
        return "GRID";
    }

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    GridDistribution(); // disable default constructor as it has no size

    common::Grid mGlobalGrid;
    common::Grid mLocalGrid;

    common::Grid mProcGrid;

    // block distribution in each dimension

    IndexType mRank[SCAI_GRID_MAX_DIMENSION];        //!< rank of this processor in each dimension
    IndexType mBlockSize[ SCAI_GRID_MAX_DIMENSION];  //!< block size of each dimension
    IndexType mLB[ SCAI_GRID_MAX_DIMENSION];
    IndexType mUB[ SCAI_GRID_MAX_DIMENSION];
};

} /* end namespace dmemo */

} /* end namespace scai */
