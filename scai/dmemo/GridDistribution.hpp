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
 *  For a one-dimensional grid, this distribution is exactly the same as a block distribution.
 *  Compared to a GeneralDistribution that might also be used for the distributon of matrices and vectors
 *  it offers some advantages:
 *
 *  - ownership can be calculated by closed formulas, i.e. each processor can determine directly without
 *    communication the owning processor of each grid element
 *  - This distribution delivers a 'local' grid with the local dimensions on each processor that might
 *    be used to allocate local data for distributed grids
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

    /** Construct a grid distribution
     *
     *  @param[in] globalGrid specifies the problem grid to be distributed
     *  @param[in] communicator used for the partitions onto which grid is distributed
     *
     *  Note: Here the processor grid is determined by a good factorization of the available processors in the communicator.
     *  The size of the grid can be set by the environment variable SCAI_NP, e.g. SCAI_NP=3x3x3
     */
    GridDistribution( const common::Grid& globalGrid, const CommunicatorPtr communicator );

    virtual ~GridDistribution();

    virtual bool isLocal( const IndexType index ) const;

    /** Overrides the default implementation of Distribution::findOwner 
     *
     *  Note: In contrary to the default implemantion here no communication is needed. 
     */
    virtual PartitionId findOwner( const IndexType globalIndex ) const;

    /** Implementation of pure method Distribution::getLocalSize */

    virtual IndexType getLocalSize() const;

    /** Get access to the local grid, helpful for traversing */

    const common::Grid& getLocalGrid() const;

    /** Get access to the global grid, helpful for traversing */

    const common::Grid& getGlobalGrid() const;

    /** Override default implementation Distribution::getMaxLocalSize() */

    virtual IndexType getMaxLocalSize() const;

    /** Implementation of pure method Distribution::local2global */

    virtual IndexType local2global( const IndexType localIndex ) const;

    /** This method does the local2global calculation with the grid positions. */

    void local2global( IndexType globalGridPos[], const IndexType localGridPos[] ) const;

    /** Implementation of pure method Distribution::global2local */

    virtual IndexType global2local( const IndexType globalIndex ) const;

    /** This method does the global to local calculation with the grid positions. 
     *
     *  @param[in] globalGridPos global position in the grid
     *  @param[out] localGridPos  local position in the grid (undefined if not local)
     *  @return     true if the global position is owned by this processor
     */

    bool global2local( IndexType localGridPos[], const IndexType globalGridPos[] ) const;

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

    /** Implementation of pure method Distribution::getKind for this class.  */

    virtual const char* getKind() const;

    /** static method to get kind of this class. */

    static const char* getId();

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    GridDistribution(); // disable default constructor as it has no size

    void localize();    // help routine to compute block sizes and local dimensions

    common::Grid mGlobalGrid;  // class keeps an own copy with the shape of the global grid
    common::Grid mLocalGrid;   // class keeps an own copy with the shape of the local grid

    common::Grid mProcGrid;    // class keeps an own copy with the shape of the processor grid

    // block distribution in each dimension

    IndexType mRank[SCAI_GRID_MAX_DIMENSION];        //!< rank of this processor in each dimension
    IndexType mBlockSize[ SCAI_GRID_MAX_DIMENSION];  //!< block size of each dimension
    IndexType mLB[ SCAI_GRID_MAX_DIMENSION];
    IndexType mUB[ SCAI_GRID_MAX_DIMENSION];
};

const common::Grid& GridDistribution::getLocalGrid() const
{
    return mLocalGrid;
}

const common::Grid& GridDistribution::getGlobalGrid() const
{
    return mGlobalGrid;
}

const char* GridDistribution::getKind() const
{
    return getId();
}

const char* GridDistribution::getId()
{
    return "GRID";
}

} /* end namespace dmemo */

} /* end namespace scai */