/**
 * @file GridDistribution.cpp
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
 * @brief Implementation of methods for grid distribution class.
 * @author Thomas Brandes
 * @date 18.03.2011
 */

// hpp
#include <scai/dmemo/GridDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

// std
#include <fstream>

namespace scai
{

using common::Grid;
using common::Grid1D;
using common::Grid2D;
using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GridDistribution::logger, "Distribution.GridDistribution" )

/* ---------------------------------------------------------------------- */

GridDistribution::GridDistribution( const Grid& globalGrid, const CommunicatorPtr communicator, const Grid& procGrid ) :

    Distribution( globalGrid.size(), communicator ),

    mGlobalGrid( globalGrid ),
    mLocalGrid( globalGrid ),
    mProcGrid( procGrid )
{
    SCAI_ASSERT_EQ_ERROR( procGrid.nDims(), globalGrid.nDims(), "Global grid and processor grid must have same number of dims" )
    SCAI_ASSERT_LE_ERROR( procGrid.size(), communicator->getSize(), "processor array too big" )

    localize();
}

/* ---------------------------------------------------------------------- */

GridDistribution::GridDistribution( const Grid& globalGrid, const CommunicatorPtr communicator ) :

    Distribution( globalGrid.size(), communicator ),

    mGlobalGrid( globalGrid ),
    mLocalGrid( globalGrid ),
    mProcGrid( globalGrid.nDims() )
{
    // Build a correct processor grid by the communicator
    // Note: border type of the processor grid does not matter

    IndexType nDims = globalGrid.nDims();

    if ( nDims == 1 )
    {
        mProcGrid.setSize( 0, communicator->getSize() );
    }
    else if ( nDims == 2 )
    {
        PartitionId procGrid[2];

        double size = static_cast<double>( globalGrid.size( 0 ) + globalGrid.size( 1 ) );

        double w1 = globalGrid.size( 0 ) / size;
        double w2 = globalGrid.size( 1 ) / size;

        communicator->factorize2( procGrid, w1, w2 );

        mProcGrid.setSize( 0, procGrid[0] );
        mProcGrid.setSize( 1, procGrid[1] );
    }
    else
    {
        PartitionId procGrid[3];

        double size = static_cast<double>( globalGrid.size( 0 ) + globalGrid.size( 1 ) + globalGrid.size( 2 ) );

        double w1 = globalGrid.size( 0 ) / size;
        double w2 = globalGrid.size( 1 ) / size;
        double w3 = globalGrid.size( 2 ) / size;

        communicator->factorize3( procGrid, w1, w2, w3 );

        mProcGrid.setSize( 0, procGrid[0] );
        mProcGrid.setSize( 1, procGrid[1] );
        mProcGrid.setSize( 2, procGrid[2] );

        // set the remaining dimensions to 1

        for ( IndexType i = 3; i < nDims; ++i )
        {
            mProcGrid.setSize( i, 1 );
        }
    }

    localize(); // compute mBlockSize, mLB, mUB, mRank, m
}

/* ---------------------------------------------------------------------- */

void GridDistribution::setLocalGrid( Grid& localGrid, const PartitionId rank ) const
{
    localGrid = mGlobalGrid;

    if ( rank < mProcGrid.size() )
    {
        IndexType procPos[SCAI_GRID_MAX_DIMENSION];

        mProcGrid.gridPos( procPos, rank );

        for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
        {
            IndexType lb, ub;

            BlockDistribution::getLocalRange( lb, ub, mGlobalGrid.size( idim ), procPos[idim], mProcGrid.size( idim ) );

            if ( lb <= ub )
            {
                localGrid.setSize( idim, ub - lb );
            }
            else
            {
                localGrid.setSize( idim, 0 );
            }
        }
    }
    else
    {
        for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
        {
            localGrid.setSize( idim, 0 );
        }
    }
}

/* ---------------------------------------------------------------------- */

void GridDistribution::localize()
{
    // compute the block sizes for each grid dimension

    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        IndexType nProcs = mProcGrid.size( idim );  // number of procs for this dimension

        mBlockSize[idim] = ( mGlobalGrid.size( idim ) + nProcs - 1 ) / nProcs;
    }

    PartitionId pRank = mCommunicator->getRank();

    // determine my grid position in procGrid

    if ( pRank < mProcGrid.size() )
    {
        mProcGrid.gridPos( mRank, pRank );

        // compute block sizes and local ranges

        for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
        {
            if ( mProcGrid.size( idim ) == 1 )
            {
                // this dimension is not distributed at all

                mLB[idim] = 0;
                mUB[idim] = mGlobalGrid.size( idim );

                // note: border types of the local grid do not change in this dimension
            }
            else
            {
                BlockDistribution::getLocalRange( mLB[idim], mUB[idim], mGlobalGrid.size( idim ), mRank[idim], mProcGrid.size( idim ) );

                // Be careful: size = 5, np = 4 ->  0:2, 2:4, 4:5, 6:5

                if ( mLB[idim] <= mUB[idim] )
                {
                    mLocalGrid.setSize( idim, mUB[idim] - mLB[idim] );

                }
                else
                {
                    mLocalGrid.setSize( idim, 0 );
                }

                // only reflecting dimensions are kept at the borders, all other will be ABSORBIC

                Grid::BorderType left  = Grid::BORDER_ABSORBING;
                Grid::BorderType right = Grid::BORDER_ABSORBING;

                mLocalGrid.setBorderType ( idim, left, right );
            }
        }
    }
    else
    {
        // set some convenient default values for rank, lb, ub, and mLocalGrid

        for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
        {
            mRank[ idim ] = invalidIndex;
            mLB[ idim ]   = 0;
            mUB[ idim ]   = 0;

            mLocalGrid.setSize( idim, 0 );
        }
    }
}

/* ---------------------------------------------------------------------- */

GridDistribution::~GridDistribution()
{
    SCAI_LOG_DEBUG( logger, "~GridDistribution" )
}

/* ---------------------------------------------------------------------- */

bool GridDistribution::isLocal( const IndexType globalIndex ) const
{
    //  check that this processor belongs to processor array

    if ( mCommunicator->getRank() >= mProcGrid.size() )
    {
        return false;
    }

    IndexType globalGridPos[ SCAI_GRID_MAX_DIMENSION ];

    mGlobalGrid.gridPos( globalGridPos, globalIndex );

    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        if ( globalGridPos[idim] < mLB[idim] )
        {
            return false;
        }

        if ( globalGridPos[idim] >= mUB[idim] )
        {
            return false;
        }
    }

    return true;
}

/* ---------------------------------------------------------------------- */

PartitionId GridDistribution::findOwner( const IndexType globalIndex ) const
{
    // ownership can be computed locally without communication

    IndexType globalGridPos[ SCAI_GRID_MAX_DIMENSION ];

    mGlobalGrid.gridPos( globalGridPos, globalIndex );

    IndexType gridOwner[ SCAI_GRID_MAX_DIMENSION ];

    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        gridOwner[idim] = globalGridPos[idim] / mBlockSize[ idim ];
    }

    IndexType owner = mProcGrid.linearPos( gridOwner );

    return owner;
}

/* ---------------------------------------------------------------------- */

IndexType GridDistribution::getLocalSize() const
{
    //  check that this processor belongs to processor array

    return mLocalGrid.size();
}

/* ---------------------------------------------------------------------- */

IndexType GridDistribution::getMaxLocalSize() const
{
    // multiply maximal block sizes for each dimension

    IndexType maxSize = 1;

    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        maxSize *= mBlockSize[idim];
    }

    return maxSize;
}

/* ---------------------------------------------------------------------- */

IndexType GridDistribution::getBlockDistributionSize() const
{
    // only for one-dimensional grid we have a general block distribution

    bool isBlocked = true;

    for ( IndexType idim = 1; idim < mProcGrid.nDims(); ++idim )
    {
        if ( mProcGrid.size( idim ) > 1 )
        {
            isBlocked = false;
        }
    }

    if ( isBlocked )
    {
        return mLocalGrid.size();
    }
    else
    {
        return invalidIndex;
    }
}

/* ---------------------------------------------------------------------- */

void GridDistribution::local2global( IndexType globalGridPos[], const IndexType localGridPos[] ) const
{
    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        globalGridPos[idim] = mLB[idim] + localGridPos[idim];
    }
}

/* ---------------------------------------------------------------------- */

IndexType GridDistribution::local2global( const IndexType localIndex ) const
{
    IndexType localGridPos[SCAI_GRID_MAX_DIMENSION];

    mLocalGrid.gridPos( localGridPos, localIndex );

    IndexType globalGridPos[SCAI_GRID_MAX_DIMENSION];

    local2global( globalGridPos, localGridPos );

    return mGlobalGrid.linearPos( globalGridPos );
}

/* ---------------------------------------------------------------------- */

bool GridDistribution::global2local( IndexType localGridPos[], const IndexType globalGridPos[] ) const
{
    bool isLocal = true;

    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        if ( globalGridPos[idim] < mLB[idim ] || globalGridPos[idim] >= mUB[idim] )
        {
            isLocal = false;
            break;
        }

        localGridPos[idim] = globalGridPos[idim] - mLB[idim];
    }

    return isLocal;
}

/* ---------------------------------------------------------------------- */

IndexType GridDistribution::global2local( const IndexType globalIndex ) const
{
    IndexType globalGridPos[SCAI_GRID_MAX_DIMENSION];

    mGlobalGrid.gridPos( globalGridPos, globalIndex );

    IndexType localGridPos[SCAI_GRID_MAX_DIMENSION];

    bool isLocal = global2local( localGridPos, globalGridPos );

    if ( !isLocal )
    {
        SCAI_LOG_DEBUG( logger, "globalIndex = " << globalIndex << " not local, dist = " << *this )
        return invalidIndex;
    }

    return mLocalGrid.linearPos( localGridPos );
}

/* ---------------------------------------------------------------------- */

void GridDistribution::computeOwners( HArray<PartitionId>& owners, const HArray<IndexType>& indexes ) const
{
    ContextPtr ctx = Context::getHostPtr();    // currently only available @ Host

    const IndexType n = indexes.size();

    ReadAccess<IndexType> rIndexes( indexes, ctx );
    WriteOnlyAccess<PartitionId> wOwners( owners, ctx, n );

    // ToDo: call a kernel and allow arbitrary context

    for ( IndexType i = 0; i < n; i++ )
    {
        wOwners[i] = GridDistribution::findOwner( rIndexes[i] );
    }
}

/* ---------------------------------------------------------------------- */

void GridDistribution::getOwnedIndexes( hmemo::HArray<IndexType>& myGlobalIndexes ) const
{
    const IndexType nLocal  = getLocalSize();

    SCAI_LOG_INFO( logger, getCommunicator() << ": getOwnedIndexes, have " << nLocal << " of " << mGlobalSize )

    WriteOnlyAccess<IndexType> wGlobalIndexes( myGlobalIndexes, nLocal );

    for ( IndexType i = 0; i < nLocal; ++i )
    {
        wGlobalIndexes[i] = local2global( i );
    }
}

/* ---------------------------------------------------------------------- */

bool GridDistribution::hasAnyAddressing() const
{
    return true;
}

/* ---------------------------------------------------------------------- */

void GridDistribution::enableAnyAddressing() const
{
    // nothing to do, we have closed formulas
}

IndexType GridDistribution::getAnyLocalSize( const PartitionId rank ) const
{
    Grid localGrid( mProcGrid );

    setLocalGrid( localGrid, rank );

    return localGrid.size();
}

PartitionId GridDistribution::getAnyOwner( const IndexType globalIndex ) const
{
    // same as findOwner

    IndexType globalGridPos[ SCAI_GRID_MAX_DIMENSION ];

    mGlobalGrid.gridPos( globalGridPos, globalIndex );

    IndexType gridOwner[ SCAI_GRID_MAX_DIMENSION ];

    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        gridOwner[idim] = globalGridPos[idim] / mBlockSize[ idim ];
    }

    IndexType owner = mProcGrid.linearPos( gridOwner );

    return owner;
}

IndexType GridDistribution::getAnyLocalIndex( const IndexType globalIndex, const PartitionId owner ) const
{
    IndexType gridPos[ SCAI_GRID_MAX_DIMENSION ];

    mGlobalGrid.gridPos( gridPos, globalIndex );

    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        gridPos[idim] = gridPos[idim] % mBlockSize[ idim ];
    }

    Grid localGrid( mGlobalGrid );

    setLocalGrid( localGrid, owner );

    IndexType localIndex = localGrid.linearPos( gridPos );

    return localIndex;
}

/* ---------------------------------------------------------------------- */

IndexType GridDistribution::getAnyGlobalIndex( const IndexType localIndex, const PartitionId owner ) const
{
    Grid localGrid( mGlobalGrid );

    setLocalGrid( localGrid, owner );

    IndexType gridPos[ SCAI_GRID_MAX_DIMENSION ];
    IndexType gridRank[ SCAI_GRID_MAX_DIMENSION ];

    localGrid.gridPos( gridPos, localIndex );
    mProcGrid.gridPos( gridRank, owner );

    for ( IndexType idim = 0; idim < mGlobalGrid.nDims(); ++idim )
    {
        gridPos[idim] = gridPos[idim] + mBlockSize[idim] * gridRank[idim];
    }

    return mGlobalGrid.linearPos( gridPos );
}

/* ---------------------------------------------------------------------- */

bool GridDistribution::isEqual( const Distribution& other ) const
{
    bool isSame = false;

    bool proven = proveEquality( isSame, other );

    if ( proven )
    {
        return isSame;
    }

    if ( other.getKind() == getKind() )
    {
        const GridDistribution& gridOther = reinterpret_cast<const GridDistribution&>( other );

        isSame = ( mGlobalGrid == gridOther.mGlobalGrid ) && ( mProcGrid == gridOther.mProcGrid );
    }

    // we know already that global size and communicator are equal

    return isSame;
}

/* ---------------------------------------------------------------------- */

void GridDistribution::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "GridDistribution( comm = " << *mCommunicator << ", grid = " << mGlobalGrid
           << " onto " << mProcGrid << ", local grid = " << mLocalGrid << " )";
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string GridDistribution::createValue()
{
    return getId();
}

DistributionPtr GridDistribution::create( const DistributionArguments arg )
{
    SCAI_LOG_INFO( logger, "create" )

    // Note: weight argument is not used here

    IndexType globalSize = arg.globalSize;

    IndexType size1 = 1;
    IndexType size2 = globalSize;

    // try a factorization in 2 arguments

    for ( IndexType i1 = 2; i1 < globalSize; i1++ )
    {
        IndexType i2 = globalSize / i1;

        if ( i2 < i1 )
        {
            break;   // makes sure that loop runs only up to sqrt( globalSize )
        }

        if ( i1* i2  != globalSize )
        {
            continue;
        }

        // i1 * i2 == size, i1 <= i2 , is a good factorization

        size1 = i1;
        size2 = i2;
    }

    if ( size1 == 1 )
    {
        return std::make_shared<GridDistribution>( Grid1D( globalSize ), arg.communicator );
    }
    else
    {
        return std::make_shared<GridDistribution>( Grid2D( size1, size2 ), arg.communicator );
    }
}

} /* end namespace dmemo */

} /* end namespace scai */
