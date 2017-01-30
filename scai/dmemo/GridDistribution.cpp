/**
 * @file GridDistribution.cpp
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
using namespace hmemo;

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( GridDistribution::logger, "Distribution.GridDistribution" )

GridDistribution::GridDistribution( const Grid& globalGrid, const CommunicatorPtr communicator, const Grid& procGrid ) :
 
    Distribution( globalGrid.size(), communicator ),

    mGlobalGrid( globalGrid ),
    mLocalGrid( globalGrid ),
    mProcGrid( procGrid )
{
    SCAI_ASSERT_EQ_ERROR( procGrid.ndims(), globalGrid.ndims(), "Global grid and processor grid must have same number of dims" )
    SCAI_ASSERT_LE_ERROR( procGrid.size(), communicator->getSize(), "processor array too big" )

    // determine my rank

    PartitionId pRank = communicator->getRank();

    for ( int idim = 0; idim < globalGrid.ndims(); ++idim )
    {
        mBlockSize[idim] = ( globalGrid.size( idim ) + procGrid.size( idim ) - 1 ) / procGrid.size( idim );
    }

    if ( pRank < procGrid.size() )
    {
        mProcGrid.gridPos( mRank, pRank );
    
        // compute block sizes and local ranges

        for ( int idim = 0; idim < globalGrid.ndims(); ++idim )
        {
            BlockDistribution::getLocalRange( mLB[idim], mUB[idim], globalGrid.size( idim ), mRank[idim], procGrid.size( idim ) );

            if ( mLB[idim] <= mUB[idim] )
            {
                mLocalGrid.setSize( idim, mUB[idim] - mLB[idim] );
            }
            else
            {
                mLocalGrid.setSize( idim, 0 );
            }
        }
    }
    else
    {
         // set some convenient default values for rank, lb, ub, and mLocalGrid

        for ( int idim = 0; idim < globalGrid.ndims(); ++idim )
        {
            mRank[ idim ] = nIndex;
            mLB[ idim ]   = 0;
            mUB[ idim ]   = 0;

            mLocalGrid.setSize( idim, 0 );
        }
    }

    // be careful: size of local grid might be 0
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
 
    for ( IndexType idim = 0; idim < mGlobalGrid.ndims(); ++idim )
    {
        if ( globalGridPos[idim] < mLB[idim] ) return false;
        if ( globalGridPos[idim] >= mUB[idim] ) return false;
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

    for ( IndexType idim = 0; idim < mGlobalGrid.ndims(); ++idim )
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
    
    for ( IndexType idim = 0; idim < mGlobalGrid.ndims(); ++idim )
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

    for ( IndexType idim = 1; idim < mProcGrid.ndims(); ++idim )
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
        return nIndex;
    }
}

/* ---------------------------------------------------------------------- */

IndexType GridDistribution::local2global( const IndexType localIndex ) const
{
    IndexType localGridPos[SCAI_GRID_MAX_DIMENSION];

    mLocalGrid.gridPos( localGridPos, localIndex );

    IndexType globalGridPos[SCAI_GRID_MAX_DIMENSION];

    for ( IndexType idim = 0; idim < mGlobalGrid.ndims(); ++idim )
    {
        globalGridPos[idim] = mLB[idim] + localGridPos[idim];
    }

    return mGlobalGrid.linearPos( globalGridPos );
}

/* ---------------------------------------------------------------------- */

IndexType GridDistribution::global2local( const IndexType globalIndex ) const
{
    IndexType globalGridPos[SCAI_GRID_MAX_DIMENSION];

    mGlobalGrid.gridPos( globalGridPos, globalIndex );

    IndexType localGridPos[SCAI_GRID_MAX_DIMENSION];

    bool nonLocal = false;

    for ( IndexType idim = 0; idim < mGlobalGrid.ndims(); ++idim )
    {   
        if ( globalGridPos[idim] < mLB[idim ] || globalGridPos[idim] >= mUB[idim] )
        {
            nonLocal = true;
            break;
        }

        localGridPos[idim] = globalGridPos[idim] - mLB[idim];
    }

    if ( nonLocal )
    {
        return nIndex;
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

    stream << "GridDistribution( comm = " << *mCommunicator << ", grid )";
}

/* ---------------------------------------------------------------------------------*
 *   static create methods ( required for registration in distribution factory )    *
 * ---------------------------------------------------------------------------------*/

std::string GridDistribution::createValue()
{
    return getId();
}

Distribution* GridDistribution::create( const DistributionArguments arg )
{
    SCAI_LOG_INFO( logger, "create" )

    // Note: weight argument is not used here

    return new GridDistribution( Grid( arg.globalSize ), arg.communicator, Grid( arg.communicator->getSize() ) );
}

} /* end namespace dmemo */

} /* end namespace scai */
