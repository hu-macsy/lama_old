/**
 * @file common/Grid.cpp
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
 * @brief Implementation of methods of class Grid
 * @author Thomas Brandes
 * @date 30.01.2017
 */
 
#include <scai/common/Grid.hpp>

namespace scai
{

namespace common
{

/* ------------------------------------------------------------------------------------ */
/*   Implementations of constructors                                                    */
/* ------------------------------------------------------------------------------------ */

Grid::Grid()
{
    init( 0 );
}

Grid::Grid( const IndexType nDims )
{
    init( nDims );
}

Grid1D::Grid1D( const IndexType n1 )
{
    init ( 1 );

    mSize[0] = n1;
}

Grid2D::Grid2D( const IndexType n1, const IndexType n2 ) 
{
    init( 2 );

    mSize[0] = n1;
    mSize[1] = n2;
}

Grid3D::Grid3D( const IndexType n1, const IndexType n2, const IndexType n3 ) 
{
    init( 3 );

    mSize[0] = n1;
    mSize[1] = n2;
    mSize[2] = n3;
}

Grid4D::Grid4D( const IndexType n1, const IndexType n2, const IndexType n3, const IndexType n4 )
{
    init( 4 );

    mSize[0] = n1;
    mSize[1] = n2;
    mSize[2] = n3;
    mSize[3] = n4;
}

Grid::Grid( const IndexType nDims, const IndexType size[] )
{
    init( nDims );

    for ( IndexType i = 0; i < nDims; ++i )
    {
        mSize[i] = size[i];
    }
}

/* ------------------------------------------------------------------------------------ */
/*   Implementations of methods                                                         */
/* ------------------------------------------------------------------------------------ */

void Grid::init( const IndexType nDims )
{
    mNDims = nDims;

    // for convenience make sure that all dims are set
    // this allows to cast a grid to a more dimensional grid

    for ( IndexType iDim = 0; iDim < SCAI_GRID_MAX_DIMENSION; iDim++ )
    {
        mSize[iDim] = 1;
        mBorder[iDim][0] = BorderType::ABSORBING;
        mBorder[iDim][1] = BorderType::ABSORBING;
    }
}

/* ------------------------------------------------------------------------------------ */

/*  Help routine to get the position in the last dimension */

static void inline getDimPos( IndexType& dimpos, IndexType& pos, const IndexType dimSize )
{
    IndexType k = pos / dimSize;
    dimpos      = pos - k * dimSize;  // modulo part
    pos         = k;
}

void Grid::gridPos( IndexType pos[], const IndexType gridPos ) const
{
    if ( mNDims == 1 )
    {
        pos[0] = gridPos;
    }
    else
    {
        IndexType p = gridPos;

        for ( IndexType i = mNDims - 1; i > 0; --i )
        {
            getDimPos( pos[i], p, mSize[i] );

            // p is now linear index in Grid with one dim less
        }
    
        pos[0] = p;
    }
}

bool Grid::operator==( const Grid& other ) const
{
    if ( mNDims != other.mNDims )
    {
        return false;
    }

    for ( IndexType i = 0; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        if ( i >= mNDims ) 
        {
            break;
        }

        if ( mSize[i] != other.mSize[i] ) 
        {
            return false;
        }
    }

    return true;
}

static inline const char* codeBorder( BorderType border, bool leftFlag )
{
    switch ( border )
    {
        case BorderType::PERIODIC:
            if ( leftFlag ) 
                return "P ";
            else
                return " P";
        default:
            return "";
    }
}

std::ostream& operator<<( std::ostream& stream, const Grid& grid )
{
    const BorderType* borders = grid.borders();

    stream << "Grid" << grid.nDims() << "D( "
           << codeBorder( borders[0], true ) 
           << grid.size( 0 ) 
           << codeBorder( borders[1], false );

    for ( IndexType i = 1; i < grid.nDims(); ++i )
    {
        stream << " x "
               << codeBorder( borders[ 2 * i ], true ) 
               << grid.size(i) 
               << codeBorder( borders[ 2 * i + 1 ], false );
    }

    // print size values not equal 1, might be helpful for debugging

    for ( IndexType i = grid.nDims(); i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        if ( grid.size( i ) != 1 )
        {
            stream << " x " << i << " : " << grid.size(i);
        }
    }

    stream << " )";

    return stream;
}

static bool inline getBorderPosL( IndexType& pos, const IndexType offset, const IndexType size, const BorderType border )
{
    bool valid = true;

    if ( pos >= offset )
    {
        pos = pos - offset;  // is a valid pos
    }
    else if ( border == BorderType::ABSORBING )
    {
        valid = false;
    }
    else if ( border == BorderType::PERIODIC )
    {
        pos = ( pos + size ) - offset;
    }

    return valid;
}

/*  Alternative implementation using previous routine:

    pos = size - 1 - pos;
    bool valid = getBorderPosL( pos, offset, size, border );
    pos = size - 1 - pos;
    return valid;
*/

static bool inline getBorderPosR( IndexType& pos, const IndexType offset, const IndexType size, const BorderType border )
{
    bool valid = true;

    if ( pos + offset < size )
    {
        pos += offset;  // is a valid pos
    }
    else if ( border == BorderType::ABSORBING )
    {
        valid = false;
    }
    else if ( border == BorderType::PERIODIC )
    {
        pos = ( pos + offset ) - size;
    }
    return valid;
}

bool Grid::getOffsetPos( IndexType pos[], 
                         const int offsets[],
                         const IndexType sizes[], 
                         const BorderType borders[], 
                         const IndexType nDims )
{
    bool valid = true;

    for ( IndexType iDim = 0; iDim < nDims; ++iDim )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos[iDim], sizes[iDim], "pos out of grid @ dim = " << iDim )

        if ( offsets[iDim] < 0 )
        {
            valid = getBorderPosL( pos[iDim], static_cast<IndexType>( -offsets[iDim] ), sizes[iDim],  borders[2 * iDim] );
        }
        else if ( offsets[iDim] > 0 )
        {
            valid = getBorderPosR( pos[iDim], static_cast<IndexType>( offsets[iDim] ), sizes[iDim], borders[ 2 * iDim + 1] );
        }

        if ( !valid )
        {
            break;
        }

        SCAI_ASSERT_VALID_INDEX_DEBUG( pos[iDim], sizes[iDim], "new pos out of grid @ dim = " << iDim )
    }
 
    return valid;
}

bool Grid::getOffsetPos( IndexType pos[],
                         const int offsets[],
                         const IndexType sizes[],
                         const IndexType borders[],
                         const IndexType nDims )
{
    bool valid = true;

    for ( IndexType iDim = 0; iDim < nDims; ++iDim )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos[iDim], sizes[iDim], "pos out of grid @ dim = " << iDim )

        if ( offsets[iDim] < 0 )
        {
            valid = getBorderPosL( pos[iDim], static_cast<IndexType>( -offsets[iDim] ), sizes[iDim], static_cast<BorderType>( borders[2 * iDim] ) );
        }
        else if ( offsets[iDim] > 0 )
        {
            valid = getBorderPosR( pos[iDim], static_cast<IndexType>( offsets[iDim] ), sizes[iDim], static_cast<BorderType>( borders[ 2 * iDim + 1] ) );
        }

        if ( !valid )
        {
            break;
        }

        SCAI_ASSERT_VALID_INDEX_DEBUG( pos[iDim], sizes[iDim], "new pos out of grid @ dim = " << iDim )
    }

    return valid;
}

bool Grid::getOffsetPos( IndexType pos[], const int offsets[] ) const
{
    // call the static method with my member variables

    return getOffsetPos( pos, offsets, mSize, mBorder[0], mNDims );
}


} /* end namespace common */

} /* end namespace scai */
