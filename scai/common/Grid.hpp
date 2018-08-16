/**
 * @file common/Grid.hpp
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
 * @brief Definition of a class for multidimensional grids
 * @author Thomas Brandes
 * @date 30.01.2017
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/Utils.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>

// std
#include <stdint.h>
#include <iostream>

namespace scai
{

namespace common
{

/** Maximal number of supported dimensions for a grid. 
 *
 *  Might be increased as required. 
 */

#define SCAI_GRID_MAX_DIMENSION ( (IndexType) 4 )

/** Data structure for an n-dimensional grid.
 *
 *  A grid stands for a rectangular grid with certain extensions in each dimension.
 *
 *  This class supports also linearization of indexes (row-major)
 */
class COMMON_DLL_IMPORTEXPORT Grid
{
public:

    /** BorderType describes for each side of the grid which value is taken
     *  if the corresponding neighbored pos is not available. 
     */
    typedef enum
    {
        BORDER_ABSORBING,    //!   factor 0 if neighbored pos is not available
        BORDER_PERIODIC,     //!   take it from other side
    } BorderType;

    /** Default constructor provides an empty grid with zero dimensions.
     *
     */
    Grid();

    /** Constructor of a n-dimensional grid.
     *
     *  @param[in] nDims is the number of dimensions
     *
     *  All sizes are set to 1.
     */
    Grid( const IndexType nDims );

    /** Constructor of a n-dimensional grid.
     *
     *  @param[in] nDims is the number of dimensions
     *  @param[in] sizes contains the dimensions of the grid
     */
    Grid( const IndexType nDims, const IndexType sizes[] );

    /** Set explicitly the size in a dim */

    inline void setSize( IndexType dim, IndexType size );

    /** Return number of dimensions for the grid. */

    inline IndexType nDims() const;

    /** Return the total number of grid points. */

    inline IndexType size() const;

    /** Size of the grid in a certain dimension */

    inline IndexType size( const IndexType dim ) const;

    /** Size of the grid in all dimensions as pointer to array */

    inline const IndexType* sizes() const;

    /** Borders of the grid in all dimensions */
  
    inline const BorderType* borders() const;

    /** Return an array with sizes of all grid dimensions */

    inline void getSizes( IndexType sizes[] ) const;

    /** Return an array with distances between grid elements for linear positions.
     *
     *  linearPos( gridPos[] ) = gridPos[0] * distances[0] + gridPos[1] * distances[1] + ... 
     */

    inline void getDistances( IndexType distances[] ) const;

    /** Predicate to check for a correct position in the grid, all position indexes must be valid */

    inline bool validPos( const IndexType gridPos[] ) const;

    /** Linearize grid position */

    inline IndexType linearPos( const IndexType gridPos[] ) const;

    /** Get grid position from a linearized index.
     *
     *  Due to divide/modulo operation this routine is rather expensive.
     *  Grids should be traversed by multidimensional indexes if possible.
     */

    void gridPos( IndexType pos[], const IndexType linearPos ) const;

    /** Compares two grids for equality. */

    bool operator==( const Grid& other ) const;

    /** Compares two grids for inequality. */

    inline bool operator!=( const Grid& other ) const;

    /** This method allows to set a border type for one dimension at both sides  */

    inline void setBorderType( const IndexType dim, const BorderType type );

    /** This method allows to set individual border types for the grid. */

    inline void setBorderType( const IndexType dim, const BorderType left, const BorderType right );

    /** This method determines a neighbored position of a point in the grid.
     *
     *  @param[in,out] pos contains the position as input and will contain the new position
     *  @param[in]     offsets specifies the offsets in each dim, can be positive or negative
     *  @param[in]     sizes   sizes for each dim needed for checking boundaries
     *  @param[in]     borders specify the boundary types of the grid, size is 2 times nDims
     *  @param[in]     nDims   number of dimensions for the position
     */
    static bool getOffsetPos( 
        IndexType pos[], 
        const int offsets[],
        const IndexType sizes[], 
        const BorderType borders[], 
        const IndexType nDims );

    /** This method determines a neighbored position of a point in the grid.
     *
     *  @param[in,out] pos contains the position as input and will contain the new position
     *  @param[in]     offsets specifies the offsets in each dim, can be positive or negative
     *
     *  At the boundaries of the grid the grid border types play an important role how the
     *  corresponding dimension is determined.
     */
    bool getOffsetPos( IndexType pos[], const int offsets[] ) const;

protected:

    inline void init( const IndexType nDims );

    IndexType mNDims;

    IndexType mSize[ SCAI_GRID_MAX_DIMENSION];

    BorderType mBorder[SCAI_GRID_MAX_DIMENSION][2];
};

/** 1-dimensional Grid */

class COMMON_DLL_IMPORTEXPORT Grid1D : public Grid
{
public:

    /** Constructor of a three-dimensional grid */

    Grid1D( const IndexType n1 );

    /** More convenient call of linearPos for one-dimensional grid.
     *  This routine is not really necessary as it is just the identity.
     */

    inline IndexType linearPos( const IndexType pos ) const;

    using Grid::linearPos;
};

/** 2-dimensional Grid */

class COMMON_DLL_IMPORTEXPORT Grid2D : public Grid
{
public:

    /** Constructor of a three-dimensional grid */

    Grid2D( const IndexType n1, const IndexType n2 );

    /** More convenient call of linearPos for two-dimensional grids */

    inline IndexType linearPos( const IndexType pos1, const IndexType pos2 ) const;

    using Grid::linearPos;
};

/** 3-dimensional Grid */

class COMMON_DLL_IMPORTEXPORT Grid3D : public Grid
{
public:

    /** Constructor of a three-dimensional grid */

    Grid3D( const IndexType n1, const IndexType n2, const IndexType n3 );

    /** More convenient call of linearPos for tree-dimensional grids */

    inline IndexType linearPos( const IndexType pos1, const IndexType pos2, const IndexType pos3 ) const;

    using Grid::linearPos;
};

/** 4-dimensional Grid */

class COMMON_DLL_IMPORTEXPORT Grid4D : public Grid
{
public:

    /** Constructor of a four-dimensional grid */

    Grid4D( const IndexType n1, const IndexType n2, const IndexType n3, const IndexType n4 );

    /** More convenient call of linearPos for tree-dimensional grids */

    inline IndexType linearPos( const IndexType pos1, const IndexType pos2, const IndexType pos3, const IndexType pos4 ) const;

    using Grid::linearPos;
};

/*
 * Output of Grid in stream by writing its extensions
 *
 * \code
 *    cout << grid;   // e.g. prints Grid-3( 2 x 3 x 2 )
 * \endcode
 */
COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const Grid& grid );

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline methods                                                  */
/* ------------------------------------------------------------------------------------ */

IndexType Grid::nDims() const
{
    return mNDims;
}

void Grid::setSize( IndexType dim, IndexType size )
{
    mSize[dim] = size;
}


void Grid::getSizes( IndexType sizes[] ) const
{
    for ( IndexType iDim = 0; iDim < mNDims; ++iDim )
    {
        sizes[iDim] = mSize[iDim];
    }
}

void Grid::getDistances( IndexType distances[] ) const
{
    IndexType distance = 1;

    for ( IndexType iDim = mNDims; iDim-- > 0; )
    {
        distances[iDim] = distance;
        distance *= mSize[iDim];
    }
}

bool Grid::validPos( const IndexType gridPos[] ) const
{
    bool isValid = true;

    for ( IndexType i = 0; i < mNDims; ++i )
    {
        isValid = Utils::validIndex( gridPos[i], mSize[i] );

        if ( !isValid )
        {
            break;
        }
    }

    return isValid;
}

/* ------------------------------------------------------------------------------------ */

void Grid::setBorderType( const IndexType dim, const BorderType type )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( dim, mNDims, "illegal dim used" )

    mBorder[dim][0] = type;
    mBorder[dim][1] = type;
}

void Grid::setBorderType( const IndexType dim, const BorderType left, const BorderType right )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( dim, mNDims, "illegal dim used" )

    mBorder[dim][0] = left;
    mBorder[dim][1] = right;
}

/* ------------------------------------------------------------------------------------ */

IndexType Grid::linearPos( const IndexType gridPos[] ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( gridPos[0], mSize[0], "grid index out of range" )

    IndexType pos = gridPos[0];

    for ( IndexType i = 1; i < mNDims; ++i )
    {  
        pos = pos * mSize[i] + gridPos[i];
        SCAI_ASSERT_VALID_INDEX_DEBUG( gridPos[i], mSize[i], "grid index out of range" )
    }

    return pos;
}

IndexType Grid1D::linearPos( const IndexType pos0 ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos0, mSize[0], "grid index out of range" )

    return pos0;
}

IndexType Grid2D::linearPos( const IndexType pos0, const IndexType pos1 ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos0, mSize[0], "grid index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos1, mSize[1], "grid index out of range" )

    return pos0 * mSize[1] + pos1;
}

IndexType Grid3D::linearPos( const IndexType pos0, const IndexType pos1, const IndexType pos2 ) const
{
    // 3-dimensional grid has always 3 dimensions, no assert required
  
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos0, mSize[0], "grid index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos1, mSize[1], "grid index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos2, mSize[2], "grid index out of range" )

    return ( pos0 * mSize[1] + pos1 ) * mSize[2] + pos2;
}

IndexType Grid4D::linearPos( const IndexType pos0, const IndexType pos1, const IndexType pos2, const IndexType pos3 ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos0, mSize[0], "grid index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos1, mSize[1], "grid index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos2, mSize[2], "grid index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos3, mSize[3], "grid index out of range" )

    // 4-dimensional grid has always 4 dimensions, no assert required

    return ( ( pos0 * mSize[1] + pos1 ) * mSize[2] + pos2 ) * mSize[3] + pos3;
}

/* ------------------------------------------------------------------------------------ */

bool Grid::operator!=( const Grid& other ) const
{
    return ! operator==( other );
}

const IndexType* Grid::sizes() const
{
    return mSize;
}

const Grid::BorderType* Grid::borders() const
{
    return mBorder[0];
}

IndexType Grid::size( IndexType dim ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( dim, SCAI_GRID_MAX_DIMENSION, "illegal dim" )
 
    // ndims <= dim < SCAI_GRID_MAX_DIMENSION is allowed, returns 1

    return mSize[dim];
}

IndexType Grid::size() const
{
    if ( mNDims == 0 )
    {
        return 0;
    }

    IndexType s = mSize[0];

    for ( IndexType i = 1; i < mNDims; ++i )
    {
        s *= mSize[i];
    }

    return s;
}

} /* end namespace common */

} /* end namespace scai */
