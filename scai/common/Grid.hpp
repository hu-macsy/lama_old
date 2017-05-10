/**
 * @file common/Grid.hpp
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

#define SCAI_GRID_MAX_DIMENSION 7

/** Data structure for an n-dimensional grid.
 *
 *  A grid stands for a rectangular grid with certain extensions in each dimension.
 *
 *  This class supports also linearization of indexes (row-major)
 */
class COMMON_DLL_IMPORTEXPORT Grid
{
public:

    /** Constructor of a n-dimensional grid.
     *
     *  @param[in] nDims is the number of dimensions
     *  @param[in] sizes contains the dimensions of the grid
     */
    inline Grid( const IndexType nDims, const IndexType sizes[] );

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

    /** Get grid position from a linearized index. */

    inline void gridPos( IndexType pos[], const IndexType linearPos ) const;

    /** Compares two grids for equality. */

    inline bool operator==( const Grid& other ) const;

    /** Compares two grids for inequality. */

    inline bool operator!=( const Grid& other ) const;

protected:

    inline Grid();

    IndexType mNDims;

    IndexType mSize[ SCAI_GRID_MAX_DIMENSION];
};

/** 1-dimensional Grid */

class COMMON_DLL_IMPORTEXPORT Grid1D : public Grid
{
public:

    /** Constructor of a three-dimensional grid */

    inline Grid1D( const IndexType n1 );

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

    inline Grid2D( const IndexType n1, const IndexType n2 );

    /** More convenient call of linearPos for two-dimensional grids */

    inline IndexType linearPos( const IndexType pos1, const IndexType pos2 ) const;

    using Grid::linearPos;
};

/** 3-dimensional Grid */

class COMMON_DLL_IMPORTEXPORT Grid3D : public Grid
{
public:

    /** Constructor of a three-dimensional grid */

    inline Grid3D( const IndexType n1, const IndexType n2, const IndexType n3 );

    /** More convenient call of linearPos for tree-dimensional grids */

    inline IndexType linearPos( const IndexType pos1, const IndexType pos2, const IndexType pos3 ) const;

    using Grid::linearPos;
};

/** 4-dimensional Grid */

class COMMON_DLL_IMPORTEXPORT Grid4D : public Grid
{
public:

    /** Constructor of a four-dimensional grid */

    inline Grid4D( const IndexType n1, const IndexType n2, const IndexType n3, const IndexType n4 );

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
inline COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const Grid& grid );

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline constructors                                             */
/* ------------------------------------------------------------------------------------ */

Grid::Grid()
{
    mNDims = 0;

    // initialize all sizes with a default, just in case it is reinterpreted

    for ( IndexType i = 0; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        mSize[i] = 1;
    }
}

Grid1D::Grid1D( const IndexType n1 )
{
    mNDims = 1;
    mSize[0] = n1;

    for ( IndexType i = 1; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        mSize[i] = 1;
    }
}

Grid2D::Grid2D( const IndexType n1, const IndexType n2 ) 
{
    mNDims = 2;
    mSize[0] = n1;
    mSize[1] = n2;

    for ( IndexType i = 2; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        mSize[i] = 1;
    }
}

Grid3D::Grid3D( const IndexType n1, const IndexType n2, const IndexType n3 ) 
{
    mNDims = 3;

    mSize[0] = n1;
    mSize[1] = n2;
    mSize[2] = n3;

    for ( IndexType i = 3; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        mSize[i] = 1;
    }
}

Grid4D::Grid4D( const IndexType n1, const IndexType n2, const IndexType n3, const IndexType n4 )
{
    mNDims = 4;

    mSize[0] = n1;
    mSize[1] = n2;
    mSize[2] = n3;
    mSize[3] = n4;
 
    for ( IndexType i = 4; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        mSize[i] = 1;
    }
}

Grid::Grid( const IndexType nDims, const IndexType size[] )
{
    mNDims = nDims;

    for ( IndexType i = 0; i < nDims; ++i )
    {
        mSize[i] = size[i];
    }

    // make sure that the other sizes are valid

    for ( IndexType i = nDims; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        mSize[i] = 1;
    }
}

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
    for ( IndexType idim = 0; idim < mNDims; ++idim )
    {
        sizes[idim] = mSize[idim];
    }
}

void Grid::getDistances( IndexType distances[] ) const
{
    IndexType distance = 1;

    for ( IndexType idim = mNDims; idim-- > 0; )
    {
        distances[idim] = distance;
        distance *= mSize[idim];
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

IndexType Grid::linearPos( const IndexType gridPos[] ) const
{
    IndexType pos = gridPos[0];

    for ( IndexType i = 1; i < mNDims; ++i )
    {  
        pos = pos * mSize[i] + gridPos[i];
    }

    return pos;
}

IndexType Grid1D::linearPos( const IndexType pos1 ) const
{
    return pos1;
}

IndexType Grid2D::linearPos( const IndexType pos1, const IndexType pos2 ) const
{
    return pos1 * mSize[1] + pos2;
}

IndexType Grid3D::linearPos( const IndexType pos1, const IndexType pos2, const IndexType pos3 ) const
{
    // 3-dimensional grid has always 3 dimensions, no assert required
  
    return ( pos1 * mSize[1] + pos2 ) * mSize[2] + pos3;
}

IndexType Grid4D::linearPos( const IndexType pos1, const IndexType pos2, const IndexType pos3, const IndexType pos4 ) const
{
    // 4-dimensional grid has always 3 dimensions, no assert required

    return ( ( pos1 * mSize[1] + pos2 ) * mSize[2] + pos3 ) * mSize[3] + pos4;
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

bool Grid::operator!=( const Grid& other ) const
{
    return ! operator==( other );
}

const IndexType* Grid::sizes() const
{
    return mSize;
}

IndexType Grid::size( IndexType dim ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( dim, SCAI_GRID_MAX_DIMENSION, "illegal dim" )

    return mSize[dim];
}

IndexType Grid::size() const
{
    if ( mNDims == 1 )
    {
        return mSize[0];
    }
    else
    {
        IndexType s = mSize[0];

        for ( IndexType i = 1; i < mNDims; ++i )
        {
            s *= mSize[i];
        }

        return s;
    }
}

std::ostream& operator<<( std::ostream& stream, const Grid& grid )
{
    stream << "Grid" << grid.nDims() << "D( " << grid.size(0) ;

    for ( IndexType i = 1; i < grid.nDims(); ++i )
    {
        stream << " x " << grid.size(i);
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

} /* end namespace common */

} /* end namespace scai */
