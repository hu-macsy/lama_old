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

namespace scai
{

namespace common
{

/** Maximal number of supported dimensions for a grid. 
 *
 *  Might be increased as required. 
 */

#define SCAI_GRID_MAX_DIMENSION 4

/** Data structure for an n-dimensional grid.
 *
 *  A grid stands for a rectangular grid with certain extensions in each dimension.
 *
 *  This class supports also linearization of indexes (row-major)
 */
class COMMON_DLL_IMPORTEXPORT Grid
{
public:

    /** Constructor of a one-dimensional grid. */

    Grid( const IndexType n1 );

    /** Constructor of a two-dimensional grid. */

    Grid( const IndexType n1, const IndexType n2 );

    /** Constructor of a three-dimensional grid. */

    Grid( const IndexType n1, const IndexType n2, const IndexType n3 );

    /** Constructor of a n-dimensioal grid.
     *
     *  @param[in] ndims is the number of dimensions
     *  @param[in] sizes contains the dimensions of the grid
     */
    Grid( const IndexType ndims, const IndexType sizes[] );

    /** Set explicitly the size in a dim */

    void setSize( IndexType dim, IndexType size );

    /** Return number of dimensions for the grid. */

    IndexType ndims() const;

    /** Return the total number of grid points. */

    IndexType size() const;

    /** Size of the grid in a certain dimension */

    IndexType size( const IndexType dim ) const;

    /** Predicate to check for a correct position in the grid, all position indexes must be valid */

    bool validPos( const IndexType gridPos[] ) const;

    /** Linearize grid position */

    IndexType linearPos( const IndexType gridPos[] ) const;

    /** More convenient call of linearPos on two-dimensional grids */

    IndexType linearPos( const IndexType pos1, const IndexType pos2 ) const;

    /** More convenient call of linearPos for tree-dimensional grids */

    IndexType linearPos( const IndexType pos1, const IndexType pos2, const IndexType pos3 ) const;

    /** Get grid position from a linearized index. */

    void gridPos( IndexType pos[], const IndexType linearPos ) const;

    /** Compares two grids for equality. */

    bool operator==( const Grid& other ) const;

private:

    Grid();

    IndexType mNDims;

    IndexType mSize[ SCAI_GRID_MAX_DIMENSION];
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
/*   Implementations of inline constructors                                             */
/* ------------------------------------------------------------------------------------ */

Grid::Grid( const IndexType n1 )
{
    mNDims = 1;
    mSize[0] = n1;
}

Grid::Grid( const IndexType n1, const IndexType n2 )
{
    mNDims = 2;
    mSize[0] = n1;
    mSize[1] = n2;
}

Grid::Grid( const IndexType n1, const IndexType n2, const IndexType n3 )
{
    mNDims = 3;

    mSize[0] = n1;
    mSize[1] = n2;
    mSize[2] = n3;
}

Grid::Grid( const IndexType nDims, const IndexType size[] )
{
    mNDims = nDims;

    for ( IndexType i = 0; i < nDims; ++i )
    {
        mSize[i] = size[i];
    }
}

/* ------------------------------------------------------------------------------------ */
/*   Implementations of inline methods                                                  */
/* ------------------------------------------------------------------------------------ */

IndexType Grid::ndims() const
{
    return mNDims;
}

void Grid::setSize( IndexType dim, IndexType size )
{
    mSize[dim] = size;
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

IndexType Grid::linearPos( const IndexType gridPos[] ) const
{
    if ( mNDims == 1 )
    {
        return gridPos[0];
    }
    else 
    {
        IndexType pos = gridPos[0];

        for ( IndexType i = 1; i < mNDims; ++i )
        {  
            pos = pos * mSize[i] + gridPos[i];
        }

        return pos;
    }
}

IndexType Grid::linearPos( const IndexType pos0, const IndexType pos1 ) const

{
    SCAI_ASSERT_EQ_DEBUG( 2, mNDims, "linearPos with 2 values, but grid has " << mNDims << " dimensions" )

    return pos0 * mSize[1] + pos1;
}

IndexType Grid::linearPos( const IndexType pos0, const IndexType pos1, const IndexType pos2 ) const
{
    SCAI_ASSERT_EQ_DEBUG( 3, mNDims, "linearPos with 3 values, but grid has " << mNDims << " dimensions" )

    return ( pos0 * mSize[1] + pos1 ) * mSize[2] + pos2;
}

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
    if ( other.ndims() != mNDims ) 
    {
        return false;
    }

    // same number of dimensions, now we can compare the sizes

    for ( IndexType i = 0; i < mNDims; ++i )
    {
        if ( other.size( i ) != mSize[i] )
        {
            return false;
        }
    }

    return true;
}

IndexType Grid::size( IndexType dim ) const
{
    if ( dim < mNDims )
    {
        return mSize[dim];
    }
    else
    {
        // just helpful for handling two-dimensional grids in three-dimensional context
        return 1;
    }
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
    stream << "Grid-" << grid.ndims() << "( " << grid.size(0) ;

    for ( IndexType i = 1; i < grid.ndims(); ++i )
    {
        stream << " x " << grid.size(i);
    }

    stream << " )";

    return stream;
}

} /* end namespace common */

} /* end namespace scai */
