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
#include <scai/common/SCAITypes.hpp>

// std
#include <stdint.h>

namespace scai
{

namespace common
{

#define SCAI_GRID_MAX_DIMENSION 3

class COMMON_DLL_IMPORTEXPORT Grid
{
public:

    /** Constructor of a one-dimensional grid. */

    Grid( const IndexType n1 );

    /** Constructor of a two-dimensional grid. */

    Grid( const IndexType n1, const IndexType n2 );

    /** Constructor of a three-dimensional grid. */

    Grid( const IndexType n1, const IndexType n2, const IndexType n3 );

    /** Return the total number of grid points. */

    IndexType size() const;

    /** Size of the grid in a certain dimension */

    IndexType size( const IndexType dim ) const;

    /** Linearize grid position */

    IndexType linearPos( const IndexType gridPos[] ) const;

    void gridPos( IndexType pos[], const IndexType linearPos ) const;

private:

    Grid();

    IndexType mNDims;

    IndexType mSize[ SCAI_GRID_MAX_DIMENSION];
};

IndexType Grid::linearPos( const IndexType gridPos[] ) const
{
    if ( mNDims == 1 )
    {
        return gridPos[0];
    }
    else if ( mNDims == 2 )
    {
        return gridPos[0] * mSize[1] + gridPos[1];
    }
    else if ( mNDims == 3 )
    {
        return gridPos[0] * mSize[1] * mSize[2] + gridPos[1] * mSize[2] + gridPos[2];
    }
    else
    {
        return nIndex;
    }
}

static void inline split( IndexType& p1, IndexType& p2, const IndexType p, const IndexType s )
{
    p1 = p / s;
    p2 = p - p1 * s; 
}

void Grid::gridPos( IndexType pos[], const IndexType gridPos ) const
{
    if ( mNDims == 1 )
    {
        pos[0] = gridPos;
    }
    else if ( mNDims == 2 )
    {
        split( pos[0], pos[1], gridPos, mSize[1] );
    }
    else if ( mNDims == 3 )
    {
        IndexType k;
        split( pos[0], k, gridPos, mSize[1] * mSize[2] );
        split( pos[1], pos[2], k, mSize[2] );
    }
}

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

IndexType Grid::size( IndexType dim ) const
{
    return mSize[dim];
}

IndexType Grid::size() const
{
    if ( mNDims == 1 )
    {
        return mSize[0];
    }
    else if ( mNDims == 2 )
    {
        return mSize[0] * mSize[1];
    }
    else if ( mNDims == 3 )
    {
        return mSize[0] * mSize[1] * mSize[2];
    }
    else
    {
        return nIndex;
    }
}

} /* end namespace common */

} /* end namespace scai */
