;/**
 * @file GridTest.cpp
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
 * @brief Test routines for class Grid
 * @author Thomas Brandes
 * @date 30.01.2017
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Grid.hpp>
#include <scai/common/test/TestMacros.hpp>

using scai::IndexType;
using scai::invalidIndex;

using scai::common::Grid;
using scai::common::Grid1D;
using scai::common::Grid2D;
using scai::common::Grid3D;
using scai::common::Grid4D;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( GridTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    IndexType n1 = 10;
    Grid1D grid1( n1 );
    BOOST_CHECK_EQUAL( n1, grid1.size() );
    n1 = 2;
    IndexType n2 = 5;
    Grid2D grid2( n1, n2 );
    BOOST_CHECK_EQUAL( n1 * n2 , grid2.size() );
    n1 = 2;
    n2 = 3;
    IndexType n3 = 2;
    Grid3D grid3( n1, n2, n3 );
    BOOST_CHECK_EQUAL( n1 * n2 * n3, grid3.size() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( grid1Test )
{
    IndexType n1 = 10;

    Grid1D grid1( n1 );
    Grid1D grid2( grid1 );

    BOOST_CHECK_EQUAL( grid1.nDims(), grid2.nDims() );

    for ( IndexType i = 0; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        BOOST_CHECK_EQUAL( grid1.size( i ), grid2.size( i ) );
    }

    Grid1D grid3( 0 );
    grid3 = grid1;

    BOOST_CHECK_EQUAL( grid1.nDims(), grid3.nDims() );

    for ( IndexType i = 0; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        BOOST_CHECK_EQUAL( grid1.size( i ), grid2.size( i ) );
    }

    Grid gridX( 0, NULL );
    gridX = grid1;

    BOOST_CHECK_EQUAL( grid1.nDims(), gridX.nDims() );

    for ( IndexType i = 0; i < SCAI_GRID_MAX_DIMENSION; ++i )
    {
        BOOST_CHECK_EQUAL( gridX.size( i ), gridX.size( i ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( validTest )
{
    Grid3D grid3( 2, 3, 2 );

    // Define points with full size to avoid warning messages of overcautious compiler

    IndexType p1[] = { 1, 2, 1, invalidIndex };
    IndexType p2[] = { 1, 2, 2, invalidIndex };

    BOOST_CHECK( grid3.validPos( p1 ) );
    BOOST_CHECK( !grid3.validPos( p2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( pos2Test )
{
    const IndexType n1 = 10;
    const IndexType n2 = 15;

    Grid2D grid( n1, n2 );

    const IndexType n = grid.size();

    // Define points with full size to avoid warning messages of overcautious compiler

    IndexType p1[] = { 3, 5, invalidIndex, invalidIndex };
    IndexType p2[] = { 3, 6, invalidIndex, invalidIndex };
    IndexType p3[] = { 4, 5, invalidIndex, invalidIndex };

    BOOST_CHECK_EQUAL( grid.linearPos( p1), grid.linearPos( p1[0], p1[1] ) );

    // verify row-major ordering

    BOOST_CHECK_EQUAL( grid.linearPos( p1 ) + 1, grid.linearPos( p2 ) );
    BOOST_CHECK_EQUAL( grid.linearPos( p1 ) + n2, grid.linearPos( p3 ) );

    for ( IndexType i = 0; i < n; ++i )
    {
        IndexType pos[2];
        grid.gridPos( pos, i );
        BOOST_CHECK_EQUAL( i, grid.linearPos( pos ) );
    }

    for ( IndexType i1 = 0; i1 < n1; ++i1 )
    {
        for ( IndexType i2 = 0; i2 < n2; ++i2 )
        {
            IndexType linearPos = grid.linearPos( i1, i2 );
            IndexType newPos[SCAI_GRID_MAX_DIMENSION];
            grid.gridPos( newPos, linearPos );
            BOOST_CHECK_EQUAL( i1, newPos[0] );
            BOOST_CHECK_EQUAL( i2, newPos[1] );

            // std::cout << "pos = " << i1 << ", " << i2 << " -> linearPos " << linearPos 
            //          << " -> " << newPos[0] << ", " << newPos[1] << std::endl;
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( pos3Test )
{
    const IndexType n1 = 5;
    const IndexType n2 = 3;
    const IndexType n3 = 4;

    Grid3D grid( n1, n2, n3 );

    const IndexType n = grid.size();

    IndexType p0[] = { 3, 1, 2, invalidIndex };
    IndexType p1[] = { 3, 1, 3, invalidIndex };
    IndexType p2[] = { 3, 2, 2, invalidIndex };
    IndexType p3[] = { 4, 1, 2, invalidIndex };

    // verify linearPos( x, y, z ) == linearPos( { x, y, z } )

    BOOST_CHECK_EQUAL( grid.linearPos( p1), grid.linearPos( p1[0], p1[1], p1[2] ) );

    // verify row-major ordering

    BOOST_CHECK_EQUAL( grid.linearPos( p0 ) + 1, grid.linearPos( p1 ) );
    BOOST_CHECK_EQUAL( grid.linearPos( p0 ) + n3, grid.linearPos( p2 ) );
    BOOST_CHECK_EQUAL( grid.linearPos( p0 ) + n2 * n3, grid.linearPos( p3 ) );

    for ( IndexType i = 0; i < n; ++i )
    {
        IndexType pos[3];
        grid.gridPos( pos, i );
        BOOST_CHECK_EQUAL( i, grid.linearPos( pos ) );
    }

    for ( IndexType i1 = 0; i1 < n1; ++i1 )
    {
        for ( IndexType i2 = 0; i2 < n2; ++i2 )
        {
            for ( IndexType i3 = 0; i3 < n3; ++i3 )
            {
                IndexType linearPos = grid.linearPos( i1, i2, i3 );
                IndexType newPos[3];
                grid.gridPos( newPos, linearPos );
                BOOST_CHECK_EQUAL( i1, newPos[0] );
                BOOST_CHECK_EQUAL( i2, newPos[1] );
                BOOST_CHECK_EQUAL( i3, newPos[2] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( pos4Test )
{
    const IndexType n1 = 5;
    const IndexType n2 = 3;
    const IndexType n3 = 4;
    const IndexType n4 = 2;

    Grid4D grid( n1, n2, n3, n4 );

    const IndexType n = grid.size();

    IndexType p0[] = { 3, 1, 2, 1 };
    IndexType p1[] = { 3, 1, 3, 1 };
    IndexType p2[] = { 3, 2, 2, 1  };
    IndexType p3[] = { 4, 1, 2, 1 };

    // verify linearPos( x, y, z ) == linearPos( { x, y, z } )

    BOOST_CHECK_EQUAL( grid.linearPos( p1), grid.linearPos( p1[0], p1[1], p1[2], p1[3] ) );

    // verify row-major ordering

    BOOST_CHECK_EQUAL( grid.linearPos( p0 ) + n4, grid.linearPos( p1 ) );
    BOOST_CHECK_EQUAL( grid.linearPos( p0 ) + n3 * n4, grid.linearPos( p2 ) );
    BOOST_CHECK_EQUAL( grid.linearPos( p0 ) + n2 * n3 * n4, grid.linearPos( p3 ) );

    for ( IndexType i = 0; i < n; ++i )
    {
        IndexType pos[4];
        grid.gridPos( pos, i );
        BOOST_CHECK_EQUAL( i, grid.linearPos( pos ) );
    }

    for ( IndexType i1 = 0; i1 < n1; ++i1 )
    {
        for ( IndexType i2 = 0; i2 < n2; ++i2 )
        {
            for ( IndexType i3 = 0; i3 < n3; ++i3 )
            {
                for ( IndexType i4 = 0; i4 < n4; ++i4 )
                {
                    IndexType linearPos = grid.linearPos( i1, i2, i3, i4 );
                    IndexType newPos[4];
                    grid.gridPos( newPos, linearPos );
                    BOOST_CHECK_EQUAL( i1, newPos[0] );
                    BOOST_CHECK_EQUAL( i2, newPos[1] );
                    BOOST_CHECK_EQUAL( i3, newPos[2] );
                    BOOST_CHECK_EQUAL( i4, newPos[3] );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( posNTest )
{
    const IndexType ndims = SCAI_GRID_MAX_DIMENSION;

    IndexType sizes[ SCAI_GRID_MAX_DIMENSION ];

    for ( IndexType i = 0; i < ndims; ++i )
    {
        // sizes = { 2, 3, 2, 3, 2, 3, .. };

        sizes[i] = ( i % 2 ) == 0 ? 2 : 3;
    }

    Grid grid( ndims, sizes );

    IndexType n = grid.size();

    for ( IndexType i = 0; i < n; ++i )
    {
        IndexType pos[SCAI_GRID_MAX_DIMENSION];
        grid.gridPos( pos, i );
        BOOST_CHECK_EQUAL( i, grid.linearPos( pos ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( equalTest )
{
    Grid1D grid1( 5 );
    Grid2D grid2a( 5, 3 );
    Grid2D grid2b( 3, 5 );
    Grid2D grid2c( 5, 3 );

    BOOST_CHECK( ! ( grid1 == grid2a ) );
    BOOST_CHECK( grid2a != grid2b );
    BOOST_CHECK_EQUAL( grid2a, grid2c );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeTest )
{
    std::ostringstream f;

    Grid3D grid( 5, 2, 4 );

    f << grid;

    const std::string& fstr = f.str();

    BOOST_CHECK( fstr.length() > 7 );

    BOOST_CHECK( fstr.find( "5" ) != std::string::npos );
    BOOST_CHECK( fstr.find( "3" ) != std::string::npos ); // 3-dim grid
    BOOST_CHECK( fstr.find( "4" ) != std::string::npos );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

