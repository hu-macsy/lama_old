/**
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

using scai::common::Grid;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( GridTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    Grid grid1( 10 );
    BOOST_CHECK_EQUAL( 10, grid1.size() );
    Grid grid2( 2, 5 );
    BOOST_CHECK_EQUAL( 10, grid2.size() );
    Grid grid3( 2, 3, 2 );
    BOOST_CHECK_EQUAL( 12, grid3.size() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( pos2Test )
{
    const IndexType n1 = 10;
    const IndexType n2 = 15;

    Grid grid( n1, n2 );

    const IndexType n = grid.size();

    IndexType p1[] = { 3, 5 };
    IndexType p2[] = { 3, 6 };
    IndexType p3[] = { 4, 5 };

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
            IndexType pos[] = { i1, i2 };
            IndexType linearPos = grid.linearPos( pos );
            IndexType newPos[2];
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

    Grid grid( n1, n2, n3 );

    const IndexType n = grid.size();

    IndexType p0[] = { 3, 1, 2 };
    IndexType p1[] = { 3, 1, 3 };
    IndexType p2[] = { 3, 2, 2  };
    IndexType p3[] = { 4, 1, 2 };

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
                IndexType pos[] = { i1, i2, i3 };
                IndexType linearPos = grid.linearPos( pos );
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

BOOST_AUTO_TEST_SUITE_END();

