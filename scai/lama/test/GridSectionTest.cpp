/**
 * @file GridSectionTest.cpp
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
 * @brief Contains specific tests for class GridSection
 * @author Thomas Brandes
 * @date 16.05.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/lama/GridSection.hpp>
#include <scai/lama/GridVector.hpp>

using namespace scai;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( GridSectionTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.GridSectionTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( cTorTest, ValueType, scai_numeric_test_types )
{
    const IndexType n1 = 10;
    const IndexType n2 = 4;

    const common::Grid2D grid( n1, n2 );

    GridVector<ValueType> gv( grid, ValueType( 2 ) );

    gv( Range(), 0 ) += ValueType( 1 );
    gv( Range(), 0 ) *= gv( Range(), 1 );

    // ValueType x = gv( 5, 0 );
    ValueType x = 6;

    BOOST_CHECK_EQUAL( x, 6 );  // ( 2 + 1 ) * 2
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
