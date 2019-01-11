/**
 * @file test/TypeListTest.cpp
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
 * @brief Test for TypeLists
 * @author Eric Schricker
 * @date 11.04.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/mepr/TypeListUtils.hpp>

using namespace scai;
using namespace common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TypeListTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( containsTest )
{
#define TEST_LIST SCAI_TYPELIST( int, float, double )
    bool r = mepr::TypeListUtilsV<int, TEST_LIST>::contains;
    BOOST_CHECK( r );
    r = mepr::TypeListUtilsV<long long, TEST_LIST>::contains;
    BOOST_CHECK( !r );
#undef TEST_LIST
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( sizeTest )
{
#define TEST_LIST SCAI_TYPELIST( int, float, double )
    int size = mepr::TypeListUtils<TEST_LIST>::size;
    BOOST_CHECK( 3 == size );
#undef TEST_LIST
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( indexTest )
{
#define TEST_LIST SCAI_TYPELIST( int, float, double )
    int index = mepr::TypeListUtilsV<float, TEST_LIST>::index;
    BOOST_CHECK( 1 == index );
    index = mepr::TypeListUtilsV<size_t, TEST_LIST>::index;
    BOOST_CHECK( -1 == index );
#undef TEST_LIST
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
