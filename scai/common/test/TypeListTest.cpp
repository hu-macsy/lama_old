/**
 * @file TypeListTest.hpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Test for TypeLists
 * @author Eric Schricker
 * @date 11.04.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/mepr/TypeListUtils.hpp>

using namespace scai;
using namespace common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TypeListTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( containsTest )
{
#define TEST_LIST TYPELIST( 3, int, float, double )

    bool r = mepr::TypeListUtilsV<int, TEST_LIST>::contains;

    BOOST_CHECK( r );

    r = mepr::TypeListUtilsV<long long, TEST_LIST>::contains;
    BOOST_CHECK( !r );

#undef TEST_LIST
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( sizeTest )
{
#define TEST_LIST TYPELIST( 3, int, float, double )

    IndexType size = mepr::TypeListUtils<TEST_LIST>::size;

    BOOST_CHECK( 3 == size );

#undef TEST_LIST
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( indexTest )
{
#define TEST_LIST TYPELIST( 3, int, float, double )

    IndexType index = mepr::TypeListUtilsV<float, TEST_LIST>::index;

    BOOST_CHECK( 1 == index );

    index = mepr::TypeListUtilsV<size_t, TEST_LIST>::index;

    BOOST_CHECK( -1 == index );
#undef TEST_LIST
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
