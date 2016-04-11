/**
 * @file ComplexTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Test routines for class Complex
 *
 * @author Eric Schricker
 * @date 24.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/common/Complex.hpp>
#include <scai/common/Constants.hpp>

using namespace scai::common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ConstantsTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( zeroTest, ValueType, scai_arithmetic_test_types )
{
    ValueType zero( 0 );

    BOOST_CHECK( zero == constants::ZERO );

    BOOST_CHECK( constants::ZERO == zero );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( oneTest, ValueType, scai_arithmetic_test_types )
{
    ValueType one( 1 );

    BOOST_CHECK( one == constants::ONE );

    BOOST_CHECK( constants::ONE == one );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( randomTest, ValueType, scai_arithmetic_test_types )
{
    ValueType val;

    Math::random( val );

    BOOST_CHECK( val - val == constants::ZERO );

    BOOST_CHECK( constants::ZERO == val - val );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
