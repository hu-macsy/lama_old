/**
 * @file MathTest.cpp
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
 * @brief Test routines for class Math
 *
 * @author Eric Schricker
 * @date 07.04.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Complex.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/Math.hpp>

#include <scai/common/test/TestMacros.hpp>

using scai::common::Math;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MathTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sqrtTest, ValueType, scai_arithmetic_test_types )
{
    ValueType x = 9;
    ValueType sqrt_val = Math::sqrt( x );

    BOOST_CHECK( sqrt_val == 3.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absTest, ValueType, scai_arithmetic_test_types )
{
    typedef typename scai::common::TypeTraits<ValueType>::AbsType AbsType;

    ValueType x;
    AbsType abs_val;

    x = 9;
    abs_val = Math::abs( x );
    BOOST_CHECK( abs_val == x );

    x = -9;
    abs_val = Math::abs( x );
    BOOST_CHECK( abs_val == -x );
}

/* --------------------------------------------------------------------- */

// BOOST_AUTO_TEST_CASE_TEMPLATE( conjTest, ValueType, scai_arithmetic_test_types )
// is tested in ComplexTest

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( realTest, ValueType, scai_arithmetic_test_types )
{
    typedef typename scai::common::TypeTraits<ValueType>::AbsType AbsType;

    ValueType x;
    AbsType real_val;

    x = 9;
    real_val = Math::real( x );
    BOOST_CHECK( real_val == x );

    x = -9;
    real_val = Math::real( x );
    BOOST_CHECK( real_val == x );

    x = static_cast<ValueType>( ComplexFloat(3, 4) );
    real_val = Math::real( x );
    BOOST_CHECK( real_val == 3 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( imagTest, ValueType, scai_arithmetic_test_types )
{
    typedef typename scai::common::TypeTraits<ValueType>::AbsType AbsType;

    ValueType x;
    AbsType imag_val;

    x = 9;
    imag_val = Math::imag( x );
    BOOST_CHECK( imag_val == 0 );

    x = -9;
    imag_val = Math::imag( x );
    BOOST_CHECK( imag_val == 0 );

    x = static_cast<ValueType>( ComplexFloat(3, 4) );
    imag_val = Math::imag( x );

    if( scai::common::scalar::isComplex( scai::common::TypeTraits<ValueType>::stype ))
    {
        BOOST_CHECK( imag_val == 4 );
    }
    else
    {
        BOOST_CHECK( imag_val == 0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( minTest, ValueType, scai_arithmetic_test_types )
{
    ValueType x, y;
    ValueType min_val;

    x = 9;
    y = 10;
    min_val = Math::min( x, y );
    BOOST_CHECK( min_val == x );

    min_val = Math::min( y, x );
    BOOST_CHECK( min_val == x );

    y = 3;
    min_val = Math::min( x, y );
    BOOST_CHECK( min_val == y );

    min_val = Math::min( y, x );
    BOOST_CHECK( min_val == y );

    x = y;
    min_val = Math::min( x, y );
    BOOST_CHECK( min_val == x );
    BOOST_CHECK( min_val == y );

    min_val = Math::min( y, x );
    BOOST_CHECK( min_val == x );
    BOOST_CHECK( min_val == y );

}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( maxTest, ValueType, scai_arithmetic_test_types )
{
    ValueType x, y;
    ValueType max_val;

    x = 9;
    y = 10;
    max_val = Math::max( x, y );
    BOOST_CHECK( max_val == y );

    max_val = Math::max( y, x );
    BOOST_CHECK( max_val == y );

    y = 3;
    max_val = Math::max( x, y );
    BOOST_CHECK( max_val == x );

    max_val = Math::max( y, x );
    BOOST_CHECK( max_val == x );

    x = y;
    max_val = Math::max( x, y );
    BOOST_CHECK( max_val == x );
    BOOST_CHECK( max_val == y );

    max_val = Math::max( y, x );
    BOOST_CHECK( max_val == x );
    BOOST_CHECK( max_val == y );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( randomTest, ValueType, scai_arithmetic_test_types )
{
    ValueType random_val;

    Math::random( random_val );

    BOOST_CHECK( Math::abs( Math::real( random_val ) ) < 1.0 );
    BOOST_CHECK( Math::abs( Math::imag( random_val ) ) < 1.0 );
}

BOOST_AUTO_TEST_SUITE_END();

