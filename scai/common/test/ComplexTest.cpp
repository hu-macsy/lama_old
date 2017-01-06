/**
 * @file ComplexTest.cpp
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
 * @brief Test routines for class Complex
 * @author Eric Schricker
 * @date 24.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Complex.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/Math.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <string>
#include <sstream>

using namespace scai::common;

/* -------------------------------------------------------------------------------- */
/*  complex test types                                                              */
/* -------------------------------------------------------------------------------- */

typedef boost::mpl::list<ComplexFloat, ComplexDouble, ComplexLongDouble> scai_complex_test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ComplexTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructor, ValueType, scai_complex_test_types )
{
    ValueType assignee( 2 );
    ComplexFloat f( assignee );
    BOOST_CHECK_CLOSE( f.real(), 2, TypeTraits<float>::eps1() );
    BOOST_CHECK_CLOSE( f.imag(), 0, TypeTraits<float>::eps1() );
    ComplexDouble d( assignee );
    BOOST_CHECK_CLOSE( d.real(), 2, TypeTraits<double>::eps1() );
    BOOST_CHECK_CLOSE( d.imag(), 0, TypeTraits<double>::eps1() );
    ComplexLongDouble l( assignee );
    BOOST_CHECK_CLOSE( l.real(), 2, TypeTraits<long double>::eps1() );
    BOOST_CHECK_CLOSE( l.imag(), 0, TypeTraits<long double>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, ValueType, scai_complex_test_types )
{
    ValueType assignee( 2 );
    ComplexFloat f( 3 );
    f = assignee;
    BOOST_CHECK_CLOSE( f.real(), 2, TypeTraits<float>::eps1() );
    BOOST_CHECK_CLOSE( f.imag(), 0, TypeTraits<float>::eps1() );
    ComplexDouble d( 4 );
    d = assignee;
    BOOST_CHECK_CLOSE( d.real(), 2, TypeTraits<double>::eps1() );
    BOOST_CHECK_CLOSE( d.imag(), 0, TypeTraits<double>::eps1() );
    ComplexLongDouble l( 5 );
    l = assignee;
    BOOST_CHECK_CLOSE( l.real(), 2, TypeTraits<long double>::eps1() );
    BOOST_CHECK_CLOSE( l.imag(), 0, TypeTraits<long double>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( changeSign, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    ValueType x( 2, -3 );
    ValueType y = -x;
    BOOST_CHECK_CLOSE( y.real(), -2, TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( y.imag(), 3, TypeTraits<AbsType>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( mult, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    // Complex numbers
    ValueType r, x, y;
    // real parts
    AbsType a = 2;
    AbsType b = -3;
    AbsType c = 3;
    AbsType d = 1;
    AbsType re = ( a * c ) - ( b * d );
    AbsType im = ( a * d ) + ( b * c );
    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );
    // operator *
    r = x * y;
    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<AbsType>::eps1() );
    // operator *=
    x *= y;
    BOOST_CHECK( r == x );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( add, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    // Complex numbers
    ValueType r, x, y;
    // real parts
    AbsType a = 2;
    AbsType b = -3;
    AbsType c = 3;
    AbsType d = 1;
    // result
    AbsType re = a + c;
    AbsType im = b + d;
    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );
    // operator +
    r = x + y;
    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<AbsType>::eps1() );
    // operator +=
    x += y;
    BOOST_CHECK( r == x );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sub, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    // Complex numbers
    ValueType r, x, y;
    // real parts
    AbsType a = 2;
    AbsType b = -3;
    AbsType c = 3;
    AbsType d = 1;
    // result
    AbsType re = a - c;
    AbsType im = b - d;
    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );
    // operator -
    r = x - y;
    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<AbsType>::eps1() );
    // operator -=
    x -= y;
    BOOST_CHECK( r == x );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( div, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    // Complex numbers
    ValueType r, x, y;
    // real parts
    AbsType a = 2;
    AbsType b = -3;
    AbsType c = 3;
    AbsType d = 1;
    AbsType re = ( a * c + b * d ) / ( c * c + d * d );
    AbsType im = ( b * c - a * d ) / ( c * c + d * d );
    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );
    // operator /
    r = x / y;
    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<AbsType>::eps1() );
    // operator /=
    x /= y;
    BOOST_CHECK( r == x );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( output, ValueType, scai_complex_test_types )
{
    ValueType x( 2, -3 );
    std::stringstream s;
    s << x;
    BOOST_CHECK( s.str() == "2 -3" );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( input, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    ValueType x;
    // real and imaginary part set
    std::stringstream s1( "14 -893" );
    s1 >> x;
    BOOST_CHECK_CLOSE( x.real(), 14, TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( x.imag(), -893, TypeTraits<AbsType>::eps1() );
    // just real part set
    std::stringstream s2( "23" );
    s2 >> x;
    BOOST_CHECK_CLOSE( x.real(), 23, TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( x.imag(), 0, TypeTraits<AbsType>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( conj, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    ValueType x, r;
    // positiv real, negativ imag
    x = ValueType( 2, -3 );
    r = ValueType( 2, 3 );
    x = Math::conj( x );
    BOOST_CHECK_CLOSE( x.real(), r.real(), TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( x.imag(), r.imag(), TypeTraits<AbsType>::eps1() );
    // negativ real, negativ imag
    r = ValueType( -14, -3 );
    x = ValueType( -14, 3 );
    x = Math::conj( x );
    BOOST_CHECK_CLOSE( x.real(), r.real(), TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( x.imag(), r.imag(), TypeTraits<AbsType>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( abs, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    ValueType x( 3, -4 );
    AbsType abs_val = Math::abs( x );
    AbsType xa = 3;
    AbsType xb = -4;
    AbsType r = Math::sqrt( xa * xa + xb * xb );
    BOOST_CHECK_CLOSE( abs_val, r, TypeTraits<ValueType>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sqrt, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;
    ValueType x, sqrt_val, r;
    // sqrt from -1 --> i
    x = ValueType( -1, 0 );
    sqrt_val = Math::sqrt( x );
    r = ValueType( 0, 1 );
    BOOST_CHECK_CLOSE( sqrt_val.real(), r.real(), TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( sqrt_val.imag(), r.imag(), TypeTraits<AbsType>::eps1() );
    // sqrt from 1 --> 1
    x = ValueType( 1, 0 );
    sqrt_val = Math::sqrt( x );
    r = ValueType( 1, 0 );
    BOOST_CHECK_CLOSE( sqrt_val.real(), r.real(), TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( sqrt_val.imag(), r.imag(), TypeTraits<AbsType>::eps1() );
    // sqrt from 3 + 4i --> 1
    x = ValueType( 3, 4 );
    sqrt_val = Math::sqrt( x );
    r = ValueType( 2, 1 );
    BOOST_CHECK_CLOSE( sqrt_val.real(), r.real(), TypeTraits<AbsType>::eps1() );
    BOOST_CHECK_CLOSE( sqrt_val.imag(), r.imag(), TypeTraits<AbsType>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( min, ValueType, scai_complex_test_types )
{
    ValueType x, y, min_val;
    // min from -1 & 2 + i --> -1
    x = ValueType( -1, 0 );
    y = ValueType( 2, 1 );
    min_val = Math::min( x, y );
    BOOST_CHECK( min_val == x );
    min_val = Math::min( y, x );
    BOOST_CHECK( min_val == x );
    // min from 7 + i & -1 + 5i --> -1 + 5i
    x = ValueType( 7, 1 );
    y = ValueType( -1, 5 );
    min_val = Math::min( x, y );
    BOOST_CHECK( min_val == y );
    min_val = Math::min( y, x );
    BOOST_CHECK( min_val == y );
    // min from -1 + 5i & -1 + 5i
    x = ValueType( -1, 5 );
    y = ValueType( -1, 5 );
    min_val = Math::min( x, y );
    BOOST_CHECK( min_val == x );
    BOOST_CHECK( min_val == y );
    min_val = Math::min( y, x );
    BOOST_CHECK( min_val == x );
    BOOST_CHECK( min_val == y );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( max, ValueType, scai_complex_test_types )
{
    ValueType x, y, max_val;
    // min from -1 & 2 + i --> -1
    x = ValueType( -1, 0 );
    y = ValueType( 2, 1 );
    max_val = Math::max( x, y );
    BOOST_CHECK( max_val == y );
    max_val = Math::max( y, x );
    BOOST_CHECK( max_val == y );
    // min from 7 + i & -1 + 5i --> -1 + 5i
    x = ValueType( 7, 1 );
    y = ValueType( -1, 5 );
    max_val = Math::max( x, y );
    BOOST_CHECK( max_val == x );
    max_val = Math::max( y, x );
    BOOST_CHECK( max_val == x );
    // min from -1 + 5i & -1 + 5i
    x = ValueType( -1, 5 );
    y = ValueType( -1, 5 );
    max_val = Math::max( x, y );
    BOOST_CHECK( max_val == x );
    BOOST_CHECK( max_val == y );
    max_val = Math::max( y, x );
    BOOST_CHECK( max_val == x );
    BOOST_CHECK( max_val == y );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( compare, ValueType, scai_complex_test_types )
{
    ValueType w( 2,  4 );
    ValueType x( 3, -4 );
    ValueType y( 3,  4 );
    ValueType z( 3, -4 );
    // Absolute value for comparison is used
    BOOST_CHECK( w < x );
    BOOST_CHECK( x > w );
    // Elementwise comparison
    BOOST_CHECK( x != y );
    BOOST_CHECK( x == z );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
