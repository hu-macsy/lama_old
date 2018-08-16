/**
 * @file ComplexTest.cpp
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

using namespace scai;
using namespace common;

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
    BOOST_CHECK_CLOSE( f.real(), 2, 0.01 );
    BOOST_CHECK_SMALL( f.imag(), TypeTraits<float>::eps1() );
    ComplexDouble d( assignee );
    BOOST_CHECK_CLOSE( d.real(), 2, 0.00001 );
    BOOST_CHECK_SMALL( d.imag(), TypeTraits<double>::eps1() );
    ComplexLongDouble l( assignee );
    BOOST_CHECK_CLOSE( l.real(), 2, 0.000001 );
    BOOST_CHECK_SMALL( l.imag(), TypeTraits<long double>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, ValueType, scai_complex_test_types )
{
    ValueType assignee( 2 );
    ComplexFloat f( 3 );
    f = assignee;
    BOOST_CHECK_CLOSE( f.real(), 2, 0.001 );
    BOOST_CHECK_SMALL( f.imag(), TypeTraits<float>::eps1() );
    ComplexDouble d( 4 );
    d = assignee;
    BOOST_CHECK_CLOSE( d.real(), 2, 0.00001 );
    BOOST_CHECK_SMALL( d.imag(), TypeTraits<double>::eps1() );
    ComplexLongDouble l( 5 );
    l = assignee;
    BOOST_CHECK_CLOSE( l.real(), 2, 0.0000001 );
    BOOST_CHECK_SMALL( l.imag(), TypeTraits<long double>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( changeSign, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    ValueType x( 2, -3 );
    ValueType y = -x;
    BOOST_CHECK_CLOSE( y.real(), -2, TypeTraits<RealType>::eps1() );
    BOOST_CHECK_CLOSE( y.imag(), 3, TypeTraits<RealType>::eps1() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( mult, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    // Complex numbers
    ValueType r, x, y;
    // real parts
    RealType a = 2;
    RealType b = -3;
    RealType c = 3;
    RealType d = 1;
    RealType re = ( a * c ) - ( b * d );
    RealType im = ( a * d ) + ( b * c );
    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );
    // operator *
    r = x * y;
    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<RealType>::eps1() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<RealType>::eps1() );
    // operator *=
    x *= y;
    BOOST_CHECK( r == x );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( add, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    // Complex numbers
    ValueType r, x, y;
    // real parts
    RealType a = 2;
    RealType b = -3;
    RealType c = 3;
    RealType d = 1;
    // result
    RealType re = a + c;
    RealType im = b + d;
    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );
    // operator +
    r = x + y;
    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<RealType>::eps1() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<RealType>::eps1() );
    // operator +=
    x += y;
    BOOST_CHECK( r == x );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sub, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    // Complex numbers
    ValueType r, x, y;
    // real parts
    RealType a = 2;
    RealType b = -3;
    RealType c = 3;
    RealType d = 1;
    // result
    RealType re = a - c;
    RealType im = b - d;
    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );
    // operator -
    r = x - y;
    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<RealType>::eps1() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<RealType>::eps1() );
    // operator -=
    x -= y;
    BOOST_CHECK( r == x );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( div, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    // Complex numbers
    ValueType r, x, y;
    // real parts
    RealType a = 2;
    RealType b = -3;
    RealType c = 3;
    RealType d = 1;
    RealType re = ( a * c + b * d ) / ( c * c + d * d );
    RealType im = ( b * c - a * d ) / ( c * c + d * d );
    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );
    // operator /
    r = x / y;
    BOOST_CHECK_CLOSE( r.real(), re, 0.000001 );
    BOOST_CHECK_CLOSE( r.imag(), im, 0.000001 );
    // operator /=
    x /= y;
    BOOST_CHECK_CLOSE( r.real(), x.real(), 0.0000001 );
    BOOST_CHECK_CLOSE( r.imag(), x.imag(), 0.0000001 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( output, ValueType, scai_complex_test_types )
{
    ValueType x( 2, -3 );
    std::stringstream s;
    s << x;
    BOOST_CHECK_EQUAL( s.str(), "2 -3" );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( input, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    ValueType x;
    // real and imaginary part set
    std::stringstream s1( "14 -893" );
    s1 >> x;
    BOOST_CHECK_CLOSE( x.real(), 14, 0.00001 );
    BOOST_CHECK_CLOSE( x.imag(), -893, 0.00001 );
    // just real part set
    std::stringstream s2( "23" );
    s2 >> x;
    BOOST_CHECK_CLOSE( x.real(), 23, 0.00001 );
    BOOST_CHECK_EQUAL( x.imag(), RealType( 0 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( conj, ValueType, scai_complex_test_types )
{
    // positiv real, negativ imag

    ValueType x( 2, -3 );
    ValueType r( 2, 3 );

    x = Math::conj( x );

    BOOST_CHECK_CLOSE( x.real(), r.real(), 0.1 );
    BOOST_CHECK_CLOSE( x.imag(), r.imag(), 0.1 );

    // negativ real, negativ imag

    r = ValueType( -14, -3 );
    x = ValueType( -14, 3 );
    x = Math::conj( x );
    BOOST_CHECK_CLOSE( x.real(), r.real(), 0.1 );
    BOOST_CHECK_CLOSE( x.imag(), r.imag(), 0.1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( abs, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    ValueType x( 3, -4 );
    RealType abs_val = Math::abs( x );
    RealType xa = 3;
    RealType xb = -4;
    RealType r = Math::sqrt( xa * xa + xb * xb );
    BOOST_CHECK_CLOSE( abs_val, r, 0.00001 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sqrt, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    ValueType x, sqrt_val, r;
    // sqrt from -1 --> i
    x = ValueType( -1, 0 );
    sqrt_val = Math::sqrt( x );
    r = ValueType( 0, 1 );
    BOOST_CHECK_SMALL( sqrt_val.real(), TypeTraits<RealType>::small() );
    BOOST_CHECK_CLOSE( sqrt_val.imag(), r.imag(), 0.00001 );
    // sqrt from 1 --> 1
    x = ValueType( 1, 0 );
    sqrt_val = Math::sqrt( x );
    r = ValueType( 1, 0 );
    BOOST_CHECK_CLOSE( sqrt_val.real(), r.real(), 0.00001 );
    BOOST_CHECK_SMALL( sqrt_val.imag(), TypeTraits<RealType>::small() );
    // sqrt from 3 + 4i --> 1
    x = ValueType( 3, 4 );
    sqrt_val = Math::sqrt( x );
    r = ValueType( 2, 1 );
    BOOST_CHECK_CLOSE( sqrt_val.real(), r.real(), 0.000001 );
    BOOST_CHECK_CLOSE( sqrt_val.imag(), r.imag(), 0.000001 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( compare, ValueType, scai_complex_test_types )
{
    ValueType x( 3, -4 );
    ValueType y( 3,  4 );
    ValueType z( 3, -4 );
    // Elementwise comparison
    BOOST_CHECK( x != y );
    BOOST_CHECK( x == z );
    BOOST_CHECK( x == Math::conj( y ) );
    BOOST_CHECK( x != Math::conj( z ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
