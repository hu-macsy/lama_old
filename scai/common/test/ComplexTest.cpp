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
    BOOST_CHECK_CLOSE( f.real(), 2, TypeTraits<float>::getEps() );
    BOOST_CHECK_CLOSE( f.imag(), 0, TypeTraits<float>::getEps() );

    ComplexDouble d( assignee );
    BOOST_CHECK_CLOSE( d.real(), 2, TypeTraits<double>::getEps() );
    BOOST_CHECK_CLOSE( d.imag(), 0, TypeTraits<double>::getEps() );

    ComplexLongDouble l( assignee );
    BOOST_CHECK_CLOSE( l.real(), 2, TypeTraits<long double>::getEps() );
    BOOST_CHECK_CLOSE( l.imag(), 0, TypeTraits<long double>::getEps() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignment, ValueType, scai_complex_test_types )
{
    ValueType assignee( 2 );

    ComplexFloat f( 3 );
    f = assignee;
    BOOST_CHECK_CLOSE( f.real(), 2, TypeTraits<float>::getEps() );
    BOOST_CHECK_CLOSE( f.imag(), 0, TypeTraits<float>::getEps() );

    ComplexDouble d( 4 );
    d = assignee;
    BOOST_CHECK_CLOSE( d.real(), 2, TypeTraits<double>::getEps() );
    BOOST_CHECK_CLOSE( d.imag(), 0, TypeTraits<double>::getEps() );

    ComplexLongDouble l( 5 );
    l = assignee;
    BOOST_CHECK_CLOSE( l.real(), 2, TypeTraits<long double>::getEps() );
    BOOST_CHECK_CLOSE( l.imag(), 0, TypeTraits<long double>::getEps() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( changeSign, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;

    ValueType x(2, -3);

    ValueType y = -x;
    BOOST_CHECK_CLOSE( y.real(), -2, TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( y.imag(), 3, TypeTraits<AbsType>::getEps() );
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

    AbsType re = (a * c) - (b * d);
    AbsType im = (a * d) + (b * c);

    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );

    // operator *
    r = x * y;

    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<AbsType>::getEps() );

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

    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<AbsType>::getEps() );

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

    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<AbsType>::getEps() );

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

    AbsType re = (a*c + b*d)/(c*c+d*d);
    AbsType im = (b*c - a*d)/(c*c+d*d);

    // x = a + bi, y = c + di
    x = ValueType( a, b );
    y = ValueType( c, d );

    // operator /
    r = x / y;

    BOOST_CHECK_CLOSE( r.real(), re, TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( r.imag(), im, TypeTraits<AbsType>::getEps() );

    // operator /=
    x /= y;

    BOOST_CHECK( r == x );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( output, ValueType, scai_complex_test_types )
{
    ValueType x(2, -3);

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
    std::stringstream s1("14 -893");
    s1 >> x;

    BOOST_CHECK_CLOSE( x.real(), 14, TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( x.imag(), -893, TypeTraits<AbsType>::getEps() );

    // just real part set
    std::stringstream s2("23");
    s2 >> x;

    BOOST_CHECK_CLOSE( x.real(), 23, TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( x.imag(), 0, TypeTraits<AbsType>::getEps() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( conj, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;

    ValueType x, r;

    // positiv real, negativ imag
    x = ValueType(2, -3);
    r = ValueType(2, 3);

    x = Math::conj( x );

    BOOST_CHECK_CLOSE( x.real(), r.real(), TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( x.imag(), r.imag(), TypeTraits<AbsType>::getEps() );

    // negativ real, negativ imag
    r = ValueType(-14, -3);
    x = ValueType(-14, 3);

    x = Math::conj( x );

    BOOST_CHECK_CLOSE( x.real(), r.real(), TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( x.imag(), r.imag(), TypeTraits<AbsType>::getEps() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( abs, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;

    ValueType x(3, -4);

    AbsType abs_val = Math::abs( x );

    AbsType xa = 3;
    AbsType xb = -4;

    AbsType r = Math::sqrt( xa * xa + xb * xb );

    BOOST_CHECK_CLOSE( abs_val, r, TypeTraits<ValueType>::getEps());
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sqrt, ValueType, scai_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::AbsType AbsType;

    ValueType x, sqrt_val, r;

    // sqrt from -1 --> i
    x = ValueType(-1, 0);

    sqrt_val = Math::sqrt( x );

    r = ValueType(0, 1);

    BOOST_CHECK_CLOSE( sqrt_val.real(), r.real(), TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( sqrt_val.imag(), r.imag(), TypeTraits<AbsType>::getEps() );

    // sqrt from 1 --> 1
    x = ValueType(1, 0);

    sqrt_val = Math::sqrt( x );

    r = ValueType(1, 0);

    BOOST_CHECK_CLOSE( sqrt_val.real(), r.real(), TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( sqrt_val.imag(), r.imag(), TypeTraits<AbsType>::getEps() );


    // sqrt from 3 + 4i --> 1
    x = ValueType(3, 4);

    sqrt_val = Math::sqrt( x );

    r = ValueType(2, 1);

    BOOST_CHECK_CLOSE( sqrt_val.real(), r.real(), TypeTraits<AbsType>::getEps() );
    BOOST_CHECK_CLOSE( sqrt_val.imag(), r.imag(), TypeTraits<AbsType>::getEps() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( min, ValueType, scai_complex_test_types )
{
    ValueType x, y, min_val;

    // min from -1 & 2 + i --> -1
    x = ValueType(-1, 0);
    y = ValueType( 2, 1);

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
    x = ValueType(-1, 0);
    y = ValueType( 2, 1);

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
    ValueType w(2,  4);
    ValueType x(3, -4);
    ValueType y(3,  4);
    ValueType z(3, -4);

    // Absolute value for comparison is used
    BOOST_CHECK( w < x );
//    BOOST_CHECK( w <= x );

    BOOST_CHECK( x > w );
//    BOOST_CHECK( x >= w );

//    BOOST_CHECK( x <= y );

//    BOOST_CHECK( x >= y );

    // ToDo: x<=y && x>=y should result in x == y
    //BOOST_CHECK( x == y );

    // Elementwise comparison
    BOOST_CHECK( x != y );

    BOOST_CHECK( x == z );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
