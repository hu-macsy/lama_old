/**
 * @file ScalarTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Contains the implementation of the class ScalarTest.
 * @author: Alexander Büchel
 * @date 21.06.2012
 * $Id$
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/Scalar.hpp>

#include <test/TestMacros.hpp>

using namespace lama;
using namespace boost;

typedef boost::mpl::list<float,double,long double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ScalarTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.ScalarTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( InlineCtorTest, T, test_types ) {
    typedef T ValueType;

    ValueType value = 1.0;
    Scalar s( value );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 1.0 );

    std::complex<ValueType> cvalue( 1.0, 2.0 );
    Scalar t( cvalue );
    BOOST_CHECK_EQUAL( t.getValue<std::complex<ValueType> >() , cvalue );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScalarGetTypeTest )
{
    float value_s = 2.0;
    Scalar s( value_s );
    BOOST_CHECK_EQUAL( s.getType<float>(), Scalar::FLOAT );

    double value_t = 2.0;
    Scalar t( value_t );
    BOOST_CHECK_EQUAL( t.getType<double>(), Scalar::DOUBLE );

    long double value_u = 2.0;
    Scalar u( value_u );
    BOOST_CHECK_EQUAL( u.getType<long double>(), Scalar::LONG_DOUBLE );

    std::complex<float> value_c_v( 2.0, 1.0 );
    Scalar v( value_c_v );
    BOOST_CHECK_EQUAL( v.getType<std::complex<float> >(), Scalar::COMPLEX );

    std::complex<double> value_c_w( 2.0, 1.0 );
    Scalar w( value_c_w );
    BOOST_CHECK_EQUAL( w.getType<std::complex<double> >(), Scalar::DOUBLE_COMPLEX );

    std::complex<long double> value_c_x( 2.0, 1.0 );
    Scalar x( value_c_x );
    BOOST_CHECK_EQUAL( x.getType<std::complex<long double> >(), Scalar::LONG_DOUBLE_COMPLEX );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScalarTypeSizeTest )
{
    float value_float = 1.0f;
    double value_double = 1.0;
    long double value_long_double = 1.0l;

    std::complex<float> value_c_float( 1.0, 2.0 );
    std::complex<double> value_c_double( 1.0, 2.0 );
    std::complex<long double> value_c_long_double( 1.0, 2.0 );

    Scalar s_float( value_float );
    BOOST_CHECK_EQUAL( (int) s_float.getTypeSize( Scalar::FLOAT ), 4 );

    Scalar s_double( value_double );
    BOOST_CHECK_EQUAL( (int) s_double.getTypeSize( Scalar::DOUBLE ), 8 );

    Scalar s_long_double( value_long_double );
    BOOST_CHECK_EQUAL( (int) s_long_double.getTypeSize( Scalar::LONG_DOUBLE), 16 );

    Scalar s_float_c( value_c_float );
    BOOST_CHECK_EQUAL( (int) s_float.getTypeSize( Scalar::COMPLEX ), 8 );

    Scalar s_double_c( value_c_double );
    BOOST_CHECK_EQUAL( (int) s_double_c.getTypeSize( Scalar::DOUBLE_COMPLEX ), 16 );

    Scalar s_long_double_c( value_c_long_double );
    BOOST_CHECK_EQUAL( (int) s_long_double_c.getTypeSize( Scalar::LONG_DOUBLE_COMPLEX), 32 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( AdditionTest, T, test_types ) {
    typedef T ValueType;

    Scalar s ( 2.0 );
    Scalar t ( 3.0 );

    Scalar u = s + t;
    s += t;

    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), 5.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 5.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MultiplicationTest, T, test_types ) {
    typedef T ValueType;

    Scalar s ( 2.0 );
    Scalar t ( 3.0 );

    Scalar u = s * t;
    s *= t;

    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), 6.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 6.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SubtractionTest, T, test_types ) {
    typedef T ValueType;

    Scalar s ( 2.0 );
    Scalar t ( 3.0 );

    Scalar u = s - t;
    s -= t;

    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), -1.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), -1.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( DivisionTest, T, test_types ) {
    typedef T ValueType;

    Scalar s ( 2.0 );
    Scalar t ( 3.0 );

    Scalar u = s / t;
    s /= t;

    BOOST_CHECK_CLOSE( u.getValue<ValueType>(), 0.6666, 1 );
    BOOST_CHECK_CLOSE( s.getValue<ValueType>(), 0.6666, 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( IsRealTest, T, test_types ) {
    typedef T ValueType;

    std::complex<ValueType> cvalue( 1.0, 2.0 );
    Scalar s( cvalue );
    BOOST_CHECK( s.isReal() == false );

    Scalar t( 2.0 );
    BOOST_CHECK( t.isReal() == true );

    std::complex<ValueType> cvalue2( 1.0, 0.0 );
    Scalar u( cvalue2 );
    BOOST_CHECK( u.isReal() == true );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( NegativNumberTest, T, test_types ) {
    typedef T ValueType;

    Scalar s( 2.0 );
    Scalar t( -s );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 2.0 );
    BOOST_CHECK_EQUAL( t.getValue<ValueType>(), -2.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( EqualityTest, T, test_types ) {
    typedef T ValueType;

    std::complex<ValueType> c_value( 1.0, 2.0 );
    Scalar s ( 2.0 );
    Scalar t ( 2.0 );
    Scalar u ( 3.0 );

    BOOST_CHECK( s == t );
    BOOST_CHECK( s != u );
    BOOST_CHECK( s < u );
    BOOST_CHECK( u > s );
    BOOST_CHECK( s <= t );
    BOOST_CHECK( s <= u );
    BOOST_CHECK( u >= s );
    BOOST_CHECK( u >= t );

    LAMA_CHECK_THROW( c_value < s, Exception );
    LAMA_CHECK_THROW( c_value > s, Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( MiscTests )
{
    Scalar s( 4.0 );
    Scalar t( 9.0 );
    Scalar u( -2.0 );

    BOOST_CHECK( sqrt( s ) == 2.0 );
    BOOST_CHECK( sqrt( t ) == 3.0 );

    BOOST_CHECK( abs( u ) == 2.0 );
    BOOST_CHECK( abs( t ) == 9.0 );

    BOOST_CHECK( max( s, t ) == 9.0 );
    BOOST_CHECK( min( s, t ) == 4.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    Scalar s = 2.0;
    LAMA_WRITEAT_TEST( s );
}

/* --------------------------------------------------------------------- */

void printtestmethod( std::string string, Scalar::ScalarType type )
{
    std::stringstream mStream;
    operator<<( mStream, type );
    std::string mString = mStream.str();
    BOOST_CHECK_EQUAL( string, mString );
}

BOOST_AUTO_TEST_CASE( printTest )
{
    printtestmethod( "float", Scalar::FLOAT );
    printtestmethod( "double", Scalar::DOUBLE );
    printtestmethod( "long double", Scalar::LONG_DOUBLE );
    printtestmethod( "complex<float>", Scalar::COMPLEX );
    printtestmethod( "complex<double>", Scalar::DOUBLE_COMPLEX );
    printtestmethod( "complex<long double>", Scalar::LONG_DOUBLE_COMPLEX );
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
