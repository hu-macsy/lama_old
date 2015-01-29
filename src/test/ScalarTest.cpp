/**
 * @file ScalarTest.cpp
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
 * @brief Contains the implementation of the class ScalarTest.
 * @author: Jiri Kraus, Eric Stricker
 * @date 21.06.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/Scalar.hpp>
#include <lama/Complex.hpp>

#include <test/TestMacros.hpp>

#include <complex>

using namespace lama;
using namespace boost;

// Scalar can be tested for all LAMA arithmetic types even if LAMA matrices
// and vectors have not been instantiated for these types

typedef boost::mpl::list<float, double, long double, 
                         ComplexFloat, ComplexDouble, ComplexLongDouble> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ScalarTest )

/* --------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( logger, "Test.ScalarTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScalarGetTypeTest )
{
    float value_s = 2.0;
    Scalar s( value_s );
    BOOST_CHECK_EQUAL( s.getType<float>(), Scalar::FLOAT );

    double value_t = 2.0;
    Scalar t( value_t );
    BOOST_CHECK_EQUAL( t.getType<double>(), Scalar::DOUBLE );

    LongDouble value_u = 2.0;
    Scalar u( value_u );
    BOOST_CHECK_EQUAL( u.getType<LongDouble>(), Scalar::LONG_DOUBLE );

    ComplexFloat value_c_v( 2.0, 1.0 );
    Scalar v( value_c_v );
    BOOST_CHECK_EQUAL( v.getType<ComplexFloat>(), Scalar::COMPLEX );

    ComplexDouble value_c_w( 2.0, 1.0 );
    Scalar w( value_c_w );
    BOOST_CHECK_EQUAL( w.getType<ComplexDouble>(), Scalar::DOUBLE_COMPLEX );

    ComplexLongDouble value_c_x( 2.0, 1.0 );
    Scalar x( value_c_x );
    BOOST_CHECK_EQUAL( x.getType<ComplexLongDouble>(), Scalar::LONG_DOUBLE_COMPLEX );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScalarTypeSizeTest )
{
    // Note: this test might be misleading but you should be aware
    //       that for a Scalar s you cannot determine how it has been constructed

    float value_float = 1.0f;
    double value_double = 1.0;
    LongDouble value_long_double = 1.0l;

    ComplexFloat value_c_float( 1.0, 2.0 );
    ComplexDouble value_c_double( 1.0, 2.0 );
    ComplexLongDouble value_c_long_double( 1.0, 2.0 );

    Scalar s_float( value_float );
    size_t size = s_float.getTypeSize( Scalar::FLOAT );
    BOOST_CHECK_EQUAL( size, sizeof( float ) );

    Scalar s_double( value_double );
    size = s_double.getTypeSize( Scalar::DOUBLE );
    BOOST_CHECK_EQUAL( size, sizeof( double ) );

    Scalar s_long_double( value_long_double );
    size = s_long_double.getTypeSize( Scalar::LONG_DOUBLE );
    BOOST_CHECK_EQUAL( size, sizeof( LongDouble) );

    Scalar s_float_c( value_c_float );
    size = s_float.getTypeSize( Scalar::COMPLEX );
    BOOST_CHECK_EQUAL( size, sizeof( ComplexFloat ) );

    Scalar s_double_c( value_c_double );
    size = s_double_c.getTypeSize( Scalar::DOUBLE_COMPLEX );
    BOOST_CHECK_EQUAL( size, sizeof( ComplexDouble ) );

    Scalar s_long_double_c( value_c_long_double );
    size = s_long_double_c.getTypeSize( Scalar::LONG_DOUBLE_COMPLEX );
    BOOST_CHECK_EQUAL( size, sizeof( ComplexLongDouble ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( AdditionTest, T, test_types )
{
    typedef T ValueType;

    Scalar s ( 2.0 );
    Scalar t ( 3.0 );

    Scalar u = s + t;
    s += t;

    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), 5.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 5.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MultiplicationTest, T, test_types )
{
    typedef T ValueType;

    Scalar s ( 2.0 );
    Scalar t ( 3.0 );

    Scalar u = s * t;
    s *= t;

    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), 6.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 6.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SubtractionTest, T, test_types )
{
    typedef T ValueType;

    Scalar s ( 2.0 );
    Scalar t ( 3.0 );

    Scalar u = s - t;
    s -= t;
    
    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), -1.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), -1.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( DivisionTest, T, test_types )
{
    Scalar s ( 2.0 );
    Scalar t ( 3.0 );

    Scalar u = s / t;
    s /= t;

    BOOST_CHECK_CLOSE( u.getValue<float>(), 0.6666, 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( IsRealTest )
{
    Scalar complexScalar = ComplexFloat( 1.0, 2.0 );
    BOOST_CHECK( !complexScalar.isReal() );

    Scalar realScalar( 2 );
    BOOST_CHECK( realScalar.isReal() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( NegativNumberTest, T, test_types )
{
    typedef T ValueType;

    Scalar s( 2.0 );
    Scalar t( -s );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 2.0 );
    BOOST_CHECK_EQUAL( t.getValue<ValueType>(), -2.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( EqualityTest )
{
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
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( MiscTests )
{
    Scalar s( 6.25 );
    Scalar t( 9.0 );
    Scalar u( -2.5 );

    BOOST_CHECK( sqrt( s ) == 2.5 );
    BOOST_CHECK_EQUAL( sqrt( s ), 2.5 );

    BOOST_CHECK_EQUAL( sqrt( t ), 3.0 );

    BOOST_CHECK_EQUAL( abs( u ), 2.5 );
    BOOST_CHECK_EQUAL( abs( t ), 9.0 );

    BOOST_CHECK_EQUAL( max( s, t ), 9.0  );
    BOOST_CHECK_EQUAL( min( s, t ), 6.25 );

    Scalar c1( 3.0, 4.0 );
    Scalar c2( 2.0, 2.0 );

    BOOST_CHECK_EQUAL( max( c1, c2 ), c1  );
    BOOST_CHECK_EQUAL( min( c1, c2 ), c2 );

    // Pythagoras: 3^2 + 4^2 = 5^2

    BOOST_CHECK_EQUAL( abs( c1 ), 5.0 );
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
    mStream << type;
    std::string mString = mStream.str();
    BOOST_CHECK_EQUAL( string, mString );
}

BOOST_AUTO_TEST_CASE( printTest )
{
    printtestmethod( "float", Scalar::FLOAT );
    printtestmethod( "double", Scalar::DOUBLE );
    printtestmethod( "LongDouble", Scalar::LONG_DOUBLE );
    printtestmethod( "ComplexFloat", Scalar::COMPLEX );
    printtestmethod( "ComplexDouble", Scalar::DOUBLE_COMPLEX );
    printtestmethod( "ComplexLongDouble", Scalar::LONG_DOUBLE_COMPLEX );
}

/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
