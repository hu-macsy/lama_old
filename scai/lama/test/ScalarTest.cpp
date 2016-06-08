/**
 * @file ScalarTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Contains the implementation of the class ScalarTest.
 * @author Jiri Kraus, Eric Stricker
 * @date 21.06.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/Scalar.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/lama/test/TestMacros.hpp>

//#include <complex>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::common;

// Scalar can be tested for all LAMA arithmetic types even if LAMA matrices
// and vectors have not been instantiated for these types

typedef boost::mpl::list< SCAI_ARITHMETIC_HOST> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ScalarTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.ScalarTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScalarGetTypeTest )
{
    // Challenge: some of these types were not defined in module common but later in lama
    using namespace scai::common;
    BOOST_CHECK_EQUAL( getScalarType<float>(), scalar::FLOAT );
    BOOST_CHECK_EQUAL( getScalarType<double>(), scalar::DOUBLE );
    BOOST_CHECK_EQUAL( getScalarType<LongDouble>(), scalar::LONG_DOUBLE );
#ifdef SCAI_COMPLEX_SUPPORTED
    BOOST_CHECK_EQUAL( getScalarType<ComplexFloat>(), scalar::COMPLEX );
    BOOST_CHECK_EQUAL( getScalarType<ComplexDouble>(), scalar::DOUBLE_COMPLEX );
    BOOST_CHECK_EQUAL( getScalarType<ComplexLongDouble>(), scalar::LONG_DOUBLE_COMPLEX );
#endif
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( AdditionTest, ValueType, test_types )
{
    Scalar s ( 2.0 );
    Scalar t ( 3.0 );
    Scalar u = s + t;
    s += t;
    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), 5.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 5.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MultiplicationTest, ValueType, test_types )
{
    Scalar s ( 2.0 );
    Scalar t ( 3.0 );
    Scalar u = s * t;
    s *= t;
    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), 6.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), 6.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SubtractionTest, ValueType, test_types )
{
    Scalar s ( 2.0 );
    Scalar t ( 3.0 );
    Scalar u = s - t;
    s -= t;
    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), -1.0 );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), -1.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( DivisionTest, ValueType, test_types )
{
    Scalar s ( 2.0 );
    Scalar t ( 3.0 );
    Scalar u = s / t;
    s /= t;
    SCAI_CHECK_CLOSE( u.getValue<ValueType>(), 0.6666, 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( NegativNumberTest, ValueType, test_types )
{
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
//    BOOST_CHECK( s <= t );
//    BOOST_CHECK( s <= u );
//    BOOST_CHECK( u >= s );
//    BOOST_CHECK( u >= t );
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
#ifdef SCAI_COMPLEX_SUPPORTED
    Scalar c1( ComplexFloat( 3.0, 4.0 ) );
    Scalar c2( ComplexFloat( 2.0, 2.0 ) );
    BOOST_CHECK_EQUAL( max( c1, c2 ), c1  );
    BOOST_CHECK_EQUAL( min( c1, c2 ), c2 );
    // Pythagoras: 3^2 + 4^2 = 5^2
    BOOST_CHECK_EQUAL( abs( c1 ), 5.0 );
#endif
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    Scalar s = 2.0;
    SCAI_COMMON_WRITEAT_TEST( s );
}

/* --------------------------------------------------------------------- */

void printtestmethod( std::string string, scalar::ScalarType type )
{
    std::stringstream mStream;
    mStream << type;
    std::string mString = mStream.str();
    BOOST_CHECK_EQUAL( string, mString );
}

BOOST_AUTO_TEST_CASE( printTest )
{
    printtestmethod( "float", scalar::FLOAT );
    printtestmethod( "double", scalar::DOUBLE );
    printtestmethod( "LongDouble", scalar::LONG_DOUBLE );
    printtestmethod( "ComplexFloat", scalar::COMPLEX );
    printtestmethod( "ComplexDouble", scalar::DOUBLE_COMPLEX );
    printtestmethod( "ComplexLongDouble", scalar::LONG_DOUBLE_COMPLEX );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
