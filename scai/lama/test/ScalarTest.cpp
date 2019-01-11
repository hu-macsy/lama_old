/**
 * @file ScalarTest.cpp
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
 * @brief Contains the implementation of the class ScalarTest.
 * @author Jiri Kraus, Eric Stricker
 * @date 21.06.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/Scalar.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/test/TestMacros.hpp>

//#include <complex>

using namespace scai;
using namespace lama;
using namespace common;
using intern::Scalar;

// Scalar can be tested for all LAMA arithmetic types even if LAMA matrices
// and vectors have not been instantiated for these types

typedef boost::mpl::list< SCAI_NUMERIC_TYPES_HOST> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ScalarTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.ScalarTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScalarGetTypeTest )
{
    // Challenge: some of these types were not defined in module common but later in lama
    using namespace scai::common;
    BOOST_CHECK_EQUAL( getScalarType<float>(), ScalarType::FLOAT );
    BOOST_CHECK_EQUAL( getScalarType<double>(), ScalarType::DOUBLE );
    BOOST_CHECK_EQUAL( getScalarType<LongDouble>(), ScalarType::LONG_DOUBLE );
#ifdef SCAI_COMPLEX_SUPPORTED
    BOOST_CHECK_EQUAL( getScalarType<ComplexFloat>(), ScalarType::COMPLEX );
    BOOST_CHECK_EQUAL( getScalarType<ComplexDouble>(), ScalarType::DOUBLE_COMPLEX );
    BOOST_CHECK_EQUAL( getScalarType<ComplexLongDouble>(), ScalarType::LONG_DOUBLE_COMPLEX );
#endif
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( AdditionTest, ValueType, test_types )
{
    Scalar s ( 2.0 );
    Scalar t ( 3.0 );
    Scalar u = s + t;
    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), 5.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MultiplicationTest, ValueType, test_types )
{
    Scalar s ( 2 );
    Scalar t ( 3 );
    Scalar u = s * t;
    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), ValueType( 6 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SubtractionTest, ValueType, test_types )
{
    Scalar s ( 2.0 );
    Scalar t ( 3.0 );
    Scalar u = s - t;
    BOOST_CHECK_EQUAL( u.getValue<ValueType>(), -1.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( DivisionTest, ValueType, test_types )
{
    Scalar s ( 2.0 );
    Scalar t ( 3.0 );
    Scalar u = s / t;

    ValueType v1 = u.getValue<ValueType>();
    ValueType v2 = ValueType( 2 ) / ValueType( 3 );
 
    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    BOOST_CHECK( common::Math::abs( v1 - v2 ) < eps );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( NegativNumberTest, ValueType, test_types )
{
    Scalar s( 2.0 );
    Scalar t( -s );
    BOOST_CHECK_EQUAL( s.getValue<ValueType>(), ValueType( 2 ) );
    BOOST_CHECK_EQUAL( t.getValue<ValueType>(), ValueType( -2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    Scalar s = 2.0;
    SCAI_COMMON_WRITEAT_TEST( s );
}

/* --------------------------------------------------------------------- */

void printtestmethod( std::string string, ScalarType type )
{
    std::stringstream mStream;
    mStream << type;
    std::string mString = mStream.str();
    BOOST_CHECK_EQUAL( string, mString );
}

BOOST_AUTO_TEST_CASE( printTest )
{
    printtestmethod( "float", ScalarType::FLOAT );
    printtestmethod( "double", ScalarType::DOUBLE );
    printtestmethod( "LongDouble", ScalarType::LONG_DOUBLE );
    printtestmethod( "ComplexFloat", ScalarType::COMPLEX );
    printtestmethod( "ComplexDouble", ScalarType::DOUBLE_COMPLEX );
    printtestmethod( "ComplexLongDouble", ScalarType::LONG_DOUBLE_COMPLEX );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
