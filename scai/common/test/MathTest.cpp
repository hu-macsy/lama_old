/**
 * @file MathTest.cpp
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
 * @brief Test routines for class Math
 * @author Eric Schricker
 * @date 07.04.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/Utils.hpp>

#include <scai/common/test/TestMacros.hpp>

using scai::common::Math;
using scai::common::TypeTraits;
using scai::common::Utils;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MathTest )

typedef boost::mpl::list<float, double, long double> scai_math_scalar_test_types;
#ifdef SCAI_COMPLEX_SUPPORTED
typedef boost::mpl::list<scai::common::Complex<float>, scai::common::Complex<double>, scai::common::Complex<long double> > scai_math_complex_test_types;
#endif

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sqrtTest, ValueType, scai_numeric_test_types )
{
    ValueType x = 9;
    ValueType sqrt_val = Math::sqrt( x );
    BOOST_CHECK( sqrt_val == 3.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absTest, ValueType, scai_numeric_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    ValueType x;
    RealType abs_val;
    x = 9;
    abs_val = Math::abs( x );
    BOOST_CHECK( abs_val == x );
    x = -9;
    abs_val = Math::abs( x );
    BOOST_CHECK( abs_val == -x );
}

/* --------------------------------------------------------------------- */

// BOOST_AUTO_TEST_CASE_TEMPLATE( conjTest, ValueType, scai_numeric_test_types )
// is tested in ComplexTest

/* --------------------------------------------------------------------- */

#ifdef SCAI_COMPLEX_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE( realTest, ValueType, scai_numeric_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    ValueType x;
    RealType real_val;
    x = 9;
    real_val = Math::real( x );
    BOOST_CHECK( real_val == x );
    x = -9;
    real_val = Math::real( x );
    BOOST_CHECK( real_val == x );
    x = static_cast<ValueType>( scai::ComplexFloat( 3, 4 ) );
    real_val = Math::real( x );
    BOOST_CHECK( real_val == 3 );
}

#endif

/* --------------------------------------------------------------------- */

#ifdef SCAI_COMPLEX_SUPPORTED

BOOST_AUTO_TEST_CASE_TEMPLATE( imagTest, ValueType, scai_numeric_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    ValueType x;
    RealType imag_val;
    x = 9;
    imag_val = Math::imag( x );
    BOOST_CHECK( imag_val == 0 );
    x = -9;
    imag_val = Math::imag( x );
    BOOST_CHECK( imag_val == 0 );
    x = static_cast<ValueType>( scai::ComplexFloat( 3, 4 ) );
    imag_val = Math::imag( x );

    if ( scai::common::isComplex( TypeTraits<ValueType>::stype ) )
    {
        BOOST_CHECK( imag_val == 4 );
    }
    else
    {
        BOOST_CHECK( imag_val == 0 );
    }
}

#endif

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( minTest, ValueType, scai_numeric_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealValueType;

    RealValueType x = 9;
    RealValueType y = 10;

    RealValueType min_val = Math::min( x, y );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( maxTest, ValueType, scai_numeric_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealValueType;

    RealValueType x = 9;
    RealValueType y = 10;

    RealValueType max_val = Math::max( x, y );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( randomTest, ValueType, scai_numeric_test_types )
{
    ValueType random_val = Math::random<ValueType>( 1 );

    BOOST_CHECK( Math::abs( Math::real( random_val ) ) < 1.0 );
    BOOST_CHECK( Math::abs( Math::imag( random_val ) ) < 1.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( powTest, ValueType, scai_numeric_test_types )
{
    ValueType x1 = 2;

    ValueType x1p10 = Math::pow( x1, ValueType( 10 ) );

    BOOST_CHECK( Math::abs( x1p10 - ValueType( 1024 ) ) < 0.01 );

    ValueType y = 16;

    ValueType y_r = Math::pow( y, ValueType( 0.5 ) );

    BOOST_CHECK( Math::abs( y_r - ValueType( 4 ) ) < 0.01 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( atan2Test, ValueType, scai_math_scalar_test_types )
{
    ValueType x = -10;
    ValueType y = 10;

    ValueType r = Math::atan2( y, x ) * ValueType( 180 ) / ValueType( M_PI );

    BOOST_CHECK( Math::abs( r - ValueType( 135 ) ) < 0.01 );
}

/* --------------------------------------------------------------------- */

#ifdef SCAI_COMPLEX_SUPPORTED
BOOST_AUTO_TEST_CASE_TEMPLATE( argTest, ValueType, scai_math_complex_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    ValueType x = ValueType( 3, 4 );

    RealType y = Math::arg( x );

    BOOST_CHECK_SMALL( y - RealType( 0.927295218001612 ), TypeTraits<ValueType>::small() );
}
#endif

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( signTest, ValueType, scai_numeric_test_types )
{
    BOOST_CHECK_EQUAL( Math::sign( ValueType( 0.2 ) ), ValueType( 1 ) );
    BOOST_CHECK_EQUAL( Math::sign( ValueType( -0.1 ) ), ValueType( -1 ) );
    BOOST_CHECK_EQUAL( Math::sign( ValueType( 0 ) ), ValueType( 0 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( logTest, ValueType, scai_numeric_test_types )
{
    ValueType x = ValueType( 10 );

    ValueType y = Math::log( x );

    BOOST_CHECK( Math::abs( y - ValueType( 2.30258509299405 ) ) < 0.001 );
}

/* --------------------------------------------------------------------- */

typedef boost::mpl::list<short, unsigned short, int, unsigned int, long, unsigned long> IndexTypes;

BOOST_AUTO_TEST_CASE( nextpow2Test )
{
    using scai::IndexType;

    // Here we make the test with IndexType

    BOOST_CHECK_EQUAL( Math::nextpow2( IndexType( 31 ) ), IndexType( 5 ) );
    BOOST_CHECK_EQUAL( Math::nextpow2( IndexType( 32 ) ), IndexType( 5 ) );
    BOOST_CHECK_EQUAL( Math::nextpow2( IndexType( 33 ) ), IndexType( 6 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( validIndexTest, IndexType, IndexTypes )
{
    const IndexType size = 13;
    IndexType val = 0;

    BOOST_CHECK( Utils::validIndex( val, size ) );

    val = size - 1;
    BOOST_CHECK( Utils::validIndex( val, size ) );

    val = size;
    BOOST_CHECK( !Utils::validIndex( val, size ) );

    val = static_cast<IndexType>( -1 );
    BOOST_CHECK( !Utils::validIndex( val, size ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

