/**
 * @file TypeTraitsTest.cpp
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
 * @brief Test routines for class TypeTraits
 * @author Thomas Brandes
 * @date 05.02.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Printable.hpp>
#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace common;

#include <iomanip>

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_SUITE( TypeTraitsTest );

/* -----------------------------------------------------------------------------*/

/** For type traits test use all types for which they are available */

typedef boost::mpl::list<SCAI_ALL_TYPES> TraitTypes;

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE_TEMPLATE( IdTest, ValueType, TraitTypes )
{
    ScalarType stype = TypeTraits<ValueType>::stype;
    std::ostringstream out1;
    std::ostringstream out2;
    out1 << TypeTraits<ValueType>::id();
    out2 << stype;
    BOOST_CHECK_EQUAL( out1.str(), out2.str() );
}

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE_TEMPLATE( TypeSizeTest, ValueType, TraitTypes )
{
    // this test is essential as otherwise certain copy routines might fail

    ScalarType stype = TypeTraits<ValueType>::stype;
    BOOST_CHECK_EQUAL( sizeof( ValueType ), typeSize( stype ) );
}


/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE_TEMPLATE( IsComplexTest, ValueType, TraitTypes )
{
    // this test is essential as otherwise certain copy routines might fail

    ScalarType stype = TypeTraits<ValueType>::stype;

    bool v1 = isComplex( stype );

    std::ostringstream out;

    out << stype;

    std::size_t pos = out.str().find( "omplex" );

    bool v2 = pos != std::string::npos;

    BOOST_CHECK_EQUAL( v1, v2 );
}

/* -----------------------------------------------------------------------------*/

static int goodStrings( const std::string& str1, const std::string& str2 )
{
    // returns 3 if str1 = "0.333...33" and  str2 = "0.66...67"

    int val = 0;

    size_t p1 = str1.find_first_not_of( "3", 3 );

    // okay if there is nothing else than 3 after 0.3

    if ( p1 == std::string::npos )
    {
        val += 1;
    }

    // okay if every val is 6, only last char must be 7

    size_t p2 = str2.find_first_not_of( "6", 3 );

    if ( p2 == str2.size() - 1 )
    {
        // last char should be 7

        if ( str2[p2] == '7' )
        {
            val += 2;
        }
    }

    // std::cout << "goodStrings( " << str1 << ", " << str2 << " ) = " << val << std::endl;

    return val;
}

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE_TEMPLATE( PrecisionTest, ValueType, scai_numeric_test_types )
{
    if ( isComplex( TypeTraits<ValueType>::stype ) )
    {
        // skip this test for complex types
        return;
    }

    // this test can also not be applied for int, long

    ValueType x1 = ValueType( 1 ) / ValueType( 3 );
    ValueType x2 = ValueType( 2 ) / ValueType( 3 );

    int precision = TypeTraits<ValueType>::precision();

    std::ostringstream out1;
    std::ostringstream out2;

    out1 << std::setprecision( precision - 1 ) << x1;
    out2 << std::setprecision( precision - 1 ) << x2;

    // okay:  "0.333333" or "0.3333333333" or "0.333333333333333"
    // okay:  "0.666667" or "0.6666666667" or "0.66666666666667x"

    // if this test fails, output precision is too high

    BOOST_CHECK_EQUAL( goodStrings( out1.str(), out2.str() ), 3 );

    // remove content

    out1.str( "" );
    out2.str( "" );

    // now increase precision and it should fail

    out1 << std::setprecision( precision + 2 ) << x1;
    out2 << std::setprecision( precision + 2 ) << x2;

    BOOST_CHECK( goodStrings( out1.str(), out2.str() ) != 3 );
}

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE_TEMPLATE( SmallTest, ValueType, scai_numeric_test_types )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;

    RealType min    = TypeTraits<ValueType>::getMin();
    RealType max    = TypeTraits<ValueType>::getMax();
    RealType small  = TypeTraits<ValueType>::small();
    RealType eps0   = TypeTraits<ValueType>::eps0();
    RealType eps1   = TypeTraits<ValueType>::eps1();

    // neutral elements for min and max operations

    BOOST_CHECK( min < max );

    RealType one = 1;
    RealType one1 = one + eps1;

    // Note: This test fails for LongDouble or ComplexLongDouble when using valgrind

    BOOST_CHECK( one < one1 );
    BOOST_CHECK( one1 != RealType( 0 ) );

    // there should be no value between one and one1

    RealType one2 = ( one + one1 ) * RealType( 0.5 );

    BOOST_CHECK( one2 == one1 || one2 == one );

    // compare small, eps0 against eps1

    BOOST_CHECK( RealType( 0 ) < eps0 );
    BOOST_CHECK( eps0 < eps1 );
    BOOST_CHECK( eps1 < small );
}

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_CASE_TEMPLATE( ImaginaryTest, ValueType, scai_numeric_test_types )
{
    ValueType imagUnit = TypeTraits<ValueType>::imaginaryUnit();

    ValueType imagSquare = imagUnit * imagUnit;
   
    ScalarType stype = TypeTraits<ValueType>::stype;

    if ( isComplex( stype ) )
    {
        ValueType expected = -1;
        BOOST_CHECK_EQUAL( expected, imagSquare );
    }
    else
    {
        ValueType expected = 0;
        BOOST_CHECK_EQUAL( expected, imagSquare );
    }
}

/* -----------------------------------------------------------------------------*/

BOOST_AUTO_TEST_SUITE_END();
