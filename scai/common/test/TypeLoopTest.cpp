/**
 * @file test/TypeLoopTest.cpp
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
 * @brief Test for TypeLoop
 * @author Eric Schricker
 * @date 12.04.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/macros/loop.hpp>
#include <scai/common/SCAITypes.hpp>

#ifdef SCAI_COMPLEX_SUPPORTED
#include <scai/common/Complex.hpp>

#define TEST_TYPELOOP_LIST float, double, scai::ComplexFloat
#define TEST_STRING "floatdoublescai::ComplexFloat"
#define TEST_N_TYPES 3
#else
#define TEST_TYPELOOP_LIST float, double
#define TEST_STRING "floatdouble"
#define TEST_N_TYPES 2
#endif

#include <iostream>
#include <sstream>

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TypeLoopTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( outputTest )
{
    std::stringstream s;
#define TEST_TYPELOOP_OUTPUT( type ) s << #type;
    SCAI_COMMON_LOOP( TEST_TYPELOOP_OUTPUT, TEST_TYPELOOP_LIST )
    BOOST_CHECK_EQUAL( TEST_STRING, s.str() );
#undef TEST_TYPELOOP_OUTPUT
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( countTest )
{
    int i = 0;
#define TEST_TYPELOOP_INC( type ) ++i;
    SCAI_COMMON_LOOP( TEST_TYPELOOP_INC, TEST_TYPELOOP_LIST )
    BOOST_CHECK( TEST_N_TYPES == i );
#undef TEST_TYPELOOP_INC
}

/* --------------------------------------------------------------------- */

#define TEST_TYPELOOP_DEF( type )       \
    type addOne( const type& x )        \
    {                                   \
        return x + 1 ;                  \
    }

namespace testing
{
SCAI_COMMON_LOOP( TEST_TYPELOOP_DEF, TEST_TYPELOOP_LIST )
} /* end namespace testing */

BOOST_AUTO_TEST_CASE( defTest )
{
    float f = -2.0;
    double d = 3.31;
    BOOST_CHECK_CLOSE( -1.0, testing::addOne( f ), 1e-10 );
    BOOST_CHECK_CLOSE( 4.31, testing::addOne( d ), 1e-12 );
#ifdef SCAI_COMPLEX_SUPPORTED
    scai::ComplexFloat c( 3, 3.1 );
    BOOST_CHECK_CLOSE( 4.0, ( testing::addOne( c ) ).real(), 1e-10 );
    BOOST_CHECK_CLOSE( 3.1f, ( testing::addOne( c ) ).imag(), 1e-10 );
#endif
#undef TEST_TYPELOOP_DEF
}

/* --------------------------------------------------------------------- */

#undef TEST_TYPELOOP_LIST
#undef TEST_STRING

BOOST_AUTO_TEST_SUITE_END();
