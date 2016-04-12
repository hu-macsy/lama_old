/**
 * @file TypeLoopTest.hpp
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
 * @brief Test for TypeLoop
 * @author Eric Schricker
 * @date 12.04.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/macros/typeloop.hpp>
#include <scai/common/Complex.hpp>
#include <scai/common/SCAITypes.hpp>

#include <iostream>
#include <sstream>

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TypeLoopTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( outputTest )
{
#define TEST_TYPELOOP_LIST float, double, ComplexFloat

    std::stringstream s;

#define TEST_TYPELOOP_OUTPUT( type ) s << #type;

    SCAI_COMMON_TYPELOOP( 3, TEST_TYPELOOP_OUTPUT, TEST_TYPELOOP_LIST )

    BOOST_CHECK( "floatdoubleComplexFloat" == s.str() );

#undef TEST_TYPELOOP_OUTPUT
#undef TEST_TYPELOOP_LIST
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( countTest )
{
#define TEST_TYPELOOP_LIST float, double, ComplexFloat

    int i = 0;

#define TEST_TYPELOOP_INC( type ) ++i;

    SCAI_COMMON_TYPELOOP( 3, TEST_TYPELOOP_INC, TEST_TYPELOOP_LIST )

    BOOST_CHECK( 3 == i );

#undef TEST_TYPELOOP_INC
#undef TEST_TYPELOOP_LIST
}

/* --------------------------------------------------------------------- */

#define TEST_TYPELOOP_LIST float, double, ComplexFloat

#define TEST_TYPELOOP_DEF( type )       \
    type addOne( const type& x )        \
    {                                   \
        return x + 1 ;                  \
    }

namespace testing
{
    SCAI_COMMON_TYPELOOP( 3, TEST_TYPELOOP_DEF, TEST_TYPELOOP_LIST )
} /* end namespace testing */

BOOST_AUTO_TEST_CASE( defTest )
{
    float f = -2.0;
    double d = 3.31;
    ComplexFloat c = ComplexFloat( 3, 3.1 );


    BOOST_CHECK_CLOSE( -1.0, testing::addOne( f ), 1e-10 );
    BOOST_CHECK_CLOSE( 4.31, testing::addOne( d ), 1e-12 );
    BOOST_CHECK_CLOSE( 4.0, (testing::addOne( c )).real(), 1e-10 );
    BOOST_CHECK_CLOSE( 3.1f, (testing::addOne( c )).imag(), 1e-10 );

}
#undef TEST_TYPELOOP_DEF
#undef TEST_TYPELOOP_LIST

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
