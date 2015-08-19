/**
 * @file L1NormTest.hpp
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
 * @brief Contains the implementation of the class L1NormTest
 * @author Alexander Büchel
 * @date 21.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/norm/L1Norm.hpp>

#include <test/NormTest.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( L1NormTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.L1NormTest" )

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( L1NormVectorTests, ValueType, test_types )
{
    IndexType n = 4;
    ValueType val = 5.0;
    DenseVector<ValueType> vec( n, val );
    L1Norm l1norm;
    ValueType expected = n * val;
    BOOST_CHECK_EQUAL( expected, l1norm( vec ) );
    WriteAccess<ValueType> hwa( vec.getLocalValues() );
    hwa[0] = 1.0;
    hwa[1] = -2.0;
    hwa[2] = 3.0;
    hwa[3] = -4.0;
    hwa.release();
    expected = 10.0;
    BOOST_CHECK_EQUAL( expected, l1norm( vec ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( L1NormScalarTests )
{
    Scalar scalar( -4.0 );
    L1Norm l1norm;
    BOOST_CHECK_EQUAL( Scalar( 4.0 ), l1norm( scalar ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    L1Norm l1norm;
    NormTest normtest( l1norm );

    if ( base_test_case )
    {
        SCAI_LOG_INFO( logger, "Run method " << testcase << " in L1NormTest." );
        NORMTEST_COMMONTESTCASES( normtest );
    }
    else
    {
        normtest.runTests();
    }
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
