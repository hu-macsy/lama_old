/**
 * @file L2NormTest.cpp
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
 * @brief Contains the implementation of the class L2NormTest
 * @author Alexander Büchel, Micha
 * @date 03.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>
#include <lama/norm/L2Norm.hpp>

#include <test/NormTest.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( L2NormTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.L2NormTest" )

typedef boost::mpl::list<float,double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( L2NormVectorTests, T, test_types )
{
    typedef T ValueType;

    IndexType n = 4;
    ValueType val = 5.0;

    DenseVector<ValueType> vec( n, val );
    L2Norm l2norm;

    ValueType expected = std::sqrt( n * val * val );

    BOOST_CHECK_EQUAL( expected, l2norm( vec ) );

    HostWriteAccess<ValueType> hwa( vec.getLocalValues() );
    hwa[0] = 1.0;
    hwa[1] = -2.0;
    hwa[2] = 3.0;
    hwa[3] = -4.0;
    hwa.release();

    expected = static_cast<ValueType>( 5.47722 );

    Scalar s = l2norm( vec );

    BOOST_CHECK_CLOSE( expected, s.getValue<ValueType>(), 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( L2NormScalarTests )
{
    Scalar scalar( -4.0 );
    L2Norm l2norm;

    BOOST_CHECK_EQUAL( Scalar( 4.0 ), l2norm( scalar ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    L2Norm l2norm;
    NormTest normtest( l2norm );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in L2NormTest." );
        NORMTEST_COMMONTESTCASES( normtest );
    }
    else
    {
        normtest.runTests();
    }
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
