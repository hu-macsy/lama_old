/**
 * @file MaxNormTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Contains the implementation of the class MaxNormTest
 * @author Alexander Büchel
 * @date 21.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/lama/test/NormTest.hpp>
#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MaxNormTest )
;

SCAI_LOG_DEF_LOGGER( logger, "Test.MaxNormTest" );

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MaxNormVectorTests, ValueType, test_types )
{
    IndexType n = 4;
    ValueType val = 5.0;
    DenseVector<ValueType> vec( n, val );
    MaxNorm maxnorm;
    ValueType expected = 5.0;
    BOOST_CHECK_EQUAL( expected, maxnorm( vec ) );
    WriteAccess<ValueType> hwa( vec.getLocalValues() );
    hwa[0] = 1.0;
    hwa[1] = -2.0;
    hwa[2] = 3.0;
    hwa[3] = -4.5;
    hwa.release();
    expected = 4.5;
    BOOST_CHECK_EQUAL( expected, maxnorm( vec ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( MaxNormScalarTests )
{
    Scalar scalar( -4.0 );
    MaxNorm maxnorm;
    BOOST_CHECK_EQUAL( Scalar( 4.0 ), maxnorm( scalar ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    MaxNorm maxnorm;
    NormTest normtest( maxnorm );

    if ( base_test_case )
    {
        SCAI_LOG_INFO( logger, "Run test method " << testcase << " in MaxNormTest." );
        NORMTEST_COMMONTESTCASES( normtest );
    }
    else
    {
        normtest.runTests();
    }
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
