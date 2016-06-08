/**
 * @file L1NormTest.cpp
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
 * @endlicense
 *
 * @brief Contains the implementation of the class L1NormTest
 * @author Alexander Büchel
 * @date 21.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/norm/L1Norm.hpp>

#include <scai/lama/test/NormTest.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( L1NormTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.L1NormTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( L1NormVectorTests, ValueType, scai_arithmetic_test_types )
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
