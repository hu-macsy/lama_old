/**
 * @file MaxNormTest.cpp
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
 * @brief Contains specific tests for derived norm class MaxNorm
 * @author Thomas Brandes
 * @date 21.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MaxNormTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.MaxNormTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MaxNormVectorTests, ValueType, scai_arithmetic_test_types )
{
    IndexType n = 4;
    ValueType val = 5.0;
    DenseVector<ValueType> vec( n, val );
    MaxNorm maxnorm;
    ValueType expected = 5.0;
    BOOST_CHECK_EQUAL( expected, maxnorm( vec ) );
    BOOST_CHECK_EQUAL( vec.maxNorm(), maxnorm( vec ) );
    WriteAccess<ValueType> hwa( vec.getLocalValues() );
    hwa[0] = 1.0;
    hwa[1] = -2.0;
    hwa[2] = 3.0;
    hwa[3] = -4.5;
    hwa.release();
    expected = 4.5;
    BOOST_CHECK_EQUAL( expected, maxnorm( vec ) );
    BOOST_CHECK_EQUAL( vec.maxNorm(), maxnorm( vec ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( MaxNormScalarTests )
{
    Scalar scalar( -4.0 );
    MaxNorm maxnorm;
    BOOST_CHECK_EQUAL( Scalar( 4.0 ), maxnorm( scalar ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
