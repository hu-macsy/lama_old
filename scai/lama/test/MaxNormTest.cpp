/**
 * @file MaxNormTest.cpp
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
 * @brief Contains specific tests for derived norm class MaxNorm
 * @author Thomas Brandes
 * @date 21.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MaxNormTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.MaxNormTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MaxNormVectorTests, ValueType, scai_numeric_test_types )
{
    IndexType n = 4;
    ValueType val = 5.0;

    auto vec = fillDenseVector<ValueType>( n, val );

    MaxNorm<ValueType> maxnorm;
    ValueType expected = val;
    BOOST_CHECK_EQUAL( expected, maxnorm( vec ) );
    BOOST_CHECK_EQUAL( vec.maxNorm(), maxnorm( vec ) );

    DenseVector<ValueType> vec2( HArray<ValueType>( { 1.0, -2.0, 3.0, -4.5 } ) );
    expected = 4.5;
    BOOST_CHECK_EQUAL( expected, maxnorm( vec2 ) );
    BOOST_CHECK_EQUAL( vec2.maxNorm(), maxnorm( vec2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MaxNormScalarTests, ValueType, scai_numeric_test_types )
{
    ValueType scalar( -4 );
    MaxNorm<ValueType> maxnorm;
    BOOST_CHECK_EQUAL( ValueType( 4 ), maxnorm( scalar ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
