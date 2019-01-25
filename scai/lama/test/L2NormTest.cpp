/**
 * @file L2NormTest.cpp
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
 * @brief Contains specific tests for derived norm class L2Norm
 * @author Thomas Brandes
 * @date 03.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( L2NormTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.L2NormTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( L2NormVectorTests, ValueType, scai_numeric_test_types )
{
    IndexType n = 4;
    ValueType val = 5.0;
    auto vec = denseVectorFill<ValueType>( n, val );
    L2Norm<ValueType> l2norm;
    ValueType expected = common::Math::sqrt( ValueType( n ) * val * val );
    BOOST_CHECK_EQUAL( expected, l2norm( vec ) );
    BOOST_CHECK_EQUAL( vec.l2Norm(), l2norm( vec ) );
    WriteAccess<ValueType> hwa( vec.getLocalValues() );
    hwa[0] = 1.0;
    hwa[1] = -2.0;
    hwa[2] = 3.0;
    hwa[3] = -4.0;
    hwa.release();
    expected = static_cast<ValueType>( 5.47722 );
    ValueType s = l2norm( vec );
    RealType<ValueType> eps = 0.001; // 
    BOOST_CHECK( common::Math::abs( expected - s ) < eps );
    BOOST_CHECK_EQUAL( vec.l2Norm(), l2norm( vec ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( L2NormScalarTests, ValueType, scai_numeric_test_types )
{
    ValueType scalar( -4 );
    L2Norm<ValueType> l2norm;
    BOOST_CHECK_EQUAL( ValueType( 4 ), l2norm( scalar ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
