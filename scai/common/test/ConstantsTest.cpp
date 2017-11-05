/**
 * @file test/ConstantsTest.cpp
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
 * @brief Test routines for class Complex
 * @author Eric Schricker
 * @date 24.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Constants.hpp>

using namespace scai::common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ConstantsTest )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( zeroTest, ValueType, scai_numeric_test_types )
{
    ValueType zero( 0 );
    BOOST_CHECK( zero == Constants::ZERO );
    BOOST_CHECK( Constants::ZERO == zero );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( oneTest, ValueType, scai_numeric_test_types )
{
    ValueType one( 1 );
    BOOST_CHECK( one == Constants::ONE );
    BOOST_CHECK( Constants::ONE == one );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( randomTest, ValueType, scai_numeric_test_types )
{
    ValueType val = Math::random<ValueType>( 1000 );   // random val in 0..1000
    BOOST_CHECK( val - val == Constants::ZERO );
    BOOST_CHECK( Constants::ZERO == val - val );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
