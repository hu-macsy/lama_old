/**
 * @file COOUtilsTest.cpp
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
 * @brief Test of methods for COOUtils
 * @author Thomas Brandes
 * @date 12.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/COOUtils.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

using namespace scai;
using namespace utilskernel;
using namespace hmemo;
using namespace lama;

using scai::common::Exception;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( COOUtilsTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.COOUtilsTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( sortTest, ValueType, scai_numeric_test_types )
{
    // Note: we sort coordinates and not values, so this test works also for complex numbers

    // Here is the unsorted COO data

    HArray<IndexType> ia(     { 2, 1, 0, 2, 1, 1, 0, 0, 0 } );
    HArray<IndexType> ja(     { 2, 1, 0, 1, 1, 1, 0, 2, 1 } );
    HArray<ValueType> values( { 1, 2, 3, 4, 5, 6, 7, 8, 9 } );

    // This is the expected COO data

    HArray<IndexType> ia1(     { 0, 0, 0, 0, 1, 1, 1, 2, 2 } );
    HArray<IndexType> ja1(     { 0, 0, 1, 2, 1, 1, 1, 1, 2 } );
    HArray<ValueType> values1( { 3, 7, 9, 8, 2, 5, 6, 4, 1 } );

    COOUtils::sort( ia, ja, values );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( ia1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( ja1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( values1 ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( uniqueTest, ValueType, scai_numeric_test_types )
{
    // Here is the sorted COO data with double entries

    HArray<IndexType> ia(     { 0, 0, 0, 0, 1, 1, 1, 2, 2 } );
    HArray<IndexType> ja(     { 0, 0, 1, 2, 1, 1, 1, 1, 2 } );
    HArray<ValueType> values( { 3, 7, 9, 8, 2, 5, 6, 4, 1 } );

    // This is the unique COO data

    HArray<IndexType> ia1(     { 0, 0, 0, 1, 2, 2 } );
    HArray<IndexType> ja1(     { 0, 1, 2, 1, 1, 2 } );
    HArray<ValueType> values1( { 7, 9, 8, 6, 4, 1 } );

    COOUtils::unique( ia, ja, values );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( ia1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( ja1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( values1 ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( uniqueOpTest, ValueType, scai_numeric_test_types )
{
    // Here is the sorted COO data with double entries

    HArray<IndexType> ia(     { 0, 0, 0, 0, 1, 1, 1, 2, 2 } );
    HArray<IndexType> ja(     { 0, 0, 1, 2, 1, 1, 1, 1, 2 } );
    HArray<ValueType> values( { 3, 7, 9, 8, 2, 5, 6, 4, 1 } );

    // This is the unique COO data

    HArray<IndexType> ia1(     {  0, 0, 0,  1, 2, 2 } );
    HArray<IndexType> ja1(     {  0, 1, 2,  1, 1, 2 } );
    HArray<ValueType> values1( { 10, 9, 8, 13, 4, 1 } );

    COOUtils::unique( ia, ja, values, common::BinaryOp::ADD );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( ia1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( ja1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( values1 ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
