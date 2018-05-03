/**
 * @file DIAStorageTest.cpp
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
 * @brief Contains the implementation of the class DIAStorageTest
 * @author Thomas Brandes
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/utilskernel.hpp>

#include <scai/lama/test/storage/StorageTemplateTests.hpp>

using namespace scai;
using namespace lama;
using namespace utilskernel;
using namespace hmemo;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DIAStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.DIAStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
{
    // Test the full DIAStorge constructor and the individual getter routines of DIA storage

    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    // test storage =  ( ( 0.5, 0.2, 0.0 ), ( 0.0, 0.5, 0.1 ), ( 0.0, 0.0, 0.3 ) )

    HArray<IndexType> diaOffsets( { 0, 1 } );
    HArray<ValueType> diaValues( { 0.5, 0.5, 0.3, 0.2, 0.1, 0.0 } );

    // just verify correct sizes of the arrays

    SCAI_ASSERT_EQ_ERROR( diaValues.size(), numRows * diaOffsets.size(), "illegal DIA data" )

    DIAStorage<ValueType> diaStorage( numRows, numColumns, diaOffsets, diaValues );

    BOOST_REQUIRE_EQUAL( numRows, diaStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, diaStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( diaOffsets.size(), diaStorage.getNumDiagonals() );

	BOOST_TEST( hostReadAccess( diaValues ) == hostReadAccess( diaStorage.getValues() ), boost::test_tools::per_element() );

    HArray<IndexType> csrIa;
    HArray<IndexType> csrJa;
    HArray<ValueType> csrValues;

    HArray<IndexType> expectedIA( { 0, 2, 4, 5 } );
    HArray<IndexType> expectedJA( { 0, 1, 1, 2, 2 } );
    HArray<ValueType> expectedValues( { 0.5, 0.2, 0.5, 0.1, 0.3 } );

    diaStorage.buildCSRData( csrIa, csrJa, csrValues );

	BOOST_TEST( hostReadAccess( expectedIA ) == hostReadAccess( csrIa ), boost::test_tools::per_element() );
	BOOST_TEST( hostReadAccess( expectedJA ) == hostReadAccess( csrJa ), boost::test_tools::per_element() );
	BOOST_TEST( hostReadAccess( expectedValues ) == hostReadAccess( csrValues ), boost::test_tools::per_element() );

    DIAStorage<ValueType> diaStorage1;
    diaStorage1.setCSRData( numRows, numColumns, csrIa, csrJa, csrValues );

    auto diaRead = hostReadAccess( diaValues );
    auto diaRead1 = hostReadAccess( diaStorage1.getValues() );

	// BOOST_TEST( hostReadAccess( diaValues ) == hostReadAccess( diaStorage1.getValues() ), boost::test_tools::per_element() );

    BOOST_CHECK_EQUAL_COLLECTIONS( diaRead.begin(), diaRead.end(), diaRead1.begin(), diaRead1.end() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_numeric_test_types )
{
    SCAI_LOG_INFO( logger, "typeNameTest for DIAStorage<" << common::TypeTraits<ValueType>::id() << ">" )
    storageTypeNameTest<DIAStorage<ValueType> >( "DIA" );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( DIACopyTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    copyStorageTest<DIAStorage<ValueType> >();
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
