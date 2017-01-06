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
#include <scai/utilskernel/LArray.hpp>

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
    const IndexType offsets[] =
    { 0, 1 };
    const ValueType values[] =
    { 0.5f, 0.5f, 0.3f, 0.2f, 0.1f, 0.0f };
    const IndexType numValues = sizeof( values ) / sizeof( ValueType );
    const IndexType numDiagonals = sizeof( offsets ) / sizeof( IndexType );
    LArray<IndexType> diaOffsets( 2, offsets );
    LArray<ValueType> diaValues( numValues, values );
    DIAStorage<ValueType> diaStorage( numRows, numColumns, numDiagonals, diaOffsets, diaValues );
    BOOST_REQUIRE_EQUAL( numRows, diaStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, diaStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numDiagonals, diaStorage.getNumDiagonals() );
    {
        ReadAccess<IndexType> diaOffsets( diaStorage.getOffsets() );
        ReadAccess<ValueType> diaValues( diaStorage.getValues() );

        // DIA keeps values in same order

        for ( IndexType i = 0; i < numDiagonals; ++i )
        {
            BOOST_CHECK_EQUAL( offsets[i], diaOffsets[i] );
        }

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( values[i], diaValues[i] );
        }
    }
    // copy constructor (TODO Bea: on all available locations, after host runs without mistake
    // DIAStorage<ValueType> diaStorageCopy( diaStorage, loc );
//    TODO ThoBra: weiter unten die "Nachbildung" des Copy-Constructors laeuft fehlerfrei,
//    hier knallt es beim letzten BOOST_CHECK_EQUAL
//    DIAStorage<ValueType> diaStorageCopy( diaStorage );
//    SCAI_LOG_INFO( logger, "copy constructor" );
//
//    BOOST_REQUIRE_EQUAL( numRows, diaStorageCopy.getNumRows() );
//    BOOST_REQUIRE_EQUAL( numColumns, diaStorageCopy.getNumColumns() );
//    BOOST_REQUIRE_EQUAL( numDiagonals, diaStorageCopy.getNumDiagonals() );
//
//    {
//        ReadAccess<ValueType> diaValues( diaStorageCopy.getValues() );
//        BOOST_CHECK_EQUAL( diaValues.size(), numValues );
//
//        BOOST_CHECK_EQUAL( values[0], diaValues[0] );
//        BOOST_CHECK_EQUAL( values[1], diaValues[1] );
//        BOOST_CHECK_EQUAL( values[2], diaValues[2] );
//        BOOST_CHECK_EQUAL( values[3], diaValues[3] );
//        BOOST_CHECK_EQUAL( values[4], diaValues[4] );
//        BOOST_CHECK_EQUAL( values[5], diaValues[5] );
//    }
    LArray<IndexType> csrIa;
    LArray<IndexType> csrJa;
    LArray<ValueType> csrValues;
    // const IndexType csrIaResult[] = { 0, 2, 4, 5 };
    // const IndexType csrJaResult[] = { 0, 1, 1, 2, 2 };
    // const ValueType csrValuesResult[] = { 0.5, 0.2, 0.5, 0.1, 0.3 };
    diaStorage.buildCSRData( csrIa, csrJa, csrValues );
    IndexType numRowsDia = diaStorage.getNumRows();
    IndexType numColumnsDia = diaStorage.getNumColumns();
    IndexType numValuesCSR = csrJa.size();
    diaStorage.setCSRData( numRowsDia, numColumnsDia, numValuesCSR, csrIa, csrJa, csrValues );
    {
        ReadAccess<ValueType> diaValues( diaStorage.getValues() );
        SCAI_LOG_INFO( logger, "diaValues.size() = " << diaValues.size() );
        BOOST_CHECK_EQUAL( diaValues.size(), numValues );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( values[0], diaValues[0] );
            BOOST_CHECK_EQUAL( values[1], diaValues[1] );
            BOOST_CHECK_EQUAL( values[2], diaValues[2] );
            BOOST_CHECK_EQUAL( values[3], diaValues[3] );
            BOOST_CHECK_EQUAL( values[4], diaValues[4] );
            BOOST_CHECK_EQUAL( values[5], diaValues[5] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, scai_numeric_test_types )
{
    // use template storage test
    storageSwapTest<DIAStorage<ValueType> >();
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
