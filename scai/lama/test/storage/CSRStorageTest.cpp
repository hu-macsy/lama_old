/**
 * @file CSRStorageTest.cpp
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
 * @brief Test cases for CSRStorage( only specific ones )
 * @author Thomas Brandes
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/lama/test/storage/TestStorages.hpp>
#include <scai/lama/test/storage/StorageTemplateTests.hpp>

using namespace scai;
using namespace lama;
using namespace utilskernel;
using namespace hmemo;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( CSRStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRStorageTest" );

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_arithmetic_test_types )
{
    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType ia[] =
    { 0, 1, 2, 4 };
    // Note: ja, values are stored column-major order
    const IndexType ja[] =
    { 0, 1, 2, 2 };
    const ValueType values[] =
    { 0.5f, 0.5f, 0.3f, 0.2f };
    const IndexType numValues = ia[numRows];
    const IndexType sizeJA = sizeof( ja ) / sizeof( IndexType );
    const IndexType sizeValues = sizeof( values ) / sizeof( ValueType );
    BOOST_CHECK_EQUAL( numValues, sizeJA );
    BOOST_CHECK_EQUAL( numValues, sizeValues );
    LArray<IndexType> csrIA( numRows + 1, ia );
    LArray<IndexType> csrJA( numValues, ja );
    LArray<ValueType> csrValues( numValues, values );
    CSRStorage<ValueType> csrStorage( numRows, numColumns, numValues, csrIA, csrJA, csrValues );
    BOOST_REQUIRE_EQUAL( numRows, csrStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, csrStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, csrStorage.getNumValues() );
    BOOST_CHECK( csrStorage.hasDiagonalProperty() );
    {
        ReadAccess<IndexType> csrIA( csrStorage.getIA() );
        ReadAccess<IndexType> csrJA( csrStorage.getJA() );
        ReadAccess<ValueType> csrValues( csrStorage.getValues() );

        // CSR keeps values in same order

        for ( IndexType i = 0; i < numRows + 1; ++i )
        {
            BOOST_CHECK_EQUAL( ia[i], csrIA[i] );
        }

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( ja[i], csrJA[i] );
            BOOST_CHECK_EQUAL( values[i], csrValues[i] );
        }
    }
    // copy constructor on all available locations
    CSRStorage<ValueType> csrStorageCopy( csrStorage, context );
    BOOST_REQUIRE_EQUAL( numRows, csrStorageCopy.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, csrStorageCopy.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, csrStorageCopy.getNumValues() );
    BOOST_CHECK( csrStorageCopy.hasDiagonalProperty() );
    {
        ReadAccess<IndexType> csrIA( csrStorageCopy.getIA() );
        ReadAccess<IndexType> csrJA( csrStorageCopy.getJA() );
        ReadAccess<ValueType> csrValues( csrStorageCopy.getValues() );

        // CSR keeps values in same order

        for ( IndexType i = 0; i < numRows + 1; ++i )
        {
            BOOST_CHECK_EQUAL( ia[i], csrIA[i] );
        }

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( ja[i], csrJA[i] );
            BOOST_CHECK_EQUAL( values[i], csrValues[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( compressTest, ValueType, scai_arithmetic_test_types )
{
    return; // TODO: fails

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType ia[] = { 0, 1, 2, 4 };
    const IndexType ja[] = { 0, 1, 1, 2 };
    const ValueType values[] = { 1, 1, 0, 1 };
    const IndexType numValues = ia[numRows];
    const IndexType sizeJA     = sizeof( ja ) / sizeof( IndexType );
    const IndexType sizeValues = sizeof( values ) / sizeof( ValueType );

    BOOST_CHECK_EQUAL( numValues, sizeJA );
    BOOST_CHECK_EQUAL( numValues, sizeValues );

    LArray<IndexType> csrIA( numRows + 1, ia );
    LArray<IndexType> csrJA( numValues, ja );
    LArray<ValueType> csrValues( numValues, values );

    CSRStorage<ValueType> csr( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    csr.setContextPtr( context );

    BOOST_CHECK_EQUAL( numValues, csr.getNumValues() );

    csr.compress();

    // one zero element (not diagonal) is removed by compress

    BOOST_CHECK_EQUAL( numValues - 1, csr.getNumValues() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, scai_arithmetic_test_types )
{
    SCAI_LOG_INFO( logger, "swapTest for CSRStorage<" << common::TypeTraits<ValueType>::id() << ">" )

    // use template storage test

    storageSwapTest<CSRStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_arithmetic_test_types )
{
    SCAI_LOG_INFO( logger, "typeNameTest for CSRStorage<" << common::TypeTraits<ValueType>::id() << ">" )

    storageTypeNameTest<CSRStorage<ValueType> >( "CSR" );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( sortRowTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 4;
    const IndexType numColumns = 4;

    const IndexType ia[] = { 0, 2, 4, 6, 8 };
    const IndexType ja[] = { 1, 0, 2, 1, 3, 2, 0, 3 };
    const IndexType sorted_ja[] = { 0, 1, 1, 2, 2, 3, 0, 3 };
    const IndexType values[] = { 1, 0, 2, 1, 3, 2, 0, 3 };

    const IndexType numValues = ia[numRows];

    LArray<IndexType> csrIA( numRows + 1, ia, context );
    LArray<IndexType> csrJA( numValues, ja, context );
    LArray<ValueType> csrValues( numValues, values, context );

    CSRStorage<ValueType> csrStorage;

    csrStorage.setContextPtr( context );
    csrStorage.allocate( numRows, numColumns );

    csrStorage.swap( csrIA, csrJA, csrValues );

    BOOST_CHECK_EQUAL( 0, csrJA.size() );
    BOOST_CHECK_EQUAL( 0, csrValues.size() );
    BOOST_CHECK_EQUAL( numValues, csrStorage.getNumValues() );

    bool diagonalProperty = csrStorage.hasDiagonalProperty();

    csrStorage.sortRows( diagonalProperty );

    const LArray<IndexType>& sortedJA = csrStorage.getJA();
    const LArray<ValueType>& sortedVals = csrStorage.getValues();

    BOOST_REQUIRE_EQUAL( numValues, sortedJA.size() );
    BOOST_REQUIRE_EQUAL( numValues, sortedVals.size() );

    for ( IndexType i = 0; i < numValues; ++i )
    {
        BOOST_CHECK_EQUAL( sorted_ja[i], sortedJA[i] );
        BOOST_CHECK_EQUAL( ValueType( sorted_ja[i] ), sortedVals[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( CSRCopyTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here

    copyStorageTest<CSRStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
