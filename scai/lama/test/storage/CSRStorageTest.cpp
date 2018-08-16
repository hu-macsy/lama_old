/**
 * @file CSRStorageTest.cpp
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

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
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
    HArray<IndexType> csrIA( numRows + 1, ia );
    HArray<IndexType> csrJA( numValues, ja );
    HArray<ValueType> csrValues( numValues, values );
    CSRStorage<ValueType> csrStorage( numRows, numColumns, csrIA, csrJA, csrValues, context );
    BOOST_REQUIRE_EQUAL( numRows, csrStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, csrStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, csrStorage.getNumValues() );
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

    // copy constructor

    CSRStorage<ValueType> csrStorageCopy( csrStorage );
    BOOST_REQUIRE_EQUAL( numRows, csrStorageCopy.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, csrStorageCopy.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, csrStorageCopy.getNumValues() );
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

template<typename ValueType>
static inline const ValueType* getPointer( const HArray<ValueType>& a, ContextPtr ctx )
{
    ReadAccess<ValueType> rA( a, ctx );
    return rA.get();
}

BOOST_AUTO_TEST_CASE( moveConstructorTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType numValues  = 4;

    HArray<IndexType> csrIA( { 0, 1, 2, 4 } );
    HArray<IndexType> csrJA( { 0, 1, 2, 2 } );
    HArray<ValueType> csrValues( { 5, 5, 3, 2 } );

    const IndexType* ptrIA = getPointer( csrIA, context );
    const IndexType* ptrJA = getPointer( csrJA, context );
    const ValueType* ptrValues = getPointer( csrValues, context );

    SCAI_LOG_INFO( logger, "call moveConstructor with csr arrays" )

    CSRStorage<ValueType> csrStorage( numRows, numColumns, std::move( csrIA ), csrJA, std::move( csrValues ) );

    BOOST_REQUIRE_EQUAL( numRows, csrStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, csrStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, csrStorage.getNumValues() );

    // verify that move was okay

    BOOST_CHECK_EQUAL( ptrIA, getPointer( csrStorage.getIA(), context ) );
    BOOST_CHECK( ptrJA != getPointer( csrStorage.getJA(), context ) );
    ptrJA = getPointer( csrStorage.getJA(), context );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( csrStorage.getValues(), context ) );

    BOOST_CHECK_EQUAL( csrIA.size(), 0 );
    BOOST_CHECK_EQUAL( csrJA.size(), numValues );
    BOOST_CHECK_EQUAL( csrValues.size(), 0 );

    SCAI_LOG_INFO( logger, "call moveConstructor of CSRStorage" )

    CSRStorage<ValueType> csrStorage1( std::move( csrStorage ) );

    // verify that move was okay

    BOOST_CHECK_EQUAL( ptrIA, getPointer( csrStorage1.getIA(), context ) );
    BOOST_CHECK_EQUAL( ptrJA, getPointer( csrStorage1.getJA(), context ) );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( csrStorage1.getValues(), context ) );

    BOOST_CHECK_EQUAL( csrStorage.getIA().size(), 0 );
    BOOST_CHECK_EQUAL( csrStorage.getJA().size(), 0 );
    BOOST_CHECK_EQUAL( csrStorage.getValues().size(), 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( moveAssignTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    HArray<IndexType> csrIA( { 0, 1, 2, 4 } );
    HArray<IndexType> csrJA( { 0, 1, 2, 2 } );
    HArray<ValueType> csrValues( { 5, 5, 3, 2 } );

    const IndexType numValues = csrJA.size();

    CSRStorage<ValueType> csrStorage( numRows, numColumns, csrIA, csrJA, csrValues );

    const IndexType* ptrIA = getPointer( csrStorage.getIA(), context );
    const IndexType* ptrJA = getPointer( csrStorage.getJA(), context );
    const ValueType* ptrValues = getPointer( csrStorage.getValues(), context );

    CSRStorage<ValueType> csrStorage1;
    csrStorage1 = std::move( csrStorage );

    BOOST_CHECK_EQUAL( csrStorage1.getNumRows(), numRows );
    BOOST_CHECK_EQUAL( csrStorage1.getNumColumns(), numColumns );
    BOOST_CHECK_EQUAL( csrStorage1.getNumValues(), numValues );

    // verify that move was okay

    BOOST_CHECK_EQUAL( ptrIA, getPointer( csrStorage1.getIA(), context ) );
    BOOST_CHECK_EQUAL( ptrJA, getPointer( csrStorage1.getJA(), context ) );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( csrStorage1.getValues(), context ) );

    BOOST_CHECK_EQUAL( csrStorage.getIA().size(), 0 );
    BOOST_CHECK_EQUAL( csrStorage.getJA().size(), 0 );
    BOOST_CHECK_EQUAL( csrStorage.getValues().size(), 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( splitUpTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows   = 3;
    const IndexType numColumns = 3;

    HArray<IndexType> csrIA( { 0, 1, 2, 4 } );
    HArray<IndexType> csrJA( { 0, 1, 2, 2 } );
    HArray<ValueType> csrValues( { 5, 5, 3, 2 } );

    const IndexType* ptrIA = getPointer( csrIA, context );
    const IndexType* ptrJA = getPointer( csrJA, context );
    const ValueType* ptrValues = getPointer( csrValues, context );

    CSRStorage<ValueType> csrStorage( numRows, numColumns, std::move( csrIA ), std::move( csrJA ), std::move( csrValues ) );

    IndexType outNumRows;
    IndexType outNumColumns;
    HArray<IndexType> outIA;
    HArray<IndexType> outJA;
    HArray<ValueType> outValues;

    csrStorage.splitUp( outNumRows, outNumColumns, outIA, outJA, outValues );

    BOOST_CHECK_EQUAL( outNumRows, numRows );
    BOOST_CHECK_EQUAL( outNumColumns, numColumns );

    BOOST_CHECK_EQUAL( ptrIA, getPointer( outIA, context ) );
    BOOST_CHECK_EQUAL( ptrJA, getPointer( outJA, context ) );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( outValues, context ) );

    BOOST_CHECK_EQUAL( csrStorage.getIA().size(), 0 );
    BOOST_CHECK_EQUAL( csrStorage.getJA().size(), 0 );
    BOOST_CHECK_EQUAL( csrStorage.getValues().size(), 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( compressTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    /*   matrix:

          1  -   -
          -  1   -
          -  0   1
     */

    const IndexType ia[]     = { 0, 1, 2,    4 };
    const IndexType ja[]     = { 0, 1, 1, 2 };
    const ValueType values[] = { 1, 1, 0, 1 };

    const IndexType numValues = ia[numRows];
    const IndexType sizeJA     = sizeof( ja ) / sizeof( IndexType );
    const IndexType sizeValues = sizeof( values ) / sizeof( ValueType );
    BOOST_CHECK_EQUAL( numValues, sizeJA );
    BOOST_CHECK_EQUAL( numValues, sizeValues );
    HArray<IndexType> csrIA( numRows + 1, ia );
    HArray<IndexType> csrJA( numValues, ja );
    HArray<ValueType> csrValues( numValues, values );
    CSRStorage<ValueType> csr;
    csr.setCSRData( numRows, numColumns, csrIA, csrJA, csrValues );
    csr.setContextPtr( context );
    BOOST_CHECK_EQUAL( numValues, csr.getNumValues() );
    csr.compress();
    // one zero element (not diagonal) is removed by compress
    BOOST_CHECK_EQUAL( numValues - 1, csr.getNumValues() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_numeric_test_types )
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
    const IndexType numValues = ia[numRows];

    HArray<IndexType> csrIA( numRows + 1, ia, context );
    HArray<IndexType> csrJA( numValues, ja, context );
    HArray<ValueType> csrValues;
    HArrayUtils::assign<ValueType, IndexType>( csrValues, csrJA, context );

    CSRStorage<ValueType> csrStorage;
    csrStorage.setContextPtr( context );
    csrStorage.allocate( numRows, numColumns );

    csrStorage.swap( csrIA, csrJA, csrValues );

    BOOST_CHECK_EQUAL( IndexType( 0 ), csrJA.size() );
    BOOST_CHECK_EQUAL( IndexType( 0 ), csrValues.size() );
    BOOST_CHECK_EQUAL( numValues, csrStorage.getNumValues() );

    csrStorage.sortRows( );

    const HArray<IndexType>& sortedJA = csrStorage.getJA();
    const HArray<ValueType>& sortedVals = csrStorage.getValues();

    BOOST_REQUIRE_EQUAL( numValues, sortedJA.size() );
    BOOST_REQUIRE_EQUAL( numValues, sortedVals.size() );

    for ( IndexType i = 0; i < numValues; ++i )
    {
        BOOST_CHECK_EQUAL( sorted_ja[i], sortedJA[i] );
        BOOST_CHECK_EQUAL( ValueType( sorted_ja[i] ), sortedVals[i] );
    }
}

BOOST_AUTO_TEST_CASE( CSRCopyTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    copyStorageTest<CSRStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
