/**
 * @file ELLStorageTest.cpp
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
 * @brief Test cases for ELLStorage( only specific ones )
 * @author Thomas Brandes
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/hmemo/HArrayRef.hpp>
#include <scai/lama/storage/ELLStorage.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <scai/lama/test/storage/StorageTemplateTests.hpp>

using namespace scai;
using namespace lama;
using namespace utilskernel;
using namespace hmemo;

using common::Exception;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( ELLStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.ELLStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
{
    ContextPtr context = Context::getContextPtr();
    const IndexType numRows = 10;
    const IndexType numColumns = 15;
    ValueType zero = 0;

    ELLStorage<ValueType> ellStorage;
    ellStorage.setContextPtr( context );
    ellStorage.allocate( numRows, numColumns );

    BOOST_REQUIRE_EQUAL( numRows, ellStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, ellStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( IndexType( 0 ), ellStorage.getNumValues() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            ValueType val = ellStorage.getValue( i, j );
            BOOST_CHECK_EQUAL( zero, val );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructor1Test, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();
    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType ia[] =
    { 1, 1, 2 };
    const IndexType numValuesPerRow = 2;
    // Note: ja, values are stored column-major order
    const IndexType ja[] =
    { 0, 1, 2, 0, 0, 2 };
    const ValueType values[] =
    { 0.5f, 0.5f, 0.3f, 0.0f, 0.0f, 0.2f };
    const IndexType numValues = numRows * numValuesPerRow;
    const IndexType sizeJA = sizeof( ja ) / sizeof( IndexType );
    const IndexType sizeValues = sizeof( values ) / sizeof( ValueType );
    BOOST_CHECK_EQUAL( numValues, sizeJA );
    BOOST_CHECK_EQUAL( numValues, sizeValues );
    HArray<IndexType> ellIA( numRows, ia );
    HArray<IndexType> ellJA( numValues, ja );
    HArray<ValueType> ellValues( numValues, values );
    ELLStorage<ValueType> ellStorage( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues, loc );
    BOOST_REQUIRE_EQUAL( numRows, ellStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, ellStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValuesPerRow, ellStorage.getNumValuesPerRow() );
    BOOST_REQUIRE_EQUAL( ia[0] + ia[1] + ia[2], ellStorage.getNumValues() );
    {
        ReadAccess<IndexType> ellIA( ellStorage.getIA() );
        ReadAccess<IndexType> ellJA( ellStorage.getJA() );
        ReadAccess<ValueType> ellValues( ellStorage.getValues() );

        // ELL keeps values in same order

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_REQUIRE_EQUAL( ia[i], ellIA[i] );

            for ( IndexType jj = 0; jj < ia[i]; ++jj )
            {
                IndexType pos = i + jj * numRows;
                BOOST_CHECK_EQUAL( ja[pos], ellJA[pos] );
                BOOST_CHECK_EQUAL( values[pos], ellValues[pos] );
            }

            // values must have been filled up with 0 outside legal part

            for ( IndexType jj = ia[i]; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = i + jj * numRows;
                BOOST_CHECK_EQUAL( static_cast<ValueType>( 0 ), ellValues[pos] );
            }
        }
    }

    // copy constructor 

    ELLStorage<ValueType> ellStorageCopy( ellStorage );
    BOOST_REQUIRE_EQUAL( numRows, ellStorageCopy.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, ellStorageCopy.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValuesPerRow, ellStorageCopy.getNumValuesPerRow() );
    BOOST_REQUIRE_EQUAL( ia[0] + ia[1] + ia[2], ellStorageCopy.getNumValues() );
    {
        ReadAccess<IndexType> ellIALocal( ellStorageCopy.getIA() );
        ReadAccess<IndexType> ellJALocal( ellStorageCopy.getJA() );
        ReadAccess<ValueType> ellValuesLocal( ellStorageCopy.getValues() );

        // ELL keeps values in same order
        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_REQUIRE_EQUAL( ia[i], ellIALocal[i] );

            for ( IndexType jj = 0; jj < ia[i]; ++jj )
            {
                IndexType pos = i + jj * numRows;
                BOOST_CHECK_EQUAL( ja[pos], ellJALocal[pos] );
                BOOST_CHECK_EQUAL( values[pos], ellValuesLocal[pos] );
            }

            // values must have been filled up with 0 outside legal part

            for ( IndexType jj = ia[i]; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = i + jj * numRows;
                BOOST_CHECK_EQUAL( static_cast<ValueType>( 0 ), ellValuesLocal[pos] );
            }
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
    const IndexType numValuesPerRow = 2;
    const IndexType numValues = 4;

    HArray<IndexType> ellIA( { 1, 1, 2 } );
    HArray<IndexType> ellJA( { 0, 1, 2, 0, 0, 1 } );
    HArray<ValueType> ellValues( { 5, 5, 3, 0, 0, 2 } );

    const IndexType* ptrIA = getPointer( ellIA, context );
    const IndexType* ptrJA = getPointer( ellJA, context );
    const ValueType* ptrValues = getPointer( ellValues, context );

    SCAI_LOG_INFO( logger, "call moveConstructor with ell arrays" )

    ELLStorage<ValueType> ellStorage( numRows, numColumns, numValuesPerRow, std::move( ellIA ), ellJA, std::move( ellValues ) );

    BOOST_REQUIRE_EQUAL( numRows, ellStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, ellStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, ellStorage.getNumValues() );

    BOOST_CHECK( ellStorage.hasDiagonalProperty() );

    // verify that move was okay

    BOOST_CHECK_EQUAL( ptrIA, getPointer( ellStorage.getIA(), context ) );
    BOOST_CHECK( ptrJA != getPointer( ellStorage.getJA(), context ) );
    ptrJA = getPointer( ellStorage.getJA(), context );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( ellStorage.getValues(), context ) );

    BOOST_CHECK_EQUAL( ellIA.size(), 0 );
    BOOST_CHECK_EQUAL( ellJA.size(), numValuesPerRow * numRows );
    BOOST_CHECK_EQUAL( ellValues.size(), 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( splitUpTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType numValuesPerRow = 2;

    HArray<IndexType> ellIA( { 1, 1, 2 } );
    HArray<IndexType> ellJA( { 0, 1, 2, 0, 0, 1 } );
    HArray<ValueType> ellValues( { 5, 5, 3, 0, 0, 2 } );

    const IndexType* ptrIA = getPointer( ellIA, context );
    const IndexType* ptrJA = getPointer( ellJA, context );
    const ValueType* ptrValues = getPointer( ellValues, context );

    ELLStorage<ValueType> ellStorage( numRows, numColumns, numValuesPerRow,
                                      std::move( ellIA ), std::move( ellJA ), std::move( ellValues ) );

    IndexType outNumRows;
    IndexType outNumColumns;
    IndexType outNumValuesPerRow;
    HArray<IndexType> outIA;
    HArray<IndexType> outJA;
    HArray<ValueType> outValues;

    ellStorage.splitUp( outNumRows, outNumColumns, outNumValuesPerRow, outIA, outJA, outValues );

    BOOST_CHECK_EQUAL( outNumRows, numRows );
    BOOST_CHECK_EQUAL( outNumColumns, numColumns );
    BOOST_CHECK_EQUAL( outNumValuesPerRow, numValuesPerRow );

    BOOST_CHECK_EQUAL( ptrIA, getPointer( outIA, context ) );
    BOOST_CHECK_EQUAL( ptrJA, getPointer( outJA, context ) );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( outValues, context ) );

    BOOST_CHECK_EQUAL( ellStorage.getIA().size(), 0 );
    BOOST_CHECK_EQUAL( ellStorage.getJA().size(), 0 );
    BOOST_CHECK_EQUAL( ellStorage.getValues().size(), 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( checkTest, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    // This routine tests the check method of ELLStorage, individually for this class
    for ( int icase = 0; icase < 3; ++icase )
    {
        // build up a correct ELLPACK storage
        const IndexType numRows = 3;
        const IndexType numColumns = 3;
        const IndexType numValuesPerRow = 2;
        const IndexType numValues = numRows * numValuesPerRow;
        const IndexType ia[] = { 1, 1, 2 };
        const IndexType ja[] = { 0, 1, 2, 0, 0, 2 };
        // just make sure that ia and ja have correct sizes
        BOOST_REQUIRE_EQUAL( numRows, static_cast<IndexType>( sizeof( ia ) / sizeof( IndexType ) ) );
        BOOST_REQUIRE_EQUAL( numValues, static_cast<IndexType>( sizeof( ja ) / sizeof( IndexType ) ) );
        HArrayRef<IndexType> ellIA( numRows, ia );
        HArrayRef<IndexType> ellJA( numValues, ja );
        HArray<ValueType> ellValues( numValues, 1.0 ); // values needed, but do not matter here
        ELLStorage<ValueType> ellStorage;
        ellStorage.setContextPtr( loc );
        ellStorage.setELLData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

        if ( icase == 0 )
        {
            ellStorage.check( "test with correct values" );
        }
        else if ( icase == 1 )
        {
            //  -> invalid ia     { 1, 1, 3 }
            HArray<IndexType>& ellIA = const_cast<HArray<IndexType>&>( ellStorage.getIA() );
            HArrayUtils::setVal<IndexType>( ellIA, 2, 3 );
            BOOST_CHECK_THROW( { ellStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
        else if ( icase == 2 )
        {
            //  -> invalid ja     { 0, 1, 2, 0, 0, 2 }
            HArray<IndexType>& ellJA = const_cast<HArray<IndexType>&>( ellStorage.getJA() );
            HArrayUtils::setVal<IndexType>( ellJA, 5, 15 );
            BOOST_CHECK_THROW( { ellStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_numeric_test_types )
{
    SCAI_LOG_INFO( logger, "typeNameTest for ELLStorage<" << common::TypeTraits<ValueType>::id() << ">" )
    storageTypeNameTest<ELLStorage<ValueType> >( "ELL" );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( ELLCopyTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    copyStorageTest<ELLStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( compressTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType numValuesPerRow = 2;
    const IndexType realValues = 4;
    const IndexType storedValues = numRows * numValuesPerRow;

    /*   matrix:       saved values   saved columns   ia

          1  -   -     1 0            0 0             1
          -  1   -     1 0            1 0             1
          -  0   1     0 1            1 2             2
     */

    const IndexType ia[]     = { 1, 1, 2 };
    const IndexType ja[]     = { 0, 1, 1, 0, 0, 2 };
    const ValueType values[] = { 1, 1, 0, 0, 0, 1 };

    HArray<IndexType> ellIA( numRows, ia );
    HArray<IndexType> ellJA( storedValues, ja );
    HArray<ValueType> ellValues( storedValues, values );

    ELLStorage<ValueType> ell( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );
    ell.setContextPtr( context );

    BOOST_CHECK_EQUAL( realValues, ell.getNumValues() );

    ell.compress();
    // one zero element (not diagonal) is removed by compress
    BOOST_CHECK_EQUAL( realValues - 1, ell.getNumValues() );

    const IndexType expected_ia[]     = { 1, 1, 1 };
    const IndexType expected_ja[]     = { 0, 1, 2 };
    const IndexType expected_values[] = { 1, 1, 1 };

    HArray<IndexType> compressedIA     = ell.getIA();
    HArray<IndexType> compressedJA     = ell.getJA();
    HArray<ValueType> compressedValues = ell.getValues();

    ReadAccess<IndexType> rIA ( compressedIA );
    ReadAccess<IndexType> rJA ( compressedJA );
    ReadAccess<ValueType> rValues ( compressedValues );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        BOOST_CHECK_EQUAL( rIA[i], expected_ia[i] );
    }

    for ( IndexType i = 0; i < realValues - 1; ++i )
    {
        BOOST_CHECK_EQUAL( rJA[i], expected_ja[i] );
        BOOST_CHECK_EQUAL( rValues[i], expected_values[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
