/**
 * @file ELLStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief Test cases for ELLStorage( only specific ones )
 * @author Thomas Brandes
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LArray.hpp>
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

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_arithmetic_test_types )
{
    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 10;
    const IndexType numColumns = 15;

    ValueType zero = 0;

    ELLStorage<ValueType> ellStorage( numRows, numColumns, context );

    BOOST_REQUIRE_EQUAL( numRows, ellStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, ellStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( 0, ellStorage.getNumValues() );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( constructor1Test, ValueType, scai_arithmetic_test_types )
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
    LArray<IndexType> ellIA( numRows, ia );
    LArray<IndexType> ellJA( numValues, ja );
    LArray<ValueType> ellValues( numValues, values );
    ELLStorage<ValueType> ellStorage( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );
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
    // copy constructor on all available locations
    ELLStorage<ValueType> ellStorageCopy( ellStorage, loc );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( checkTest, ValueType, scai_arithmetic_test_types )
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
        const IndexType ia[] =
        { 1, 1, 2 };
        const IndexType ja[] =
        { 0, 1, 2, 0, 0, 2 };
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
            LArray<IndexType>& ellIA = const_cast<LArray<IndexType>&>( ellStorage.getIA() );
            HArrayUtils::setVal( ellIA, 2, 3 );
            BOOST_CHECK_THROW( { ellStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
        else if ( icase == 2 )
        {
            //  -> invalid ja     { 0, 1, 2, 0, 0, 2 }
            LArray<IndexType>& ellJA = const_cast<LArray<IndexType>&>( ellStorage.getJA() );
            HArrayUtils::setVal( ellJA, 5, 15 );
            BOOST_CHECK_THROW( { ellStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, scai_arithmetic_test_types )
{
    // use template storage test

    storageSwapTest<ELLStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_SUITE_END();
