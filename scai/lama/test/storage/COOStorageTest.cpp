/**
 * @file COOStorageTest.cpp
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
 * @brief Test cases for COOStorage( only specific ones )
 * @author Thomas Brandes
 * @date 12.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/COOStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/lama/test/storage/StorageTemplateTests.hpp>

using namespace scai;
using namespace utilskernel;
using namespace hmemo;
using namespace lama;

using scai::common::Exception;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( COOStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.COOStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();
    const IndexType numRows = 10;
    const IndexType numColumns = 15;
    // constructor doesn't exist with other location
    COOStorage<ValueType> cooStorage( numRows, numColumns );
    cooStorage.prefetch( loc );
    BOOST_REQUIRE_EQUAL( numRows, cooStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, cooStorage.getNumColumns() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float> ( cooStorage.getValue( i, j ) );
            BOOST_CHECK_SMALL( v, 1.0e-5f );
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
    {   0, 1, 2, 1};
    const IndexType ja[] =
    {   0, 1, 2, 2};
    const ValueType values[] =
    {   0.5f, 0.5f, 0.3f, 0.2f};
    const IndexType numValues = sizeof( values ) / sizeof( ValueType );
    LArray<IndexType> cooIA( numValues, ia );
    LArray<IndexType> cooJA( numValues, ja );
    LArray<ValueType> cooValues( numValues, values );
    COOStorage<ValueType> cooStorage( numRows, numColumns, cooIA, cooJA, cooValues );
    BOOST_REQUIRE_EQUAL( numRows, cooStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, cooStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, cooStorage.getNumValues() );
    {
        ReadAccess<IndexType> cooIA( cooStorage.getIA() );
        ReadAccess<IndexType> cooJA( cooStorage.getJA() );
        ReadAccess<ValueType> cooValues( cooStorage.getValues() );

        // COO keeps values in same order
        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( ia[i], cooIA[i] );
            BOOST_CHECK_EQUAL( ja[i], cooJA[i] );
            BOOST_CHECK_EQUAL( values[i], cooValues[i] );
        }
    }
    // copy constructor doesn't exist with other location
    COOStorage<ValueType> cooStorageCopy( cooStorage, loc );
    BOOST_REQUIRE_EQUAL( numRows, cooStorageCopy.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, cooStorageCopy.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, cooStorageCopy.getNumValues() );
    {
        ReadAccess<IndexType> cooIA( cooStorageCopy.getIA() );
        ReadAccess<IndexType> cooJA( cooStorageCopy.getJA() );
        ReadAccess<ValueType> cooValues( cooStorageCopy.getValues() );

        // COO keeps values in same order
        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( ia[i], cooIA[i] );
            BOOST_CHECK_EQUAL( ja[i], cooJA[i] );
            BOOST_CHECK_EQUAL( values[i], cooValues[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( checkTest, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();

    // This routine tests the check method of COOStorage, individually for this class
    for ( int icase = 0; icase < 3; ++icase )
    {
        // build up a correct COOStorage
        const IndexType numRows = 4;
        const IndexType numColumns = 4;
        const IndexType numValues = 6;
        const IndexType ia[] =
        { 0, 1, 1, 2, 2, 3 };
        const IndexType ja[] =
        { 0, 0, 1, 1, 2, 3 };
        // just make sure that ia and ja have correct sizes
        BOOST_REQUIRE_EQUAL( numValues, static_cast<IndexType>( sizeof( ia ) / sizeof( IndexType ) ) );
        BOOST_REQUIRE_EQUAL( numValues, static_cast<IndexType>( sizeof( ja ) / sizeof( IndexType ) ) );
        HArrayRef<IndexType> cooIA( numValues, ia );
        HArrayRef<IndexType> cooJA( numValues, ja );
        HArray<ValueType> cooValues( numValues, 1.0 ); // values needed, but do not matter here
        COOStorage<ValueType> cooStorage;
        cooStorage.setContextPtr( loc );
        cooStorage.setCOOData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

        if ( icase == 0 )
        {
            cooStorage.check( "test with correct values" );
        }
        else if ( icase == 1 )
        {
            //  -> invalid ia     { 4, 1, 1, 2, 2, 3 }
            HArray<IndexType>& cooIA = const_cast<HArray<IndexType>&>( cooStorage.getIA() );
            HArrayUtils::setVal( cooIA, 0, numRows );
            BOOST_CHECK_THROW( { cooStorage.check( "Expect illegal index in IA" ); }, Exception );
        }
        else if ( icase == 2 )
        {
            //  -> invalid ja     { 0, 0, 1, 1, 4, 3 }
            HArray<IndexType>& cooJA = const_cast<HArray<IndexType>&>( cooStorage.getJA() );
            HArrayUtils::setVal( cooJA, 4, numColumns );
            BOOST_CHECK_THROW( { cooStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, scai_numeric_test_types )
{
    // use template storage test
    storageSwapTest<COOStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_numeric_test_types )
{
    SCAI_LOG_INFO( logger, "typeNameTest for COOStorage<" << common::TypeTraits<ValueType>::id() << ">" )
    storageTypeNameTest<COOStorage<ValueType> >( "COO" );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( COOCopyTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    copyStorageTest<COOStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
