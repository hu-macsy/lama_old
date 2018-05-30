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
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/utilskernel.hpp>

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

BOOST_AUTO_TEST_CASE_TEMPLATE( constructor1Test, ValueType, scai_numeric_test_types )
{
    ContextPtr loc = Context::getContextPtr();
    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    HArray<IndexType> cooIA( {  0, 1, 2, 1 } );
    HArray<IndexType> cooJA( {  0, 1, 2, 2 } );
    HArray<ValueType> cooValues( { 0.5, 0.5, 0.3, 0.2 } );

    COOStorage<ValueType> cooStorage( numRows, numColumns, cooIA, cooJA, cooValues, loc );

    // COO keeps values in same order

    BOOST_TEST( hostReadAccess( cooIA ), cooStorage.getIA() );
    BOOST_TEST( hostReadAccess( cooJA ), cooStorage.getJA() );
    BOOST_TEST( hostReadAccess( cooValues ), cooStorage.getValues() );

    // copy constructor 

    COOStorage<ValueType> cooStorageCopy( cooStorage );

    BOOST_TEST( hostReadAccess( cooIA ), cooStorageCopy.getIA() );
    BOOST_TEST( hostReadAccess( cooJA ), cooStorageCopy.getJA() );
    BOOST_TEST( hostReadAccess( cooValues ), cooStorageCopy.getValues() );
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
    const IndexType numValues = 4;

    HArray<IndexType> cooIA( { 0, 1, 2, 2 } );
    HArray<IndexType> cooJA( { 0, 1, 1, 2 } );
    HArray<ValueType> cooValues( { 5, 5, 3, 3 } );

    const IndexType* ptrIA = getPointer( cooIA, context );
    const IndexType* ptrJA = getPointer( cooJA, context );
    const ValueType* ptrValues = getPointer( cooValues, context );

    SCAI_LOG_INFO( logger, "call moveConstructor with coo arrays" )

    COOStorage<ValueType> cooStorage( numRows, numColumns, std::move( cooIA ), cooJA, std::move( cooValues ) );

    BOOST_REQUIRE_EQUAL( numRows, cooStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, cooStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, cooStorage.getNumValues() );

    // verify that move was okay (no changes as COO data is sorted), ja has been copied

    BOOST_CHECK_EQUAL( ptrIA, getPointer( cooStorage.getIA(), context ) );
    BOOST_CHECK( ptrJA != getPointer( cooStorage.getJA(), context ) );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( cooStorage.getValues(), context ) );

    BOOST_CHECK_EQUAL( cooIA.size(), 0 );
    BOOST_CHECK_EQUAL( cooJA.size(), numValues );
    BOOST_CHECK_EQUAL( cooValues.size(), 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( splitUpTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows   = 3;
    const IndexType numColumns = 3;

    HArray<IndexType> cooIA( { 0, 1, 2, 2 } );
    HArray<IndexType> cooJA( { 0, 1, 1, 2 } );
    HArray<ValueType> cooValues( { 5, 5, 3, 3 } );

    const IndexType* ptrIA = getPointer( cooIA, context );
    const IndexType* ptrJA = getPointer( cooJA, context );
    const ValueType* ptrValues = getPointer( cooValues, context );

    COOStorage<ValueType> cooStorage( numRows, numColumns, std::move( cooIA ), std::move( cooJA ), std::move( cooValues ) );

    IndexType outNumRows;
    IndexType outNumColumns;
    HArray<IndexType> outIA;
    HArray<IndexType> outJA;
    HArray<ValueType> outValues;

    cooStorage.splitUp( outNumRows, outNumColumns, outIA, outJA, outValues );

    BOOST_CHECK_EQUAL( outNumRows, numRows );
    BOOST_CHECK_EQUAL( outNumColumns, numColumns );

    BOOST_CHECK_EQUAL( ptrIA, getPointer( outIA, context ) );
    BOOST_CHECK_EQUAL( ptrJA, getPointer( outJA, context ) );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( outValues, context ) );

    BOOST_CHECK_EQUAL( cooStorage.getIA().size(), 0 );
    BOOST_CHECK_EQUAL( cooStorage.getJA().size(), 0 );
    BOOST_CHECK_EQUAL( cooStorage.getValues().size(), 0 );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( matrixPlusMatrixTest, ValueType, scai_numeric_test_types )
{
    // Input data       COO1         COO2        2 * COO1 + 3 * COO2 
    //
    //                 1  2  -      -  11  -      2  37  -
    //                 -  3  4      -   -  12     -  6   44 
    //                 -  -  5      -  13  -      - 39   10

    HArray<IndexType> ia1( { 0, 0, 1, 1, 2 } );
    HArray<IndexType> ja1( { 0, 1, 1, 2, 2 } );
    HArray<ValueType> values1( { 1, 2, 3, 4, 5 } );

    HArray<IndexType> ia2( { 0, 1, 2 } );
    HArray<IndexType> ja2( { 1, 2, 1 } );
    HArray<ValueType> values2( { 11, 12, 13 } );

    // Expected result, will be sorted due to building unique values

    HArray<IndexType> ia( { 0, 0, 1, 1, 2, 2 } );
    HArray<IndexType> ja( { 0, 1, 1, 2, 1, 2 } );
    HArray<ValueType> values( { 2, 37, 6, 44, 39, 10 } );

    IndexType numRows = 3;
    IndexType numColumns = 4;

    COOStorage<ValueType> coo1( numRows, numColumns, ia1, ja1, values1 );
    COOStorage<ValueType> coo2( numRows, numColumns, ia2, ja2, values2 );

    COOStorage<ValueType> coo;
    coo.matrixPlusMatrix( 2, coo1, 3, coo2 );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( coo.getIA() ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( coo.getJA() ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( coo.getValues() ), boost::test_tools::per_element() );

    // now a bit more tricky, but with same result

    coo = COOStorage<ValueType>( numRows, numColumns, ia1, ja1, values1 );   
    auto csr2 = convert<CSRStorage<ValueType>>( coo2 );

    coo.matrixPlusMatrix( 3, csr2, 2, coo );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( coo.getIA() ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( coo.getJA() ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( coo.getValues() ), boost::test_tools::per_element() );
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
