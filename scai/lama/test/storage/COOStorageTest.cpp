/**
 * @file COOStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Test cases for COOStorage( only specific ones )
 * @author Thomas Brandes
 * @date 12.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/COOStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LArray.hpp>

using namespace scai;
using namespace utilskernel;
using namespace hmemo;
using namespace lama;

using scai::common::Exception;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( COOStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.COOStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( constructor1Test, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( checkTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_arithmetic_test_types )
{
    SCAI_LOG_INFO( logger, "typeNameTest for COOStorage<" << common::TypeTraits<ValueType>::id() << ">" )

    // context does not matter here, so runs for every context

    std::string s = COOStorage<ValueType>::typeName();

    BOOST_CHECK( s.length() > 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
