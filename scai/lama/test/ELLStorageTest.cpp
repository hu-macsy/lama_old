/**
 * @file ELLStorageTest.cpp
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
 * @brief Contains the implementation of the class ELLStorageTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <test/MatrixStorageTest.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <scai/lama/LAMAArrayUtils.hpp>
#include <scai/lama/storage/ELLStorage.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using scai::common::Exception;

extern bool base_test_case;
extern std::string testcase;

namespace scai
{
namespace lama
{
namespace ELLStorageTest
{

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{
    ELLStorage<ValueType> ellStorage;
    MatrixStorageTest<ValueType> storageTest( ellStorage );
    storageTest.mMatrixStorage.setContext( loc );

    if ( base_test_case )
    {
        MATRIXSTORAGE_COMMONTESTCASES( storageTest );
    }
    else
    {
        storageTest.runTests();
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void constructorTest( ContextPtr loc )
{
    const IndexType numRows = 10;
    const IndexType numColumns = 15;
    ELLStorage<ValueType> ellStorage( numRows, numColumns, loc->getType() );
    BOOST_REQUIRE_EQUAL( numRows, ellStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, ellStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( 0, ellStorage.getNumValues() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float>( ellStorage.getValue( i, j ) );
            BOOST_CHECK_SMALL( v, 1.0e-5f );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void constructorTest1( ContextPtr loc )
{
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
    LAMAArray<IndexType> ellIA( numRows, ia );
    LAMAArray<IndexType> ellJA( numValues, ja );
    LAMAArray<ValueType> ellValues( numValues, values );
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

template<typename ValueType>
void checkTest( ContextPtr loc )
{
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
        LAMAArrayRef<IndexType> ellIA( ia, numRows );
        LAMAArrayRef<IndexType> ellJA( ja, numValues );
        LAMAArray<ValueType> ellValues( numValues, 1.0 ); // values needed, but do not matter here
        ELLStorage<ValueType> ellStorage;
        ellStorage.setContext( loc );
        ellStorage.setELLData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

        if ( icase == 0 )
        {
            ellStorage.check( "test with correct values" );
        }
        else if ( icase == 1 )
        {
            //  -> invalid ia     { 1, 1, 3 }
            LAMAArray<IndexType>& ellIA = const_cast<LAMAArray<IndexType>&>( ellStorage.getIA() );
            LAMAArrayUtils::setVal( ellIA, 2, 3 );
            BOOST_CHECK_THROW( { ellStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
        else if ( icase == 2 )
        {
            //  -> invalid ja     { 0, 1, 2, 0, 0, 2 }
            LAMAArray<IndexType>& ellJA = const_cast<LAMAArray<IndexType>&>( ellStorage.getJA() );
            LAMAArrayUtils::setVal( ellJA, 5, 15 );
            BOOST_CHECK_THROW( { ellStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void typeNameTest()
{
    ELLStorage<ValueType> ellStorage;
    std::string s = ellStorage.typeName();
    BOOST_CHECK( s.length() > 0 );
}

} /* end namespace ELLStorageTest */

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( ELLStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.ELLStorageTest" )
LAMA_AUTO_TEST_CASE_CT( commonTestCases, ELLStorageTest, scai::lama )
LAMA_AUTO_TEST_CASE_CT( constructorTest, ELLStorageTest, scai::lama )
LAMA_AUTO_TEST_CASE_CT( constructorTest1, ELLStorageTest, scai::lama )
LAMA_AUTO_TEST_CASE_T( typeNameTest, ELLStorageTest )
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
