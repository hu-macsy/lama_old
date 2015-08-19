/**
 * @file CSRStorageTest.cpp
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
 * @brief Contains the implementation of the class CSRStorageTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/CSRStorage.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

extern bool base_test_case;
extern std::string testcase;

namespace scai
{
namespace lama
{
namespace CSRStorageTest
{

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{
    CSRStorage<ValueType> csrStorage;
    MatrixStorageTest<ValueType> storageTest( csrStorage );
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
void constructorTest()
{
    const IndexType numRows = 10;
    const IndexType numColumns = 15;
    CSRStorage<ValueType> csrStorage;
    csrStorage.allocate( numRows, numColumns );
    BOOST_REQUIRE_EQUAL( numRows, csrStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, csrStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( 0, csrStorage.getNumValues() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float>( csrStorage.getValue( i, j ) );
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
    LAMAArray<IndexType> csrIA( numRows + 1, ia );
    LAMAArray<IndexType> csrJA( numValues, ja );
    LAMAArray<ValueType> csrValues( numValues, values );
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
    CSRStorage<ValueType> csrStorageCopy( csrStorage, loc );
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

template<typename ValueType>
void typeNameTest()
{
    CSRStorage<ValueType> csrStorage;
    std::string s = csrStorage.typeName();
    BOOST_CHECK( s.length() > 0 );
}

} /* end namespace CSRStorageTest */

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( CSRStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRStorageTest" );
LAMA_AUTO_TEST_CASE_CT( commonTestCases, CSRStorageTest )
LAMA_AUTO_TEST_CASE_T( constructorTest, CSRStorageTest )
LAMA_AUTO_TEST_CASE_CT( constructorTest1, CSRStorageTest )
LAMA_AUTO_TEST_CASE_T( typeNameTest, CSRStorageTest )
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
