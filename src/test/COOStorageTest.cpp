/**
 * @file COOStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Contains the implementation of the class COOStorageTest
 * @author Alexander Büchel
 * @date 12.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/storage/COOStorage.hpp>
#include <lama/HostReadAccess.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>
#include <lama/LAMAArrayUtils.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

namespace lama
{
namespace COOStorageTest
{

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{

    COOStorage<ValueType> cooStorage;
    MatrixStorageTest<ValueType> storageTest( cooStorage );

    storageTest.mMatrixStorage.setContext( loc );

    if( base_test_case )
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

    // constructor doesn't exist with other location
    COOStorage<ValueType> cooStorage( numRows, numColumns );

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

template<typename ValueType>
void constructorTest1( ContextPtr loc )
{

    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    const IndexType ia[] =
    {   0, 1, 2, 1};
    const IndexType ja[] =
    {   0, 1, 2, 2};
    const ValueType values[] =
    {   0.5f, 0.5f, 0.3f, 0.2f};

    const IndexType numValues = sizeof( values ) / sizeof( ValueType );

    LAMAArray<IndexType> cooIA( numValues, ia );
    LAMAArray<IndexType> cooJA( numValues, ja );
    LAMAArray<ValueType> cooValues( numValues, values );

    COOStorage<ValueType> cooStorage( numRows, numColumns, cooIA, cooJA, cooValues );

    BOOST_REQUIRE_EQUAL( numRows, cooStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, cooStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numValues, cooStorage.getNumValues() );

    {
        HostReadAccess<IndexType> cooIA( cooStorage.getIA() );
        HostReadAccess<IndexType> cooJA( cooStorage.getJA() );
        HostReadAccess<ValueType> cooValues( cooStorage.getValues() );

        // COO keeps values in same order
        for ( IndexType i = 0; i < numValues; ++i)
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
        HostReadAccess<IndexType> cooIA( cooStorageCopy.getIA() );
        HostReadAccess<IndexType> cooJA( cooStorageCopy.getJA() );
        HostReadAccess<ValueType> cooValues( cooStorageCopy.getValues() );

        // COO keeps values in same order
        for ( IndexType i = 0; i < numValues; ++i)
        {
            BOOST_CHECK_EQUAL( ia[i], cooIA[i] );
            BOOST_CHECK_EQUAL( ja[i], cooJA[i] );
            BOOST_CHECK_EQUAL( values[i], cooValues[i] );
        }
    }

}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void checkTest( ContextPtr loc )
{

    // This routine tests the check method of COOStorage, individually for this class
    for( int icase = 0; icase < 3; ++icase )
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

        LAMAArrayRef<IndexType> cooIA( ia, numValues );
        LAMAArrayRef<IndexType> cooJA( ja, numValues );
        LAMAArray<ValueType> cooValues( numValues, 1.0 ); // values needed, but do not matter here

        COOStorage<ValueType> cooStorage;

        cooStorage.setContext( loc );

        cooStorage.setCOOData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

        if( icase == 0 )
        {
            cooStorage.check( "test with correct values" );
        }
        else if( icase == 1 )
        {
            //  -> invalid ia     { 4, 1, 1, 2, 2, 3 }

            LAMAArray<IndexType>& cooIA = const_cast<LAMAArray<IndexType>&>( cooStorage.getIA() );
            LAMAArrayUtils::setVal( cooIA, 0, numRows );
            BOOST_CHECK_THROW( { cooStorage.check( "Expect illegal index in IA" ); }, Exception );
        }
        else if( icase == 2 )
        {
            //  -> invalid ja     { 0, 0, 1, 1, 4, 3 }

            LAMAArray<IndexType>& cooJA = const_cast<LAMAArray<IndexType>&>( cooStorage.getJA() );
            LAMAArrayUtils::setVal( cooJA, 4, numColumns );
            BOOST_CHECK_THROW( { cooStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void typeNameTest()
{
    COOStorage<ValueType> cooStorage;
    std::string s = cooStorage.typeName();

    BOOST_CHECK( s.length() > 0 );
}

} // namespace COOStorageTest
} // namespace lama

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( COOStorageTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.COOStorageTest" )
LAMA_AUTO_TEST_CASE_CT( checkTest, COOStorageTest )
LAMA_AUTO_TEST_CASE_CT( commonTestCases, COOStorageTest )
LAMA_AUTO_TEST_CASE_T( constructorTest, COOStorageTest )
LAMA_AUTO_TEST_CASE_CT( constructorTest1, COOStorageTest )
LAMA_AUTO_TEST_CASE_T( typeNameTest, COOStorageTest )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
