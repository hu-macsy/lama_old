/**
 * @file CSRStorageTest.cpp
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
 * @brief Contains the implementation of the class CSRStorageTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/storage/CSRStorage.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CSRStorageTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.CSRStorageTest" );

typedef boost::mpl::list<float,double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, test_types )
{
    typedef T ValueType;

    CSRStorage<ValueType> csrStorage;
    MatrixStorageTest<ValueType> storageTest( csrStorage );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in CSRStorageTest." );
        MATRIXSTORAGE_COMMONTESTCASES( storageTest );
    }
    else
    {
        CONTEXTLOOP()
        {
            GETCONTEXT( context );
            storageTest.mMatrixStorage.setContext( context );
            LAMA_LOG_INFO( logger, "Using context = " << storageTest.mMatrixStorage.getContext().getType() );
            storageTest.runTests();
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, T, test_types )
{
    const IndexType numRows = 10;
    const IndexType numColumns = 15;

    CSRStorage<T> csrStorage;

    csrStorage.allocate( numRows, numColumns );

    BOOST_REQUIRE_EQUAL( numRows, csrStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, csrStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( 0, csrStorage.getNumValues() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float> ( csrStorage.getValue( i, j ) );
            BOOST_CHECK_SMALL( v , 1.0e-5f );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest1, ValueType, test_types )
{
    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    const IndexType ia[] = { 0, 1, 2, 4 };

// Note: ja, values are stored column-major order

    const IndexType ja[] = { 0, 1, 2, 2 };
    const ValueType values[] = { 0.5f, 0.5f, 0.3f, 0.2f };

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
        HostReadAccess<IndexType> csrIA( csrStorage.getIA() );
        HostReadAccess<IndexType> csrJA( csrStorage.getJA() );
        HostReadAccess<ValueType> csrValues( csrStorage.getValues() );

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

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    CSRStorage<double> csrStoraged;
    std::string s = csrStoraged.typeName();
    BOOST_CHECK_EQUAL( s, "CSRStorage<double>" );

    CSRStorage<float> csrStoragef;
    s = csrStoragef.typeName();
    BOOST_CHECK_EQUAL( s, "CSRStorage<float>" );
}
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
