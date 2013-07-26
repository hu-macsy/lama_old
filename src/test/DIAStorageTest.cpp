/**
 * @file DIAStorageTest.cpp
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
 * @brief Contains the implementation of the class DIAStorageTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/ContextFactory.hpp>

#include <lama/storage/DIAStorage.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DIAStorageTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.DIAStorageTest" )

typedef boost::mpl::list<float,double> ValueTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, ValueTypes )
{
    typedef T ValueType;

    DIAStorage<ValueType> diaStorage;
    MatrixStorageTest<ValueType> storageTest( diaStorage );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in DIAStorageTest." );
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

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, T, ValueTypes )
{
    const IndexType numRows = 10;
    const IndexType numColumns = 15;

    DIAStorage<T> diaStorage( numRows, numColumns );

    BOOST_REQUIRE_EQUAL( numRows, diaStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, diaStorage.getNumColumns() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float> ( diaStorage.getValue( i, j ) );
            BOOST_CHECK_SMALL( v , 1.0e-5f );
        }
    }

    BOOST_CHECK_EQUAL( 0, diaStorage.getNumValues() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest1, ValueType, ValueTypes )
{
// Test the full DIAStorge constructor and the individual getter routines of DIA storage

    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    const IndexType offsets[] =
    {   0, 1};
    const ValueType values[] =
    {   0.5f, 0.5f, 0.3f, 0.2f, 0.1f, 0.0f};

    const IndexType numValues = sizeof( values ) / sizeof( ValueType );
    const IndexType numDiagonals = sizeof( offsets ) / sizeof( IndexType );

    LAMAArray<IndexType> diaOffsets( 2, offsets );
    LAMAArray<ValueType> diaValues( numValues, values );

    DIAStorage<ValueType> diaStorage( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    BOOST_REQUIRE_EQUAL( numRows, diaStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, diaStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numDiagonals, diaStorage.getNumDiagonals() );

    {
        HostReadAccess<IndexType> diaOffsets( diaStorage.getOffsets() );
        HostReadAccess<ValueType> diaValues( diaStorage.getValues() );

        // DIA keeps values in same order

        for ( IndexType i = 0; i < numDiagonals; ++i )
        {
            BOOST_CHECK_EQUAL( offsets[i], diaOffsets[i] );
        }

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( values[i], diaValues[i] );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    DIAStorage<double> diaStoraged;
    std::string s = diaStoraged.typeName();
    BOOST_CHECK_EQUAL( s, "DIAStorage<double>" );

    DIAStorage<float> diaStoragef;
    s = diaStoragef.typeName();
    BOOST_CHECK_EQUAL( s, "DIAStorage<float>" );
}
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
