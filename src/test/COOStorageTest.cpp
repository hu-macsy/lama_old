/**
 * @file COOStorageTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/storage/COOStorage.hpp>
#include <lama/HostReadAccess.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( COOStorageTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.COOStorageTest" );

typedef boost::mpl::list<float,double> ValueTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, ValueTypes ) {
//TODO: Undo to Template-Testcase
    typedef T ValueType;

    COOStorage<ValueType> cooStorage;
    MatrixStorageTest<ValueType> storageTest( cooStorage );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in COOStorageTest." );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, T, ValueTypes ) {
    typedef T ValueType;

    const IndexType numRows = 10;
    const IndexType numColumns = 15;

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

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest1, T, ValueTypes ) {
    typedef T ValueType;

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
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    COOStorage<double> cooStoraged;
    std::string s = cooStoraged.typeName();
    BOOST_CHECK_EQUAL( s, "COOStorage<double>" );

    COOStorage<float> cooStoragef;
    s = cooStoragef.typeName();
    BOOST_CHECK_EQUAL( s, "COOStorage<float>" );
}
/* ------------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
