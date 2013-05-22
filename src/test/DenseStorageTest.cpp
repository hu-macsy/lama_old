/**
 * @file DenseStorageTest.cpp
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
 * @brief Contains the implementation of the class DenseStorageTest
 * @author Alexander Büchel
 * @date 12.03.2012
 * $Id$
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/storage/DenseStorage.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseStorageTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.DenseStorageTest" );

typedef boost::mpl::list<float,double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, test_types ) {
    typedef T ValueType;

    DenseStorage<ValueType> denseStorage;
    MatrixStorageTest<ValueType> storageTest( denseStorage );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in DenseStorageTest." );
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

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    DenseStorage<double> denseStoraged;
    std::string s = denseStoraged.typeName();
    BOOST_CHECK_EQUAL( s, "DenseStorage<double>" );

    DenseStorage<float> denseStoragef;
    s = denseStoragef.typeName();
    BOOST_CHECK_EQUAL( s, "DenseStorage<float>" );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setZeroTest, T, test_types ) {
    typedef T ValueType;

    const IndexType numRows = 4;
    const IndexType numColumns = 4;

    static ValueType values[] =
    {   6.0, 0.0, 0.0, 4.0,
        7.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 9.0, 4.0,
        2.0, 5.0, 0.0, 3.0
    };

    ValueType eps = static_cast<ValueType>( 1E-5 );

    DenseStorage<double> denseStorage;
    denseStorage.setRawDenseData( numRows, numColumns, values, eps );

    denseStorage.setZero();

    for ( int i = 0; i < denseStorage.getNumRows(); ++i )
    {
        for ( int j = 0; j < denseStorage.getNumColumns(); ++j )
        {
            BOOST_CHECK_EQUAL( denseStorage.getValue( i, j ), 0.0 );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
