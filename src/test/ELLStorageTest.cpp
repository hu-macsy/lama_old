/**
 * @file ELLStorageTest.cpp
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
 * @brief Contains the implementation of the class ELLStorageTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

#include <lama/storage/ELLStorage.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ELLStorageTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.ELLStorageTest" );

typedef boost::mpl::list<float,double> ValueTypes;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, ValueTypes ) {
    typedef T ValueType;

    ELLStorage<ValueType> ellStorage;
    MatrixStorageTest<ValueType> storageTest( ellStorage );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in ELLStorageTest." );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, T, ValueTypes ) {
    typedef T ValueType;

    const IndexType numRows = 10;
    const IndexType numColumns = 15;

    ELLStorage<ValueType> ellStorage( numRows, numColumns );

// TODO: to check CUDA:
// ELLStorage<T> ellStorage( numRows, numColumns, Context::CUDA );

    BOOST_REQUIRE_EQUAL( numRows, ellStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, ellStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( 0, ellStorage.getNumValues() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float> ( ellStorage.getValue( i, j ) );
            BOOST_CHECK_SMALL( v , 1.0e-5f );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest1, T, ValueTypes ) {
    typedef T ValueType;

    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    const IndexType ia[] =
    {   1, 1, 2};

    const IndexType numValuesPerRow = 2;

// Note: ja, values are stored column-major order

    const IndexType ja[] =
    {   0, 1, 2, 0, 0, 2};
    const ValueType values[] =
    {   0.5f, 0.5f, 0.3f, 0.0f, 0.0f, 0.2f};

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
        HostReadAccess<IndexType> ellIA( ellStorage.getIA() );
        HostReadAccess<IndexType> ellJA( ellStorage.getJA() );
        HostReadAccess<ValueType> ellValues( ellStorage.getValues() );

        // ELL keeps values in same order

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( ia[i], ellIA[i] );
        }

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( ja[i], ellJA[i] );
            BOOST_CHECK_EQUAL( values[i], ellValues[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    ELLStorage<double> ellStoraged;
    std::string s = ellStoraged.typeName();
    BOOST_CHECK_EQUAL( s, "ELLStorage<double>" );

    ELLStorage<float> ellStoragef;
    s = ellStoragef.typeName();
    BOOST_CHECK_EQUAL( s, "ELLStorage<float>" );
}
/* ------------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
