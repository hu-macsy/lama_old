/**
 * @file JDSStorageTest.cpp
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
 * @brief Test cases for JDSStorage( specific ones and all of MatrixStorageTest )
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/storage/JDSStorage.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( JDSStorageTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.JDSStorageTest" );

typedef boost::mpl::list<float,double> ValueTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, ValueTypes ) {
    typedef T ValueType;

    JDSStorage<ValueType> jdsStorage;
    MatrixStorageTest<ValueType> storageTest( jdsStorage );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in JDSStorageTest." );
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

    JDSStorage<ValueType> jdsStorage( numRows, numColumns );

    BOOST_REQUIRE_EQUAL( numRows, jdsStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, jdsStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( 0, jdsStorage.getNumValues() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float> ( jdsStorage.getValue( i, j ) );
            BOOST_CHECK_SMALL( v , 1.0e-5f );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest1, T, ValueTypes ) {
    typedef T ValueType;

    const IndexType numRows = 3;
    const IndexType numColumns = 3;

// Note: ja, values are stored column-major order

    const IndexType ja[] =
    {   2, 0, 1, 2,};
    const ValueType values[] =
    {   0.3f, 0.5f, 0.5f, 0.2f};
    const IndexType perm[] =
    {   2, 0, 1};
    const IndexType ilg[] =
    {   2, 1, 1};
    const IndexType dlg[] =
    {   3, 1};

    const IndexType numDiagonals = sizeof( dlg ) / sizeof( IndexType );
    const IndexType numValues = sizeof( ja ) / sizeof( IndexType );
    const IndexType sizeValues = sizeof( values ) / sizeof( ValueType );
    const IndexType sizePerm = sizeof( perm ) / sizeof( IndexType );
    const IndexType sizeILG = sizeof( ilg ) / sizeof( IndexType );

    BOOST_CHECK_EQUAL( numValues, sizeValues );
    BOOST_CHECK_EQUAL( numRows, sizePerm );
    BOOST_CHECK_EQUAL( numRows, sizeILG );

    LAMAArray<IndexType> jdsILG( numRows, ilg );
    LAMAArray<IndexType> jdsDLG( numDiagonals, dlg );
    LAMAArray<IndexType> jdsPerm( numRows, perm );
    LAMAArray<IndexType> jdsJA( numValues, ja );
    LAMAArray<ValueType> jdsValues( numValues, values );

// Call the specific constructor for JDS storage

    JDSStorage<ValueType> jdsStorage( numRows, numColumns, numValues, numDiagonals,
                                      jdsDLG, jdsILG, jdsPerm, jdsJA, jdsValues );

// Test all the getter routines, includes the specific ones for JDSStorage

    BOOST_REQUIRE_EQUAL( numRows, jdsStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, jdsStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numDiagonals, jdsStorage.getNumDiagonals() );

    {
        HostReadAccess<IndexType> jdsILG( jdsStorage.getIlg() );
        HostReadAccess<IndexType> jdsPerm( jdsStorage.getPerm() );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( ilg[i], jdsILG[i] );
            BOOST_CHECK_EQUAL( perm[i], jdsPerm[i] );
        }

        HostReadAccess<IndexType> jdsDLG( jdsStorage.getDlg() );

        for ( IndexType i = 0; i < numDiagonals; ++i)
        {
            BOOST_CHECK_EQUAL( dlg[i], jdsDLG[i] );
        }

        HostReadAccess<IndexType> jdsJA( jdsStorage.getJA() );
        HostReadAccess<ValueType> jdsValues( jdsStorage.getValues() );

        // JDS keeps values in same order

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( ja[i], jdsJA[i] );
            BOOST_CHECK_EQUAL( values[i], jdsValues[i] );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    JDSStorage<double> jdsStoraged;
    std::string s = jdsStoraged.typeName();
    BOOST_CHECK_EQUAL( s, "JDSStorage<double>" );

    JDSStorage<float> jdsStoragef;
    s = jdsStoragef.typeName();
    BOOST_CHECK_EQUAL( s, "JDSStorage<float>" );
}
/* ------------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
