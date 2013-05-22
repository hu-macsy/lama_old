/**
 * @file SparseAssemblyStorageTest.cpp
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
 * @brief Contains the implementation of the class SparseAssemblyStorageTest
 * @author Alexander Büchel
 * @date 12.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/storage/CSRStorage.hpp>
#include <lama/storage/ELLStorage.hpp>
#include <lama/storage/JDSStorage.hpp>
#include <lama/storage/DIAStorage.hpp>
#include <lama/storage/COOStorage.hpp>
#include <lama/storage/DenseStorage.hpp>

#include <test/TestSparseMatrices.hpp>
#include <test/MatrixStorageTest.hpp>

#include <lama/storage/SparseAssemblyStorage.hpp>

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SparseAssemblyStorageTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.SparseAssemblyStorageTest" );

typedef boost::mpl::list<float,double> ValueTypes;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixStorageTestWithSparseAssemblyStorage, T, ValueTypes ) {
    SparseAssemblyStorage<T> assemblyStorage;
    MatrixStorageTest<T> storageTest( assemblyStorage );
    storageTest.runTests();
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ConstructorTest, T, ValueTypes ) {
    const IndexType numRows = 10;
    const IndexType numColumns = 15;

    SparseAssemblyStorage<T> assemblyStorage( numRows, numColumns );

    BOOST_REQUIRE_EQUAL( numRows, assemblyStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, assemblyStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( 0, assemblyStorage.getNumValues() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float> ( assemblyStorage.getValue( i, j ) );
            BOOST_CHECK_SMALL( v, 1.0e-5f );
        }
    }

// fix diagonal property for each row, can be done in parallel

    #pragma omp parallel for
    for ( IndexType i = 0; i < numRows; ++i )
    {
        assemblyStorage.fixDiagonalProperty( i );
    }

// diagonal property must be recomputed

    assemblyStorage.resetDiagonalProperty();

    BOOST_CHECK( assemblyStorage.hasDiagonalProperty() );
    BOOST_CHECK_EQUAL( numRows, assemblyStorage.getNumValues() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetRowTest, T, ValueTypes ) {
    const IndexType n = 10;

    CSRStorage<T> csrStorage;
    csrStorage.setIdentity( n );

    SparseAssemblyStorage<T> assemblyStorage( n, n );

    #pragma omp parallel for
    for ( IndexType i = 0; i < n; ++i )
    {
        // Note: this test verifies also thread-safe use of LAMA arrays

        LAMAArray<IndexType> ja;
        LAMAArray<T> values;

        {
            HostWriteOnlyAccess<IndexType> wJa( ja, 1 );
            HostWriteOnlyAccess<T> wValues( values, 1 );
            wJa[0] = i;
            wValues[0] = 1.0;
        }

        assemblyStorage.setRow( i, ja, values );
    }

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < n; ++j )
        {
            BOOST_CHECK_EQUAL( csrStorage.getValue( i, j ), assemblyStorage.getValue( i, j ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    SparseAssemblyStorage<double> aStoraged;
    std::string s = aStoraged.typeName();
    BOOST_CHECK_EQUAL( s, "SparseAssemblyStorage<double>" );

    SparseAssemblyStorage<float> aStoragef;
    s = aStoragef.typeName();
    BOOST_CHECK_EQUAL( s, "SparseAssemblyStorage<float>" );
}
/* ------------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
