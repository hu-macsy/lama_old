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
 * @brief Test cases for CSRStorage( only specific ones )
 * @author Thomas Brandes
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/lama/test/storage/TestStorages.hpp>

using namespace scai;
using namespace lama;
using namespace utilskernel;
using namespace hmemo;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( CSRStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRStorageTest" );

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_arithmetic_test_types )
{
    ContextPtr context = Context::getContextPtr();

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
    LArray<IndexType> csrIA( numRows + 1, ia );
    LArray<IndexType> csrJA( numValues, ja );
    LArray<ValueType> csrValues( numValues, values );
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
    CSRStorage<ValueType> csrStorageCopy( csrStorage, context );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( compressTest, ValueType, scai_arithmetic_test_types )
{
    ContextPtr context = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType ia[] = { 0, 1, 2, 4 };
    const IndexType ja[] = { 0, 1, 1, 2 };
    const ValueType values[] = { 1, 1, 0, 1 };
    const IndexType numValues = ia[numRows];
    const IndexType sizeJA     = sizeof( ja ) / sizeof( IndexType );
    const IndexType sizeValues = sizeof( values ) / sizeof( ValueType );

    BOOST_CHECK_EQUAL( numValues, sizeJA );
    BOOST_CHECK_EQUAL( numValues, sizeValues );

    LArray<IndexType> csrIA( numRows + 1, ia );
    LArray<IndexType> csrJA( numValues, ja );
    LArray<ValueType> csrValues( numValues, values );

    CSRStorage<ValueType> csr( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    csr.setContextPtr( context );

    BOOST_CHECK_EQUAL( numValues, csr.getNumValues() );

    csr.compress();

    // one zero element (not diagonal) is removed by compress

    BOOST_CHECK_EQUAL( numValues - 1, csr.getNumValues() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, scai_arithmetic_test_types )
{
    ContextPtr context = Context::getContextPtr();

    CSRStorage<ValueType> csr1;
    CSRStorage<ValueType> csr2;

    setDenseData( csr1 );

    IndexType n = csr1.getNumRows();
    IndexType m = csr1.getNumColumns();

    LArray<ValueType> x( context );
    LArray<ValueType> y( context );
    LArray<ValueType> z1( context );
    LArray<ValueType> z2( context );
    
    ValueType alpha = 1.3;
    ValueType beta  = -0.5;

    HArrayUtils::setRandom( x, m, 1.0f );
    HArrayUtils::setRandom( y, n, 1.0f );
  
    csr1.matrixTimesVector( z1, alpha, x, beta, y );

    csr1.swap( csr2 );

    BOOST_CHECK_EQUAL( n, csr2.getNumRows() );
    BOOST_CHECK_EQUAL( m, csr2.getNumColumns() );

    BOOST_CHECK_EQUAL( 0, csr1.getNumRows() );
    BOOST_CHECK_EQUAL( 0, csr1.getNumColumns() );

    csr2.matrixTimesVector( z2, alpha, x, beta, y );

    BOOST_CHECK_EQUAL( 0, z1.maxDiffNorm( z2 ) );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_arithmetic_test_types )
{
    SCAI_LOG_INFO( logger, "typeNameTest for CSRStorage<" << common::TypeTraits<ValueType>::id() << ">" )

    // context does not matter here, so runs for every context

    std::string s = CSRStorage<ValueType>::typeName();

    BOOST_CHECK( s.length() > 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
