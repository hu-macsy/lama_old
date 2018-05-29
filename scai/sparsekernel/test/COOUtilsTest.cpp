/**
 * @file COOUtilsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Contains tests for the COOUtils interface to be tested on different devices
 * @author Thomas Brandes
 * @date 19.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/kregistry/KernelContextFunction.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>
#include <scai/sparsekernel/COOUtils.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>
#include <scai/sparsekernel/test/TestData1.hpp>
#include <scai/sparsekernel/test/TestData2.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;

using namespace hmemo;
using namespace kregistry;
using namespace sparsekernel;
using namespace utilskernel;
using common::TypeTraits;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( COOUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.COOUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( convert2CSRTest )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> cooIA( { 0, 0, 1, 1, 1, 2, 4, 4, 5, 6, 6 }, testContext );
    HArray<IndexType> expIA( { 0,    2,       5, 6, 6, 8, 9,    11 }  );

    const IndexType numRows = 7;

    SCAI_ASSERT_EQ_ERROR( expIA[numRows], cooIA.size(), "Wrong setup for test" )

    HArray<IndexType> csrIA;

    COOUtils::convertCOO2CSR( csrIA, cooIA, numRows, testContext );

    ContextPtr validContext = csrIA.getValidContext( testContext );

    if ( validContext->getType() != testContext->getType() ) 
    {
        SCAI_LOG_WARN( logger, "convertCSR was not executed on " << *testContext << " but on " << *validContext );
    }

    BOOST_TEST( hostReadAccess( expIA ) == hostReadAccess( csrIA ), per_element() );

    // now convert back the CSR offsets to COO ia indexes and compare against original values

    HArray<IndexType> cooIA1;   // new result array

    COOUtils::convertCSR2COO( cooIA1, csrIA, cooIA.size(), testContext );

    BOOST_TEST( hostReadAccess( cooIA1 ) == hostReadAccess( cooIA ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( offsets2iaTest )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "offsets2ia test for " << *testContext << " on " << *testContext )

    HArray<IndexType> csrIA( { 0, 2, 2, 3, 5 }, testContext );

    const IndexType numRows = csrIA.size() - 1;

    IndexType numValues = HArrayUtils::getVal( csrIA, numRows );

    HArray<IndexType> cooIA;  // result array

    COOUtils::convertCSR2COO( cooIA, csrIA, numValues, testContext );

    HArray<IndexType> expectedIA( { 0, 0, 2, 3, 3 } );

    BOOST_TEST( hostReadAccess( cooIA ) == hostReadAccess( expectedIA ), boost::test_tools::per_element() );

} // offsets2iaTest

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( hasDiagonalTest )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "has diagonal test for " << *testContext )

    HArray<IndexType> cooIA( { 0, 1, 2, 3 }, testContext );
    HArray<IndexType> cooJA( { 0, 1, 2, 3 }, testContext );

    BOOST_CHECK( COOUtils::hasDiagonal( 4, 4, cooIA, cooJA, testContext ) );
    BOOST_CHECK( !COOUtils::hasDiagonal( 5, 5, cooIA, cooJA, testContext ) );

    cooIA = { 0, 1, 2, 3 };
    cooJA = { 0, 1, 3, 3 };

    BOOST_CHECK( !COOUtils::hasDiagonal( 4, 4, cooIA, cooJA, testContext ) );

    cooIA.clear();
    cooJA.clear();

    BOOST_CHECK( COOUtils::hasDiagonal( 0, 0, cooIA, cooJA, testContext ) );
    BOOST_CHECK( !COOUtils::hasDiagonal( 1, 1, cooIA, cooJA, testContext ) );
} 

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> cooIA( testContext );
    HArray<IndexType> cooJA( testContext );
    HArray<ValueType> cooValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    HArray<ValueType> denseValues( testContext );

    data1::getDenseTestData( numRows, numColumns, denseValues );

    ValueType zero = 0;

    auto readDenseValues = hostReadAccess( denseValues );
    auto readSparseValues = hostReadAccess( cooValues );

    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            IndexType pos = COOUtils::getValuePos( i, j, cooIA, cooJA, testContext );

            ValueType denseVal = readDenseValues[ i * numColumns + j ];

            if ( pos == invalidIndex )
            {
                BOOST_CHECK_EQUAL( denseVal, zero );
            }
            else
            {
                BOOST_CHECK_EQUAL( denseVal, readSparseValues[pos] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getColumnPositionsTest )
{
    ContextPtr testContext = ContextFix::testContext;

    //    1.0   -   2.0
    //    0.5  0.3   -
    //     -    -   3.0

    HArray<IndexType> cooIA(  { 0, 0, 1, 1, 2 }, testContext );
    HArray<IndexType> cooJA(  { 0, 2, 0, 1, 2 }, testContext );

    HArray<IndexType> column1Pos( { 3 } );         // expected positions for column 1
    HArray<IndexType> column2Pos( { 1, 4 } );      // expected positions for column 2

    HArray<IndexType> pos;   // result for positions

    COOUtils::getColumnPositions( pos, cooJA, 1, testContext );
    BOOST_TEST( hostReadAccess( pos ) == hostReadAccess( column1Pos ), per_element() );

    COOUtils::getColumnPositions( pos, cooJA, 2, testContext );
    BOOST_TEST( hostReadAccess( pos ) == hostReadAccess( column2Pos ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getRowPositionsTest )
{
    ContextPtr testContext = ContextFix::testContext;

    //    1.0   -   2.0
    //    0.5  0.3   -
    //     -    -    -
    //     -    -   3.0

    HArray<IndexType> cooIA(  { 0, 0, 1, 1, 3 }, testContext ); 

    IndexType offset;
    IndexType n;

    COOUtils::getRowPositions( offset, n, cooIA, 1, testContext );
    BOOST_CHECK_EQUAL( offset, 2 );
    BOOST_CHECK_EQUAL( n, 2 );

    COOUtils::getRowPositions( offset, n, cooIA, 2, testContext );
    BOOST_CHECK_EQUAL( offset, invalidIndex );
    BOOST_CHECK_EQUAL( n, 0 );

    COOUtils::getRowPositions( offset, n, cooIA, 3, testContext );
    BOOST_CHECK_EQUAL( offset, 4 );
    BOOST_CHECK_EQUAL( n, 1 );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "scaleRows test for " << *testContext )

    HArray<IndexType> cooIA( testContext );
    HArray<IndexType> cooJA( testContext );
    HArray<ValueType> cooValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    HArray<ValueType> savedValues( cooValues );  // keep a copy for comparison later

    const ValueType row_factors[]   = { 2, 3, 4, 5, 1, 3, 2 };

    const IndexType n_factors = sizeof( row_factors ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( numRows, n_factors );

    HArray<ValueType> rows( n_factors, row_factors, testContext );

    COOUtils::scaleRows( cooValues, cooIA, rows, testContext );

    // prove by hand on host

    {
        auto rIA = hostReadAccess( cooIA );
        auto rRows = hostReadAccess( rows );
        auto rValues = hostReadAccess( cooValues );
        auto rSavedValues = hostReadAccess( savedValues );

        for ( IndexType k = 0; k < numValues; ++k )
        {
            ValueType f = rRows[ rIA[k] ];
            BOOST_CHECK_EQUAL( f * rSavedValues[k], rValues[k] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

// BOOST_AUTO_TEST_CASE_TEMPLATE( gemvNormalTest, ValueType, scai_numeric_test_types )
BOOST_AUTO_TEST_CASE( gemvNormalTest )
{
    typedef double ValueType;

    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> cooIA( testContext );
    HArray<IndexType> cooJA( testContext );
    HArray<ValueType> cooValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    SCAI_ASSERT_EQ_ERROR( cooIA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( cooJA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( cooValues.size(), numValues, "size mismatch" )

    const ValueType y_values[]   = { 1, -1, 2, -2, 1, 1, -1 };
    const ValueType x_values[]   = { 3, -3, 2, -2 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numColumns, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, n_y, "size mismatch" );

    HArray<ValueType> x( numColumns, x_values, testContext );
    HArray<ValueType> y( numRows, y_values, testContext );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        HArray<ValueType> res( testContext );

        SCAI_LOG_DEBUG( logger, "compute res = " << alpha << " * COO * x + " << beta << " * y "
                         << ", with x = " << x << ", y = " << y
                         << ", COO: ia = " << cooIA << ", ja = " << cooJA << ", values = " << cooValues )

        COOUtils::gemv( res, numRows, numColumns, cooIA, cooJA, cooValues, common::MatrixOp::NORMAL,
                        alpha, x, beta, y, testContext );
                   
        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( expectedRes ) == hostReadAccess( res ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTransposeTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> cooIA( testContext );
    HArray<IndexType> cooJA( testContext );
    HArray<ValueType> cooValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    SCAI_ASSERT_EQ_ERROR( cooIA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( cooJA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( cooValues.size(), numValues, "size mismatch" )

    HArray<ValueType> x( { 3, -2, -2, 3, 1, 0, 1 }, testContext );
    HArray<ValueType> y( { 1, -1, 2, -2 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numRows, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, y.size(), "size mismatch" );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        HArray<ValueType> res( testContext );

        COOUtils::gemv( res, numRows, numColumns, cooIA, cooJA, cooValues, common::MatrixOp::TRANSPOSE,
                        alpha, x, beta, y, testContext );

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext )

    HArray<IndexType> cooIA( testContext );
    HArray<IndexType> cooJA( testContext );
    HArray<ValueType> cooValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data2::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    // verify that there is at least one entry for each diagonal element 

    BOOST_CHECK( COOUtils::hasDiagonal( numRows, numColumns, cooIA, cooJA, testContext ) );

    const ValueType rhs_values[]   = { 1, -1, 2, -2 };
    const ValueType old_values[]   = { 3, -2, -2, 3 };

    HArray<ValueType> rhs( numRows, rhs_values, testContext );
    HArray<ValueType> oldSolution( numRows, old_values, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> res( testContext );

        COOUtils::jacobi( res, omega, oldSolution, rhs, cooIA, cooJA, cooValues, testContext );

        HArray<ValueType> expectedRes( testContext );

        data2::getJacobiResult( expectedRes, oldSolution, omega, rhs );

        BOOST_CHECK( HArrayUtils::maxDiffNorm( expectedRes, res ) < 0.001 );

        // ?? how to set tolerance 

        // BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( sortTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    // Note: we sort coordinates and not values, so this test works also for complex numbers

    // Here is the unsorted COO data

    HArray<IndexType> ia(     { 2, 1, 0, 2, 1, 1, 0, 0, 0 } );
    HArray<IndexType> ja(     { 2, 1, 0, 1, 1, 1, 0, 2, 1 } );
    HArray<ValueType> values( { 1, 2, 3, 4, 5, 6, 7, 8, 9 } );

    // This is the expected COO data

    HArray<IndexType> ia1(     { 0, 0, 0, 0, 1, 1, 1, 2, 2 } );
    HArray<IndexType> ja1(     { 0, 0, 1, 2, 1, 1, 1, 1, 2 } );
    HArray<ValueType> values1( { 3, 7, 9, 8, 2, 5, 6, 4, 1 } );

    COOUtils::sort( ia, ja, values, testContext );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( ia1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( ja1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( values1 ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( isSortedTest )
{
    ContextPtr testContext = ContextFix::testContext;

    // Here is some sorted COO data

    HArray<IndexType> ia(     { 0, 0, 0, 0, 1, 1, 1, 2, 2 } );
    HArray<IndexType> ja(     { 0, 1, 3, 5, 1, 3, 5, 1, 2 } );

    BOOST_CHECK( COOUtils::isSorted( ia, ja, testContext ) );

    ia     = { 2, 1, 0, 2, 1, 1, 0, 0, 0 };
    ja     = { 2, 1, 0, 1, 1, 1, 0, 2, 1 };

    BOOST_CHECK( !COOUtils::isSorted( ia, ja, testContext ) );

    ia.clear();
    ja.clear();

    BOOST_CHECK( COOUtils::isSorted( ia, ja, testContext ) );

    ia = { 1 };
    ja = { 1 };

    BOOST_CHECK( COOUtils::isSorted( ia, ja, testContext ) );

    // double elements, is also not sorted

    ia = { 0, 0, 1, 1, 3 };
    ja = { 0, 3, 1, 1, 1 };

    BOOST_CHECK( !COOUtils::isSorted( ia, ja, testContext ) );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( uniqueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    // Here is the sorted COO data with double entries

    HArray<IndexType> ia(     { 0, 0, 0, 0, 1, 1, 1, 2, 2 } );
    HArray<IndexType> ja(     { 0, 0, 1, 2, 1, 1, 1, 1, 2 } );
    HArray<ValueType> values( { 3, 7, 9, 8, 2, 5, 6, 4, 1 } );

    // This is the unique COO data

    HArray<IndexType> ia1(     { 0, 0, 0, 1, 2, 2 } );
    HArray<IndexType> ja1(     { 0, 1, 2, 1, 1, 2 } );
    HArray<ValueType> values1( { 7, 9, 8, 6, 4, 1 } );

    COOUtils::unique( ia, ja, values, common::BinaryOp::COPY, testContext );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( ia1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( ja1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( values1 ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( uniqueOpTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    // Here is the sorted COO data with double entries

    HArray<IndexType> ia(     { 0, 0, 0, 0, 1, 1, 1, 2, 2 } );
    HArray<IndexType> ja(     { 0, 0, 1, 2, 1, 1, 1, 1, 2 } );
    HArray<ValueType> values( { 3, 7, 9, 8, 2, 5, 6, 4, 1 } );

    // This is the unique COO data

    HArray<IndexType> ia1(     {  0, 0, 0,  1, 2, 2 } );
    HArray<IndexType> ja1(     {  0, 1, 2,  1, 1, 2 } );
    HArray<ValueType> values1( { 10, 9, 8, 13, 4, 1 } );

    COOUtils::unique( ia, ja, values, common::BinaryOp::ADD, testContext );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( ia1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( ja1 ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( values1 ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( diagonalTest )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    //   Input storage:
    //   -------------
    //    1.0   -   2.0  1.1  - 
    //    0.5  0.0   -    -   - 
    //     -    -   3.0   -   - 
    //     -   0.0  4.0  0.0 1.0

    HArray<IndexType> ia(     {   0,   0,   0,   1,   1,   2,   3,  3,  3,   3 }, testContext );
    HArray<IndexType> ja(     {   0,   2,   3,   0,   1,   2,   1,  2,  3,   4   },  testContext );
    HArray<ValueType> values( { 1.0, 2.0, 1.1, 0.5, 0.0, 3.0, 0.0, 4.0, 0.0, 1.0   },  testContext );

    HArray<ValueType> expDiag( { 1.0, 0.0, 3.0, 0.0 } );

    const IndexType m = 4;
    const IndexType n = 5;

    HArray<ValueType> diag;

    COOUtils::getDiagonal( diag, m, n, ia, ja, values, testContext );

    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( expDiag ), per_element() );

    HArray<ValueType> newDiag( { 1.2, 2.0, 3.3, 0.5 } );

    COOUtils::setDiagonalV( values, newDiag, m, n, ia, ja, testContext );

    COOUtils::getDiagonal( diag, m, n, ia, ja, values, testContext );

    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( newDiag ), per_element() );

    ValueType diagVal = 1;

    COOUtils::setDiagonal( values, diagVal, m, n, ia, ja, testContext );

    COOUtils::getDiagonal( diag, m, n, ia, ja, values, testContext );

    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( HArray<ValueType>( m, diagVal) ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( binaryOpTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    //       array 1                array 2
    //    1.0   -   2.0  -         1.0  0.5   -   -
    //    0.5  0.3   -   0.5        -    -    -   0.5
    //     -    -   3.0  -        2.0   -   1.0   0.5

    const IndexType m = 3;
    const IndexType n = 4;

    // csr arrays for matrix a, MUST BE sorted

    HArray<IndexType> aIA(     { 0,   0,   1,   1,   1,   2  }, testContext );
    HArray<IndexType> aJA(     { 0,   2,   0,   1,   3,   2  }, testContext );
    HArray<ValueType> aValues( { 1.0, 2.0, 0.5, 0.3, 0.5, 3.0 }, testContext );

    // csr arrays for matrix b, MUST BE sorted

    HArray<IndexType> bIA(     { 0,   0,   1,   2,   2,  2    }, testContext );
    HArray<IndexType> bJA(     { 0,   1,   3,   0,   2,  3    }, testContext );
    HArray<ValueType> bValues( { 1.0, 0.5, 0.5, 2.0, 1.0, 0.5 }, testContext );

    //       array3 = array 1 + array 2
    //
    //     2.0  0.5  2.0  -
    //     0.5  0.3   -   1.0
    //     2.0   -   4.0  0.5

    HArray<IndexType> expCIA(     { 0,   0,   0,   1,   1,   1,   2,   2,   2   } );
    HArray<IndexType> expCJA(     { 0,   1,   2,   0,   1,   3,   0,   2,   3   } );
    HArray<ValueType> expCValues( { 2.0, 0.5, 2.0, 0.5, 0.3, 1.0, 2.0, 4.0, 0.5 } );

    HArray<IndexType> cIA;
    HArray<IndexType> cJA;
    HArray<ValueType> cValues;

    COOUtils::binaryOp( cIA, cJA, cValues, aIA, aJA, aValues, bIA, bJA, bValues, m, n, common::BinaryOp::ADD, testContext );

    // Note: entries in storage C are SORTED

    BOOST_CHECK( COOUtils::isSorted( cIA, cJA, testContext ) );

    BOOST_TEST( hostReadAccess( cIA ) == hostReadAccess( expCIA ), per_element() );
    BOOST_TEST( hostReadAccess( cJA ) == hostReadAccess( expCJA ), per_element() );
    BOOST_TEST( hostReadAccess( cValues ) == hostReadAccess( expCValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()

