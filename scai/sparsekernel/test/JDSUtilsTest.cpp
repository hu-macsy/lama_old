/**
 * @file JDSUtilsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Contains tests for the class CUDAJDSUtils and OpenMPJDSUtils
 * @author Thomas Brandes
 * @date 16.10.2012
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/kregistry.hpp>
#include <scai/utilskernel.hpp>
#include <scai/sparsekernel/JDSUtils.hpp>
#include <scai/hmemo.hpp>
#include <scai/sparsekernel/test/TestMacros.hpp>

#include <scai/sparsekernel/test/TestData1.hpp>
#include <scai/sparsekernel/test/TestData2.hpp>
#include <scai/tasking/SyncToken.hpp>

#include <scai/hmemo/test/ContextFix.hpp>

#include <memory>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;
using namespace kregistry;
using common::TypeTraits;
using common::Exception;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( JDSUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.JDSUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getRowTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    const IndexType i = 2;

    const IndexType numColumns = 16;

    HArray<ValueType> values( { 1, 7, 12, 2, 8, 13, 3, 9, 14, 4, 10, 15, 5, 11, 6 }, testContext );
    HArray<IndexType> ja(  { 0, 1, 5, 2, 3, 7, 4, 5, 12, 6, 7, 15, 8, 9, 10 }, testContext );
    HArray<IndexType> dlg( { 3, 3, 3, 3, 2, 1 }, testContext );
    HArray<IndexType> ilg( { 6, 5, 4 }, testContext );
    HArray<IndexType> perm( { 1, 2, 0 }, testContext );

    HArray<ValueType> row;

    JDSUtils::getRow( row, numColumns, i, ilg, dlg, perm, ja, values, testContext );

    std::vector<ValueType> expectedValues( { 0, 7, 0, 8, 0, 9, 0, 10, 0, 11, 0, 0, 0, 0, 0, 0 } );

    BOOST_TEST( hostReadAccess( row ) == expectedValues, per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( getRowPositionsTest )
{
    ContextPtr testContext = ContextFix::testContext;

    typedef DefaultReal ValueType;

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    HArray<IndexType> pos;

    IndexType n = 0;    // count number of entries

    for ( IndexType i = 0; i < numRows; ++i )
    {
        JDSUtils::getRowPositions( pos, jdsILG, jdsDLG, jdsPerm, i, testContext );

        n += pos.size();
    }

    BOOST_CHECK_EQUAL( jdsJA.size(), n );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setRowTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    HArray<ValueType> row;

    auto op = common::BinaryOp::SUB;

    for ( IndexType i = 0; i < numRows; ++i )
    {
        JDSUtils::getRow( row, numColumns, i, jdsILG, jdsDLG, jdsPerm, jdsJA, jdsValues, testContext );
        JDSUtils::setRow( jdsValues, i, row, jdsILG, jdsDLG, jdsPerm, jdsJA, op, testContext );
    }

    // Now all values should be zero

    BOOST_TEST( hostReadAccess( jdsValues ) == std::vector<ValueType>( jdsValues.size(), 0 ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getColumnPositionsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    IndexType nTotal = 0;

    for ( IndexType j = 0; j < numColumns; ++j )
    {
        HArray<IndexType> row;
        HArray<IndexType> pos;

        JDSUtils::getColumnPositions( row, pos, jdsILG, jdsDLG, jdsPerm, jdsJA, j, testContext );
 
        BOOST_CHECK_EQUAL( row.size(), pos.size() );

        nTotal += row.size();
    }

    BOOST_CHECK_EQUAL( nTotal, jdsJA.size() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    const IndexType numColumns = 10;
    const IndexType numRows    = 5;

    // the dense storage data, stored row-wise

    ValueType expectedValues[5][10] =
    {
        { 1, 3, 0, 4, 0, 0, 7, 0, 0, 2 },
        { 0, 0, 3, 0, 2, 0, 0, 0, 0, 0 },
        { 5, 0, 0, 2, 0, 0, 0, 9, 0, 8 },
        { 0, 0, 4, 0, 0, 2, 9, 0, 0, 7 },
        { 1, 8, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    // here are the JDS arrays

    HArray<IndexType> dlg(  { 5, 5, 3, 3, 1 }, testContext );
    HArray<IndexType> ilg(  { 5, 4, 4, 2, 2 },  testContext );
    HArray<IndexType> perm( { 0, 2, 3, 1, 4 }, testContext );

    HArray<ValueType> values( { 1, 5, 4, 3, 1, 3, 2, 2, 2, 8, 4, 9, 9, 7, 8, 7, 2 }, testContext );
    HArray<IndexType> ja(     { 0, 0, 2, 2, 0, 1, 3, 5, 4, 1, 3, 7, 6, 6, 9, 9, 9 }, testContext );

    auto rValues = hostReadAccess( values );

    IndexType nnz = 0;

    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType j = 0; j < numColumns; j++ )
        {
            IndexType pos = JDSUtils::getValuePos( i, j, ilg, dlg, perm, ja, testContext );

            if ( pos == invalidIndex )
            {
                BOOST_CHECK_EQUAL( expectedValues[i][j], 0 );
            }
            else
            {
                BOOST_CHECK_EQUAL( expectedValues[i][j], rValues[pos] );
                nnz++;
            }
        }
    }

    BOOST_CHECK_EQUAL( values.size(), nnz );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    /* Example jds storage:
     *
     *      1  2   3   4   5  6    |  row 1  * 1
     *      7  8   9  10  11       |  row 2  * 2
     *     12 13  14  15           |  row 3  * 3
     */

    HArray<ValueType> values( { 1, 7, 12, 2, 8, 13, 3, 9, 14, 4, 10, 15, 5, 11, 6 }, testContext );
    HArray<IndexType> dlg( { 3, 3, 3, 3, 2, 1 }, testContext );
    HArray<IndexType> ilg( { 6, 5, 4 }, testContext );
    HArray<IndexType> perm( { 1, 2, 0 },  testContext );
    HArray<ValueType> diagonal( { 3, 1, 2 }, testContext );

    common::BinaryOp op = common::BinaryOp::MULT;

    JDSUtils::setRows( values, ilg, dlg, perm, diagonal, op, testContext );

    HArray<ValueType> expectedValues( { 1, 14, 36, 2, 16, 39, 3, 18, 42, 4, 20, 45, 5, 22, 6 } );

    BOOST_TEST( hostReadAccess( expectedValues ) == hostReadAccess( values ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( diagonalTest )
{
    ContextPtr testContext = ContextFix::testContext;

    typedef DefaultReal ValueType;

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data2::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    HArray<IndexType> diagonalPositions;

    JDSUtils::getDiagonalPositions( diagonalPositions, numRows, numColumns, jdsILG, jdsDLG, jdsPerm, jdsJA, testContext );

    HArray<ValueType> expDiagonalPositions( { 0, 1, 6, 3 } );

    BOOST_TEST( hostReadAccess( expDiagonalPositions ) == hostReadAccess( diagonalPositions ), per_element() );

    HArray<ValueType> diagonal;

    JDSUtils::getDiagonal( diagonal, numRows, numColumns, jdsILG, jdsDLG, jdsPerm, jdsJA, jdsValues, testContext );

    HArray<ValueType> expDiagonal( { 6, 8, 9, 3 } );

    BOOST_TEST( hostReadAccess( expDiagonal ) == hostReadAccess( diagonal ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( checkDiagonalPropertyTest )
{
    ContextPtr testContext = ContextFix::testContext;

    // check JDS storage without diagonal property

    {
        /*
               0  2  4   6  8  10    row 1
               1  3  5   7  9        row 2 
               5  7  12 15           row 0    */

        HArray<IndexType> ja( { 0, 1, 5, 2, 3, 7, 4, 5, 12, 6, 7, 15, 8, 9, 10 } );
        HArray<IndexType> dlg( { 3, 3, 3, 3, 2, 1 } );
        HArray<IndexType> ilg( { 6, 5, 4 } );
        HArray<IndexType> perm( { 1, 2, 0 } );

        const IndexType numRows = 3;
        const IndexType numColumns = 16;

        HArray<IndexType> diagonalPositions;

        IndexType nFound = JDSUtils::getDiagonalPositions( diagonalPositions, numRows, numColumns, ilg, dlg, perm, ja, testContext );

        BOOST_CHECK_EQUAL( nFound, 0 );
    }

    // check JDS storage with diagonal property

    {
        HArray<IndexType> ja( { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } );
        HArray<IndexType> dlg( { 3, 3, 3 } );
        HArray<IndexType> ilg( { 3, 3, 3 } );
        HArray<IndexType> perm( { 0, 1, 2 } );

        const IndexType numRows = 3;
        const IndexType numColumns = 3;

        HArray<IndexType> diagonalPositions;

        IndexType nFound = JDSUtils::getDiagonalPositions( diagonalPositions, numRows, numColumns, ilg, dlg, perm, ja, testContext );

        BOOST_CHECK_EQUAL( nFound, 3 );

        BOOST_TEST( hostReadAccess( diagonalPositions ) == std::vector<IndexType>( { 0, 1, 2 } ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( ilg2dlgTest )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> ilg( { 7, 7, 5, 4, 4, 1 }, testContext );

    HArray<IndexType> dlg;   // result array

    IndexType numValues = JDSUtils::ilg2dlg( dlg, ilg, testContext );

    BOOST_CHECK_EQUAL( numValues, HArrayUtils::reduce( dlg, common::BinaryOp::ADD ) );
    BOOST_CHECK_EQUAL( numValues, HArrayUtils::reduce( ilg, common::BinaryOp::ADD ) );

    // Expected results: 6 rows have at least 1 entry, 5 rows have at least 2, 3, 4 entries, 
    //                   3 rows have hat least 5 entries, 2 rows have at least 6, 7 entries

    std::vector<IndexType> expectedDlg( { 6, 5, 5, 5, 3, 2, 2 } );

    BOOST_TEST( hostReadAccess( dlg ) == expectedDlg, per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( convertCSRTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    /*
     * Testmatrix:
     * 0 0 5 3 0 0 4 0
     * 3 0 4 0 3 5 0 0
     * 0 2 0 8 7 9 0 5
     * 2 0 3 0 0 0 0 0
     * 5 0 0 7 0 0 0 9
     */

    HArray<IndexType> expJDSDlg( { 5, 5, 4, 2, 1 } );
    HArray<IndexType> expJDSIlg( { 5, 4, 3, 3, 2 } );
    HArray<IndexType> expJDSPerm(  { 2, 1, 0, 4, 3 } );

    HArray<IndexType> expJDSJa(     { 1, 0, 2, 0, 0, 3, 2, 3, 3, 2, 4, 4, 6, 7, 5, 5, 7 } );
    HArray<ValueType> expJDSValues( { 2, 3, 5, 5, 2, 8, 4, 3, 7, 3, 7, 3, 4, 9, 9, 5, 5 } );

    const IndexType numRows = 5;
    const IndexType numColumns = 8;

    HArray<IndexType> CSRIa( { 0, 3, 7, 12, 14, 17 }  );
    HArray<IndexType> CSRJa(  { 2, 3, 6, 0, 2, 4, 5, 1, 3, 4, 5, 7, 0, 2, 0, 3, 7 } );
    HArray<ValueType> CSRValues( { 5, 3, 4, 3, 4, 3, 5, 2, 8, 7, 9, 5, 2, 3, 5, 7, 9 } );

    // output arrays for JDS Storage

    HArray<IndexType> JDSIlg;
    HArray<IndexType> JDSDlg;
    HArray<IndexType> JDSPerm;

    HArray<IndexType> JDSJa;
    HArray<ValueType> JDSValues;

    JDSUtils::convertCSR2JDS( JDSIlg, JDSDlg, JDSPerm, JDSJa, JDSValues,
                              numRows, numColumns, CSRIa, CSRJa, CSRValues, testContext );

    BOOST_TEST( hostReadAccess( JDSIlg ) == hostReadAccess( expJDSIlg ), per_element() );
    BOOST_TEST( hostReadAccess( JDSDlg ) == hostReadAccess( expJDSDlg ), per_element() );
    BOOST_TEST( hostReadAccess( JDSPerm ) == hostReadAccess( expJDSPerm ), per_element() );

    BOOST_TEST( hostReadAccess( JDSJa ) == hostReadAccess( expJDSJa ), per_element() );
    BOOST_TEST( hostReadAccess( JDSValues ) == hostReadAccess( expJDSValues ), per_element() );

    // output arrays for new CSR Storage

    HArray<IndexType> newCSRIa;
    HArray<IndexType> newCSRJa;
    HArray<ValueType> newCSRValues;

    JDSUtils::convertJDS2CSR( newCSRIa, newCSRJa, newCSRValues, 
                              numRows, numColumns, 
                              JDSIlg, JDSDlg, JDSPerm, JDSJa, JDSValues,
                              testContext );

    BOOST_TEST( hostReadAccess( CSRIa ) == hostReadAccess( newCSRIa ), per_element() );
    BOOST_TEST( hostReadAccess( CSRJa ) == hostReadAccess( newCSRJa ), per_element() );
    BOOST_TEST( hostReadAccess( CSRValues ) == hostReadAccess( newCSRValues ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    SCAI_ASSERT_EQ_ERROR( jdsPerm.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsILG.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsDLG.size(), numDiagonals, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsJA.size(), jdsValues.size(), "size mismatch" )

    HArray<ValueType> x( { 3, -3, 2, -2 }, testContext );
    HArray<ValueType> y( { 1, -1, 2, -2, 1, 1, -1 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numColumns, x.size(), "illegal array x for gemv" );
    SCAI_ASSERT_EQ_ERROR( numRows, y.size(), "illegal array y for gemv" );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    HArray<ValueType> res( testContext );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * JDS * x + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", JDS: ilg = " << jdsILG << ", ja = " << jdsJA << ", values = " << jdsValues )

        bool async = false;

        auto op = common::MatrixOp::NORMAL;

        JDSUtils::gemv( res, alpha, x, beta, y, numRows, numColumns, 
                        jdsILG, jdsDLG, jdsPerm, jdsJA, jdsValues, op, async, testContext );

        // compare against mv product of dense matrix

        auto expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEMVTest, ValueType, scai_numeric_test_types )
{
    // JDS has no own sparse routine, also done with normalGEMV, beta = 1, y = res

    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    SCAI_ASSERT_EQ_ERROR( jdsPerm.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsILG.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsDLG.size(), numDiagonals, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsJA.size(), jdsValues.size(), "size mismatch" )

    const HArray<ValueType> x( { 3, -3, 2, -2 }, testContext );
    const HArray<ValueType> y( { 1, -1, 2, -2, 1, 1, -1 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numColumns, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, y.size(), "size mismatch" );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha; ++icase )
    {
        HArray<ValueType> res( y );

        ValueType alpha = alpha_values[icase];
        ValueType beta  = 1;

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * JDS * x + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", JDS: ilg = " << jdsILG << ", ja = " << jdsJA << ", values = " << jdsValues )

        bool async = false;

        auto op = common::MatrixOp::NORMAL;

        JDSUtils::gemv( res, alpha, x, beta, y, numRows, numColumns, 
                        jdsILG, jdsDLG, jdsPerm, jdsJA, jdsValues, op, async, testContext );

        // compare against mv product of dense matrix

        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_CHECK( HArrayUtils::maxDiffNorm( expectedRes, res ) < 0.001 );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTransposeTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    SCAI_ASSERT_EQ_ERROR( jdsPerm.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsILG.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsDLG.size(), numDiagonals, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( jdsJA.size(), jdsValues.size(), "size mismatch" )

    HArray<ValueType> x( { 3, -2, -2, 3, 1, 0, 1 }, testContext );
    HArray<ValueType> y( { 1, -1, 2, -2 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numRows, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, y.size(), "size mismatch" );

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    HArray<ValueType> res( testContext );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x * JDS + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", JDS: ilg = " << jdsILG << ", ja = " << jdsJA << ", values = " << jdsValues )

        bool async = false;

        auto op = common::MatrixOp::TRANSPOSE;

        JDSUtils::gemv( res, alpha, x, beta, y, numRows, numColumns, 
                        jdsILG, jdsDLG, jdsPerm, jdsJA, jdsValues, op, async, testContext );

        // compare against gevm product of dense matrix

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_CHECK( HArrayUtils::maxDiffNorm( expectedRes, res ) < 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext )

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data2::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    HArray<ValueType> rhs( { 1, -1, 2, -2 }, testContext );
    HArray<ValueType> oldSolution( { 3, -2, -2, 3 }, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {   
        ValueType omega  = omega_values[icase];
        
        HArray<ValueType> newSolution( testContext );
        
        bool async = false;

        JDSUtils::jacobi( newSolution, omega, oldSolution, rhs, 
                          jdsILG, jdsDLG, jdsPerm, jdsJA, jdsValues, async, testContext );
        
        HArray<ValueType> expSolution;
        
        data2::getJacobiResult( expSolution, oldSolution, omega, rhs );
        
        auto eps = common::TypeTraits<ValueType>::small();
        
        auto maxDiff = HArrayUtils::maxDiffNorm( expSolution, newSolution );
        
        BOOST_CHECK( maxDiff < eps );
        
        if ( maxDiff >= eps )
        {   
            // compare the individual values to see what went wrong

            BOOST_TEST( hostReadAccess( expSolution ) == hostReadAccess( newSolution ), per_element() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiHaloTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

    // Note: for jds we do not need the number of non empty rows, as this comes for free

    const ValueType old_values[]   = { 3, -2, -2, 3, 1, 0, 2 };
    const ValueType diag_values[]  = { 9,  8,  7, 6, 7, 8, 9 };

    const IndexType n_old_values = sizeof( old_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( n_old_values, numRows, "size mismatch" );

    HArray<ValueType> oldSolution( numRows, old_values, testContext );
    HArray<ValueType> diag( numRows, diag_values, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> solution( numRows, ValueType( 0 ), testContext );
 
        JDSUtils::jacobiHalo( solution, omega, oldSolution, diag, 
                              jdsILG, jdsDLG, jdsPerm, jdsJA, jdsValues, testContext );

        auto expectedSol = utilskernel::fillHArray<ValueType>( numRows, 0, testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        // be careful: expectedSol == oldSolution does not hold on GPUs

        BOOST_CHECK( HArrayUtils::maxDiffNorm( expectedSol, solution ) < 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
