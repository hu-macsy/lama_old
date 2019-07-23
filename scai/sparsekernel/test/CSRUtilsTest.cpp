/**
 * @file CSRUtilsTest.cpp
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
 * @brief Contains tests for the CSRUtils interface to be tested on different devices
 * @author Thomas Brandes
 * @date 05.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/kregistry.hpp>
#include <scai/utilskernel.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>
#include <scai/common/Settings.hpp>
#include <scai/sparsekernel/test/TestMacros.hpp>

#include <scai/hmemo/test/ContextFix.hpp>

#include <scai/sparsekernel/test/TestData1.hpp>
#include <scai/sparsekernel/test/TestData2.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;
using common::TypeTraits;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CSRUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( validOffsetsTest )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( { 0, 2, 5, 7 }, testContext );

    IndexType numRows = csrIA.size() - 1;
    IndexType total   = csrIA[ numRows ];

    BOOST_CHECK( CSRUtils::validOffsets( csrIA, total, testContext ) );

    csrIA.resize( numRows -  1 );
    BOOST_CHECK( !CSRUtils::validOffsets( csrIA, total, testContext ) );

    csrIA =  { 0, 5, 2, 7 };
    BOOST_CHECK( !CSRUtils::validOffsets( csrIA, total, testContext ) );

    csrIA =  {};
    BOOST_CHECK( !CSRUtils::validOffsets( csrIA, 0, testContext ) );

    csrIA =  { 0 };
    BOOST_CHECK( CSRUtils::validOffsets( csrIA, 0, testContext ) );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( offsets2sizesTest )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( { 0, 5, 11, 16, 19, 19, 21, 23 }, testContext );
    HArray<IndexType> expSizes( { 5, 6,  5,  3,  0,  2,  2 } );

    HArray<IndexType> csrSizes;

    CSRUtils::offsets2sizes( csrSizes, csrIA, testContext );
  
    BOOST_TEST( hostReadAccess( csrSizes ) == hostReadAccess( expSizes ), per_element() );

    // might also be aliased

    CSRUtils::offsets2sizes( csrIA, csrIA, testContext );
 
    BOOST_TEST( hostReadAccess( csrIA ) == hostReadAccess( expSizes ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( sizes2offsetTest )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrSizes( { 5, 6,  5,  3,  0,  2,  2 } );
    HArray<IndexType> expIA( { 0, 5, 11, 16, 19, 19, 21, 23 }, testContext );

    HArray<IndexType> csrIA;

    CSRUtils::sizes2offsets( csrIA, csrSizes, testContext );
 
    BOOST_TEST( hostReadAccess( csrIA ) == hostReadAccess( expIA ), per_element() );

    // might also be aliased

    CSRUtils::sizes2offsets( csrSizes, csrSizes, testContext );
 
    BOOST_TEST( hostReadAccess( csrSizes ) == hostReadAccess( expIA ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( gatherSizesTest )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( { 0, 5, 11, 16, 19, 19, 21, 23 }, testContext );

    HArray<IndexType> rowIndexes( { 6, 5, 4, 3, 2, 1, 0 }, testContext );
    HArray<IndexType> expSizes( { 2, 2, 0, 3, 5, 6, 5 } );

    HArray<IndexType> sizes;

    CSRUtils::gatherSizes( sizes, csrIA, rowIndexes, testContext );
    BOOST_TEST( hostReadAccess( sizes ) == hostReadAccess( expSizes ), per_element() );

    rowIndexes = { 0, 2, 4, 6 };
    expSizes = { 5, 5, 0, 2 };

    CSRUtils::gatherSizes( sizes, csrIA, rowIndexes, testContext );
    BOOST_TEST( hostReadAccess( sizes ) == hostReadAccess( expSizes ), per_element() );

    rowIndexes = {};
    expSizes = {};

    CSRUtils::gatherSizes( sizes, csrIA, rowIndexes, testContext );
    BOOST_TEST( hostReadAccess( sizes ) == hostReadAccess( expSizes ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( nonEmptyRowsTest )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( { 0, 2, 2, 5, 5, 7, 7, 9 }, testContext );
    HArray<IndexType> expRowIndexes( { 0,    2,    4,    6  } );

    HArray<IndexType> rowIndexes;  // will be set if there are less than threshold non-zero rows

    CSRUtils::nonEmptyRows( rowIndexes, csrIA, 1.0f, testContext );

    BOOST_TEST( hostReadAccess( rowIndexes ) == hostReadAccess( expRowIndexes ), per_element() );

    // set a threshold, as there are more than 20% non-empty row, no indexes are built.

    CSRUtils::nonEmptyRows( rowIndexes, csrIA, 0.2f, testContext );

    BOOST_CHECK_EQUAL( rowIndexes.size(), IndexType( 0 ) );

    csrIA = { 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2 };
    expRowIndexes = { 0, 6 };

    // here we have less than 20% non-empty rows

    CSRUtils::nonEmptyRows( rowIndexes, csrIA, 0.2f, testContext );

    BOOST_TEST( hostReadAccess( rowIndexes ) == hostReadAccess( expRowIndexes ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( hasDiagonalPropertyTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    bool hasSortedRows = false;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    bool okay = CSRUtils::hasDiagonalProperty( numRows, numColumns, csrIA, csrJA, hasSortedRows, testContext );

    BOOST_CHECK( !okay );

    // data set 2 has a square matrix with diagonal entries first

    data2::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    BOOST_REQUIRE_EQUAL( numRows, numColumns );

    okay = CSRUtils::hasDiagonalProperty( numRows, numColumns, csrIA, csrJA, hasSortedRows, testContext );

    BOOST_CHECK( okay );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( shiftDiagonalFirstTest )
{
    ContextPtr testContext = ContextFix::testContext;

    typedef SCAI_TEST_TYPE ValueType;

    // CSR data of square matrix 6 x 6, diagonal elements in row 0, 2, 4, 5
    //
    //     0   1   2  -  -  -
    //     -   -   -  -  -  -
    //     0   1   2  -  -  -
    //     0   1   2  -  -  -
    //     -   -   -  -  4  -
    //     -   -   -  -  4  5

    HArray<IndexType> csrIA(     { 0,       3, 3,       6,       9, 10,  12 }, testContext );
    HArray<IndexType> csrJA(     { 0, 1, 2,    0, 1, 2, 0, 1, 2, 4, 4, 5 }, testContext );
    HArray<ValueType> csrValues( { 0, 1, 2,    0, 1, 2, 0, 1, 2, 4, 4, 5 }, testContext );

    const IndexType numRows    = 6;
    const IndexType numColumns = 6;

    IndexType numDiags = CSRUtils::shiftDiagonalFirst( csrJA, csrValues, numRows, numColumns, csrIA, testContext );

    BOOST_CHECK_EQUAL( numDiags, 4 );

    HArray<IndexType> expJA(     { 0, 1, 2, 2, 0, 1, 0, 1, 2, 4, 5, 4 }, testContext );
    HArray<ValueType> expValues( { 0, 1, 2, 2, 0, 1, 0, 1, 2, 4, 5, 4 }, testContext );

    BOOST_TEST( hostReadAccess( csrJA ) == hostReadAccess( expJA ), per_element() );
    BOOST_TEST( hostReadAccess( csrValues ) == hostReadAccess( expValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( shiftDiagonalFirst1Test )
{
    ContextPtr testContext = ContextFix::testContext;

    typedef SCAI_TEST_TYPE ValueType;

    // CSR data of square matrix 6 x 6, diagonal elements in row 0, 2, 4, 5
    //
    //     0   1   2  -  -  -
    //     -   -   -  -  -  -
    //     0   1   2  -  -  -
    //     0   1   2  -  -  -
    //     -   -   -  -  4  -
    //     -   -   -  -  4  5

    HArray<IndexType> csrIA(     { 0,       3, 3,       6,       9, 10,  12 }, testContext );
    HArray<IndexType> csrJA(     { 0, 1, 2,    0, 1, 2, 0, 1, 2, 4, 4, 5 }, testContext );
    HArray<ValueType> csrValues( { 0, 1, 2,    0, 1, 2, 0, 1, 2, 4, 4, 5 }, testContext );

    const IndexType numRows    = 6;
    const IndexType numColumns = 6;

    HArray<IndexType> diagonals( { 2, 3, 1, 0, 5, 4 } );

    IndexType numDiags = CSRUtils::shiftDiagonalFirst( csrJA, csrValues, numRows, numColumns, csrIA, diagonals, testContext );

    BOOST_CHECK_EQUAL( numDiags, 4 );   // no shift in rows 1, 4

    HArray<IndexType> expJA(     { 2, 0, 1, 1, 0, 2, 0, 1, 2, 4, 4, 5 }, testContext );
    HArray<ValueType> expValues( { 2, 0, 1, 1, 0, 2, 0, 1, 2, 4, 4, 5 }, testContext );

    BOOST_TEST( hostReadAccess( csrJA ) == hostReadAccess( expJA ), per_element() );
    BOOST_TEST( hostReadAccess( csrValues ) == hostReadAccess( expValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxDiffValTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "absMaxDiffVal< " << TypeTraits<ValueType>::id() << "> test for " << *testContext )

    // input arrays
    //    Array1             Array2
    //
    //    1 2 3 0 0          1 2 0 0 0
    //    0 0 1 1 2          1 0 2 2 1

    const IndexType numRows = 2;
    const IndexType numColumns = 5;

    HArray<IndexType> csrIA1( { 0, 3, 6 }, testContext );
    HArray<IndexType> csrJA1( { 0, 1, 2, 2, 3, 4 }, testContext );
    HArray<ValueType> csrValues1( { 1, 2, 3, 1, 1, 2 }, testContext );
    HArray<IndexType> csrIA2( { 0, 2, 6 }, testContext );
    HArray<IndexType> csrJA2(  { 0, 1, 0, 2, 3, 4 }, testContext );
    HArray<ValueType> csrValues2(  { 1, 2, 1, 2, 2, 1 }, testContext );

    auto maxVal = CSRUtils::absMaxDiffVal( csrIA1, csrJA1, csrValues1, 
                                           csrIA2, csrJA2, csrValues2, 
                                           numRows, numColumns, false, testContext );

    BOOST_CHECK_EQUAL( 3, maxVal );

    // rows are sorted, so we can also apply sortFlag = true

    maxVal = CSRUtils::absMaxDiffVal( csrIA1, csrJA1, csrValues1, 
                                      csrIA2, csrJA2, csrValues2, 
                                      numRows, numColumns, true, testContext );

    BOOST_CHECK_EQUAL( 3, maxVal );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "setRows test for " << *testContext )

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    HArray<ValueType> savedValues( csrValues );  // keep a copy for comparison later

    HArray<ValueType> rows(  { 2, 3, 4, 5, 1, 3, 2 }, testContext );
    BOOST_REQUIRE_EQUAL( numRows, rows.size() );

    CSRUtils::setRows( csrValues, numRows, numColumns, csrIA, csrJA, rows, common::BinaryOp::MULT, testContext );

    // prove by hand

    {
        auto rIA = hostReadAccess( csrIA );

        auto rRows = hostReadAccess( rows );
        auto rSavedValues = hostReadAccess( savedValues );
        auto rValues = hostReadAccess( csrValues );

        for ( IndexType i = 0; i < numRows; ++ i )
        {
            for ( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
            {
                BOOST_CHECK_EQUAL( rRows[i] * rSavedValues[jj], rValues[jj] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setColumnsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "setColumns test for " << *testContext )

    /*   Matrix:       6  0  0  4         0  3  -    6  4  -
                       7  0  0  0         0  -  -    7  -  -
                       0  0  9  4         2  3  -    9  4  -
                       2  5  0  3         0  1  3    2  5  3
                       2  0  0  1         0  3  -    2  1  -
                       0  0  0  0         -  -  -    -  -  -
                       0  1  0  2         1  3  -    1  2  -
     */

    HArray<IndexType> csrIA( { 0,    2, 3,    5,       8,    10, 10,   12 }, testContext );
    HArray<IndexType> csrJA( { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3,     1, 3 }, testContext );
    HArray<ValueType> csrValues( { 6, 4, 7, 9, 4, 2, 5, 3, 2, 1,     1, 2 }, testContext );

    const IndexType numRows = 7;
    const IndexType numColumns = 4;

    HArray<ValueType> columns(  { 2, 1, -1, 3 }, testContext );
    BOOST_REQUIRE_EQUAL( numColumns, columns.size() );

    CSRUtils::setColumns( csrValues, numRows, numColumns, csrIA, csrJA, columns, common::BinaryOp::MULT, testContext );

    HArray<ValueType> expectedValues( { 6 * 2, 4 * 3, 7 * 2, 9 * -1, 4 * 3, 2 * 2, 5 * 1, 3 * 3, 2 * 2, 1 * 3, 1 * 1, 2 * 3 }, testContext );

    // prove by hand

    SCAI_CHECK_EQUAL_ARRAY( expectedValues, csrValues );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeSquareTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> test for " << *testContext )

    //    1.0   -   2.0       1.0  0.5   -
    //    0.5  0.3   -         -   0.3   -
    //     -    -   3.0       2.0   -   3.0

    HArray<IndexType> csrIA( { 0, 2, 4, 5 }, testContext );
    HArray<IndexType> csrJA( { 0, 2, 0, 1, 2 }, testContext );
    HArray<ValueType> csrValues( { 1.0, 2.0, 0.5, 0.3, 3.0 }, testContext );

    HArray<IndexType> expIA( { 0, 2, 3, 5 }, testContext );
    HArray<IndexType> expJA(  { 0, 1, 1, 0, 2 }, testContext );
    HArray<ValueType> expValues( { 1.0, 0.5, 0.3, 2.0, 3.0 }, testContext );

    const HArray<ValueType> values2( { 1.0, 0.5, 0.3, 2.0, 3.0 } );

    const IndexType numRows = 3;
    const IndexType numColumns = 3;

    HArray<IndexType> cscIA;
    HArray<IndexType> cscJA;
    HArray<ValueType> cscValues;

    CSRUtils::convertCSR2CSC( cscIA, cscJA, cscValues, numRows, numColumns, csrIA, csrJA, csrValues, testContext );

    BOOST_TEST( hostReadAccess( cscIA ) == hostReadAccess( expIA ), per_element() );

    // for comparison we sort the values, take care of the switched dimensions

    CSRUtils::sortRows( cscJA, cscValues, numColumns, numRows, cscIA, testContext );

    BOOST_TEST( hostReadAccess( cscJA ) == hostReadAccess( expJA ), per_element() );
    BOOST_TEST( hostReadAccess( cscValues ) == hostReadAccess( expValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( sortRowsTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr testContext = ContextFix::testContext;

    //    1.0   -   2.0  1.1
    //    0.5  0.3   -    -
    //     -    -   3.0   -
    //     -    -   4.0   1.5

    HArray<IndexType> csrIA(     {   0,             3,         5,   6,        8 }, testContext );
    HArray<IndexType> csrJA(     {   3,   2,   0,   0,   1,    2,   3,   2 }, testContext );
    HArray<ValueType> csrValues( { 1.1, 2.0, 1.0, 0.5, 0.3 , 3.0, 1.5, 4.0 }, testContext );

    HArray<IndexType> sortedJA( { 0, 2, 3, 0, 1, 2, 2, 3 } );
    HArray<ValueType> sortedValues( { 1.0, 2.0, 1.1, 0.5, 0.3, 3.0, 4.0, 1.5 } );

    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    
    CSRUtils::sortRows( csrJA, csrValues, numRows, numColumns, csrIA, testContext );

    BOOST_TEST( hostReadAccess( csrJA ) == hostReadAccess( sortedJA ), per_element() );
    BOOST_TEST( hostReadAccess( csrValues ) == hostReadAccess( sortedValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( hasSortedRowsTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr testContext = ContextFix::testContext;

    const IndexType numRows = 4;
    const IndexType numColumns = 4;

    HArray<IndexType> csrIA( { 0,       3,    5, 6,    8 }, testContext );
    HArray<IndexType> csrJA( { 3, 2, 0, 1, 0, 2, 3, 2 }, testContext );

    HArray<ValueType> csrValues( csrJA.size(), 1 );

    // csrJA is not sorted as all

    BOOST_CHECK( !CSRUtils::hasSortedRows( numRows, numColumns, csrIA, csrJA, testContext ) );

    CSRUtils::sortRows( csrJA, csrValues, numRows, numColumns, csrIA, testContext );

    // csrJA is sorted in any case

    BOOST_CHECK( CSRUtils::hasSortedRows( numRows, numColumns, csrIA, csrJA, testContext ) );

    // CSR data with double entry is also not sorted

    csrIA = HArray<IndexType>( { 0,       3,    5,    6,  8 }, testContext );
    csrJA = HArray<IndexType>( { 0, 0, 3, 1, 2, 2, 3, 2 }, testContext );

    BOOST_CHECK( !CSRUtils::hasSortedRows( numRows, numColumns, csrIA, csrJA, testContext ) );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( compressTest )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    //    1.0   -   2.0  1.1
    //    0.5  0.0   -    -
    //     -    -   3.0   -
    //     -   0.0  4.0  0.0

    // uncompressed data

    HArray<IndexType> ia(     {   0,             3,        5,   6,          9 }, testContext );
    HArray<IndexType> ja(     {   3,   2,   0,   0,   1,   2,   3,  2,  1     },  testContext );
    HArray<ValueType> values( { 1.1, 2.0, 1.0, 0.5, 0.0, 3.0, 0.0, 4.0, 0.0   },  testContext );

    // expected compress data, row   0              1     2    3   

    HArray<IndexType> expIA(     {   0,              3,     4,   5,   6 } );
    HArray<IndexType> expJA(     {   3,    2,   0,   0,     2,   2 } );
    HArray<ValueType> expValues( {  1.1, 2.0, 1.0, 0.5,   3.0, 4.0 } );

    ValueType eps = common::TypeTraits<ValueType>::small();

    CSRUtils::compress( ia, ja, values, eps, testContext );

    BOOST_TEST( hostReadAccess( ia ) == hostReadAccess( expIA ), per_element() );
    BOOST_TEST( hostReadAccess( ja ) == hostReadAccess( expJA ), per_element() );
    BOOST_TEST( hostReadAccess( values ) == hostReadAccess( expValues), per_element() );
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

    HArray<IndexType> ia(     {   0,             3,        5,   6,           10  }, testContext );
    HArray<IndexType> ja(     {   3,   2,   0,   0,   1,   2,   3,  2,  1,   4    },  testContext );
    HArray<ValueType> values( { 1.0, 2.0, 1.5, 0.5, 0.0, 3.0, 0.0, 4.0, 0.0, 1.0   },  testContext );

    HArray<ValueType> expDiag( { 1.5, 0.0, 3.0, 0.0 } );

    const IndexType m = 4;
    const IndexType n = 5;

    bool isSorted = false;

    HArray<ValueType> diag;

    CSRUtils::getDiagonal( diag, m, n, ia, ja, values, isSorted, testContext );
 
    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( expDiag ), per_element() );

    HArray<ValueType> newDiag( { 1.2, 2.0, 3.3, 0.5 } );

    bool okay = CSRUtils::setDiagonalV( values, newDiag, m, n, ia, ja, isSorted, testContext );

    BOOST_CHECK( okay );

    CSRUtils::getDiagonal( diag, m, n, ia, ja, values, isSorted, testContext );

    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( newDiag ), per_element() );

    ValueType diagVal = 1;

    okay = CSRUtils::setDiagonal( values, diagVal, m, n, ia, ja, isSorted, testContext );

    CSRUtils::getDiagonal( diag, m, n, ia, ja, values, isSorted, testContext );

    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( HArray<ValueType>( m, diagVal) ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    HArray<ValueType> denseValues( testContext );

    data1::getDenseTestData( numRows, numColumns, denseValues );

    const ValueType ZERO = 0;

    // comparison is done via accesses on the host

    auto rValues = hostReadAccess( csrValues );
    auto rDense  = hostReadAccess( denseValues );

    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            IndexType pos = CSRUtils::getValuePos( i, j, csrIA, csrJA, testContext );

            IndexType k = i * numColumns + j;

            if ( pos == invalidIndex )
            {
                BOOST_CHECK_EQUAL( rDense[ k ], ZERO );
            }
            else
            {
                BOOST_CHECK_EQUAL( rDense[ k], rValues[pos] );
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

    HArray<IndexType> csrIA( { 0, 2, 4, 5 }, testContext );
    HArray<IndexType> csrJA( { 0, 2, 0, 1, 2 }, testContext );

    HArray<IndexType> row;   // result for rowIndexes
    HArray<IndexType> pos;   // result for positions

    IndexType columnIndex = 1;   // has 1 entry

    CSRUtils::getColumnPositions( row, pos, csrIA, csrJA, columnIndex, testContext );

    BOOST_TEST( hostReadAccess( row ) == hostReadAccess( HArray<IndexType>( { 1 } ) ), per_element() );
    BOOST_TEST( hostReadAccess( pos ) == hostReadAccess( HArray<IndexType>( { 3 } ) ), per_element() );

    columnIndex = 2;

    CSRUtils::getColumnPositions( row, pos, csrIA, csrJA, columnIndex, testContext );

    //  two entries for column 2, order might be arbitrary

    BOOST_REQUIRE_EQUAL( row.size(), IndexType( 2 ) );   //  two entries for column 2, order might be arbitrary
    BOOST_REQUIRE_EQUAL( pos.size(), IndexType( 2 ) );   //  two entries for column 2, order might be arbitrary

    // Verifiy: ellJA[ pos[ i ] ] == j, row[i] == pos[i] % numRows

    HArray<IndexType> ja;
    HArrayUtils::gather( ja, csrJA, pos, common::BinaryOp::COPY, testContext );
    BOOST_TEST( hostReadAccess( ja ) == std::vector<IndexType>( ja.size(), columnIndex ), per_element() );

    // sort sparse entries by row indexes, ascending = true

    HArrayUtils::sortSparseEntries( row, pos, true );

    BOOST_TEST( hostReadAccess( row ) == std::vector<IndexType>( { 0, 2 } ), per_element() );
    BOOST_TEST( hostReadAccess( pos ) == std::vector<IndexType>( { 1, 4 } ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( transposeZeroTest )
{
    typedef DefaultReal ValueType;

    // transpose a zero matrix of size 2 x 4

    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> zero test for " << *testContext )

    HArray<IndexType> csrIA( { 0, 0, 0 }, testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows     = 2;
    IndexType numColumns  = 4;

    // CSC <- transpose CSR

    IndexType dummyValue = 4019;

    HArray<IndexType> cscIA( numColumns + 5, dummyValue );
    HArray<IndexType> cscJA;
    HArray<ValueType> cscValues;

    CSRUtils::convertCSR2CSC( cscIA, cscJA, cscValues, numRows, numColumns, csrIA, csrJA, csrValues, testContext );

    BOOST_CHECK_EQUAL( cscIA.size(), numColumns + 1 );
    BOOST_CHECK_EQUAL( cscJA.size(), 0 );
    BOOST_CHECK_EQUAL( cscValues.size(), 0 );

    HArray<IndexType> expIA( numColumns + 1, IndexType( 0 ), testContext );

    SCAI_CHECK_EQUAL_ARRAY( cscIA, expIA )
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( transposeSingleTest )
{
    typedef DefaultReal ValueType;

    // transpose a matrix of size 3 x 5 with one entry at ( 1, 3 )

    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> zero test for " << *testContext )

    HArray<IndexType> csrIA( { 0, 0, 1, 1 }, testContext );
    HArray<IndexType> csrJA( 1, IndexType( 3 ), testContext );

    // be careful: csrValues( { 5 }, testContext ) might be seen as csrValues( 5, testContext ) by Intel compiler

    HArray<ValueType> csrValues( 1, ValueType( 5 ), testContext );

    IndexType numRows     = 3;
    IndexType numColumns  = 5;

    // CSC <- transpose CSR

    HArray<IndexType> cscIA;
    HArray<IndexType> cscJA;
    HArray<ValueType> cscValues;

    CSRUtils::convertCSR2CSC( cscIA, cscJA, cscValues, numRows, numColumns, csrIA, csrJA, csrValues, testContext );

    HArray<IndexType> expIA( { 0, 0, 0, 0, 1, 1 } );
    HArray<IndexType> expJA( 1, IndexType( 1 ) );

    SCAI_CHECK_EQUAL_ARRAY( cscValues, csrValues )
    SCAI_CHECK_EQUAL_ARRAY( cscIA, expIA )
    SCAI_CHECK_EQUAL_ARRAY( cscJA, expJA )
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeNonSquareTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> non-square test for " << *testContext )

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    SCAI_ASSERT_EQ_ERROR( csrIA.size(), numRows + 1, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrJA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrValues.size(), numValues, "size mismatch" )

    HArray<IndexType> cscIA;
    HArray<IndexType> cscJA;
    HArray<ValueType> cscValues;

    // CSC <- transpose CSR

    CSRUtils::convertCSR2CSC( cscIA, cscJA, cscValues, numRows, numColumns, csrIA, csrJA, csrValues, testContext );

    //  For comparison later we sort cscJA and cscValue

    SCAI_LOG_INFO( logger, "sortRows< " << TypeTraits<ValueType>::id() << "> for " << *testContext )

    CSRUtils::sortRows( cscJA, cscValues, numColumns, numRows, cscIA, testContext );

    //  compare with the CSC test data

    IndexType numRows1;
    IndexType numColumns1;
    IndexType numValues1;

    HArray<IndexType> expIA( testContext );
    HArray<IndexType> expJA( testContext );
    HArray<ValueType> expValues( testContext );

    data1::getCSCTestData( numRows1, numColumns1, numValues1, expIA, expJA, expValues );

    // test cscIA = expJA, cscJA = expIA, cscValues = expValues

    BOOST_TEST( hostReadAccess( expJA ) == hostReadAccess( cscIA ), per_element() );
    BOOST_TEST( hostReadAccess( expIA ) == hostReadAccess( cscJA ), per_element() );
    BOOST_TEST( hostReadAccess( expValues ) == hostReadAccess( cscValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_NUMERIC_TYPES_EXT_HOST> scai_ext_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( decompositionTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    if ( common::TypeTraits<IndexType>::stype != common::ScalarType::INT )
    {
        // decomposition external, requires IndexType = INT
        return;
    }

    //       array ( 4 x 4 )        sol      rhs
    //
    //       3    4   -5   6          1      39
    //       6    5   -6   5          2      43
    //       9   -4    2   3         -2       6
    //       -    2   -3   1          3      13

    const IndexType numRows = 4;

    HArray<IndexType> csrIA(  { 0, 4, 8, 12, 15 }, testContext );
    HArray<IndexType> csrJA( { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3 }, testContext );
    HArray<ValueType> csrValues( { 3.0,  4.0, -5.0,  6.0,
                                   6.0,  5.0, -6.0, 5.0,
                                   9.0, -4.0,  2.0, 3.0,
                                   2.0, -3.0, 1.0       },  testContext );

    HArray<ValueType> rhs( { 39.0, 43.0, 6.0, 13.0 }, testContext );
    HArray<ValueType> expSolution(  { 1.0, 2.0, -2.0, 3.0 }, testContext );

    HArray<ValueType> solution;

    CSRUtils::solve( solution, rhs, numRows, numRows, csrIA, csrJA, csrValues, false, testContext );

    auto maxDiff = HArrayUtils::maxDiffNorm( solution, expSolution );
    auto eps = common::TypeTraits<ValueType>::small();
  
    if ( maxDiff >= eps )
    {
        BOOST_TEST( hostReadAccess( solution ) == hostReadAccess( expSolution ), per_element() );
    }
    else
    {
        BOOST_CHECK( maxDiff < eps );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( matMulTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "matmul< " << TypeTraits<ValueType>::id() << "> non-square test" )

    // REMARK: test with explicit zeros because cusparse often has problem with it
    //       array 1             array 2
    //    1.0   -   2.0       1.0  0.5  0.0  4.0
    //    0.5  0.3   -         -   0.3  0.0  1.5
    //     -    -   3.0       2.0   -   3.0   -
    //    4.0  1.5   -

    const IndexType m  = 4;
    const IndexType k  = 3;
    const IndexType n  = 4;

    // csr arrays for matrix a

    HArray<IndexType> aIA(     { 0,        2,        4,   5,       7 }, testContext );
    HArray<IndexType> aJA(     { 0,   2,   0,   1,   2,   0,   1 }, testContext );
    HArray<ValueType> aValues( { 1.0, 2.0, 0.5, 0.3, 3.0, 4.0, 1.5 }, testContext );

    // csr arrays for matrix b

    HArray<IndexType> bIA(     { 0,                  4,             7,       9 }, testContext );
    HArray<IndexType> bJA(     { 0,   1,   2,   3,   1,   2,   3,   0,   2 }, testContext );
    HArray<ValueType> bValues( { 1.0, 0.5, 0.0, 4.0, 0.3, 0.0, 1.5, 2.0, 3.0 }, testContext );

    // REMARK: explicit zeros are also in result, when no compress is called
    //       array3 ( 4 x 4 )
    //
    //     5.0  0.5   6.0  4.0
    //     0.5  0.34  0.0  2.45
    //     6.0   -    9.0   -
    //     4.0  2.45  0.0  18.25

    HArray<IndexType> expCIA(    { 0, 4, 8, 10, 14 } ) ;
    HArray<IndexType> expCJA(     { 0,   1,   2,   3,   0,   1,    2,   3,    0,   2,   0,   1,    2,   3 } );
    HArray<ValueType> expCValues( { 5.0, 0.5, 6.0, 4.0, 0.5, 0.34, 0.0, 2.45, 6.0, 9.0, 4.0, 2.45, 0.0, 18.25 } );

    SCAI_LOG_INFO( logger, "matrixMultiplySizes< " << TypeTraits<ValueType>::id() << "> non-square test for " << *testContext )

    HArray<IndexType> cIA;
    HArray<IndexType> cJA;
    HArray<ValueType> cValues;
  
    ValueType alpha = 1;

    CSRUtils::matrixMultiply( cIA, cJA, cValues, alpha, aIA, aJA, aValues, bIA, bJA, bValues, m, n, k, testContext );
  
    BOOST_TEST( hostReadAccess( cIA ) == hostReadAccess( expCIA ), per_element() );

    // sort the columns, otherwise comparison might fail

    CSRUtils::sortRows( cJA, cValues, m, n, cIA, testContext );

    BOOST_TEST( hostReadAccess( cJA ) == hostReadAccess( expCJA ), per_element() );

    // check values, might not be exact

    auto eps = common::TypeTraits<ValueType>::small();

    BOOST_CHECK( utilskernel::HArrayUtils::maxDiffNorm( cValues, expCValues ) < eps );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( matAddTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    //       array 1                array 2
    //    1.0   -   2.0  -         1.0  0.5   -   -
    //    0.5  0.3   -   0.5        -    -    -   0.5
    //     -    -   3.0  -        2.0   -   1.0   0.5

    const IndexType m = 3;
    const IndexType n = 4;

    // csr arrays for matrix a

    HArray<IndexType> aIA(     { 0,        2,             5,   6 }, testContext );
    HArray<IndexType> aJA(     { 0,   2,   0,   1,   3,   2  }, testContext );
    HArray<ValueType> aValues( { 1.0, 2.0, 0.5, 0.3, 0.5, 3.0 }, testContext );

    // csr arrays for matrix b

    HArray<IndexType> bIA(     { 0,        2,   3,         6 }, testContext );
    HArray<IndexType> bJA(     { 0,   1,   3,   0,   2,  3   }, testContext );
    HArray<ValueType> bValues( { 1.0, 0.5, 0.5, 2.0, 1.0, 0.5 }, testContext );

    //       array3 = array 1 + array 2
    //
    //     2.0  0.5  2.0  -
    //     0.5  0.3   -   1.0
    //     2.0   -   4.0  0.5

    HArray<IndexType> expCIA(     { 0,             3,             6,            9 } );
    HArray<IndexType> expCJA(     { 0,   1,   2,   0,   1,   3,   0,   2,   3 } );
    HArray<ValueType> expCValues( { 2.0, 0.5, 2.0, 0.5, 0.3, 1.0, 2.0, 4.0, 0.5 } );

    HArray<IndexType> cIA;
    HArray<IndexType> cJA;
    HArray<ValueType> cValues;

    ValueType alpha = 1;
    ValueType beta  = 1;

    CSRUtils::matrixAdd( cIA, cJA, cValues, alpha, aIA, aJA, aValues, beta, bIA, bJA, bValues, m, n, testContext );

    // sort the columns, otherwise comparison might fail, no diagonals

    CSRUtils::sortRows( cJA, cValues, m, n, cIA, testContext );

    BOOST_TEST( hostReadAccess( cIA ) == hostReadAccess( expCIA ), per_element() );
    BOOST_TEST( hostReadAccess( cJA ) == hostReadAccess( expCJA ), per_element() );
    BOOST_TEST( hostReadAccess( cValues ) == hostReadAccess( expCValues ), per_element() );
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

    HArray<IndexType> aIA(     { 0,        2,             5,   6 }, testContext );
    HArray<IndexType> aJA(     { 0,   2,   0,   1,   3,   2  }, testContext );
    HArray<ValueType> aValues( { 1.0, 2.0, 0.5, 0.3, 0.5, 3.0 }, testContext );

    // csr arrays for matrix b, MUST BE sorted

    HArray<IndexType> bIA(     { 0,        2,   3,         6 }, testContext );
    HArray<IndexType> bJA(     { 0,   1,   3,   0,   2,  3   }, testContext );
    HArray<ValueType> bValues( { 1.0, 0.5, 0.5, 2.0, 1.0, 0.5 }, testContext );

    //       array3 = array 1 + array 2
    //
    //     2.0  0.5  2.0  -
    //     0.5  0.3   -   1.0
    //     2.0   -   4.0  0.5

    HArray<IndexType> expCIA(     { 0,             3,             6,            9 } );
    HArray<IndexType> expCJA(     { 0,   1,   2,   0,   1,   3,   0,   2,   3 } );
    HArray<ValueType> expCValues( { 2.0, 0.5, 2.0, 0.5, 0.3, 1.0, 2.0, 4.0, 0.5 } );

    HArray<IndexType> cIA;
    HArray<IndexType> cJA;
    HArray<ValueType> cValues;

    CSRUtils::binaryOp( cIA, cJA, cValues, aIA, aJA, aValues, bIA, bJA, bValues, m, n, common::BinaryOp::ADD, testContext );

    // Note: entries in storage C are SORTED

    BOOST_TEST( hostReadAccess( cIA ) == hostReadAccess( expCIA ), per_element() );
    BOOST_TEST( hostReadAccess( cJA ) == hostReadAccess( expCJA ), per_element() );
    BOOST_TEST( hostReadAccess( cValues ) == hostReadAccess( expCValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( reduceTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    ValueType initVal = 5;

    for ( IndexType dim = 0; dim < 2; ++dim )
    {
        // IMPORTANT: result array must already be set before

        IndexType nSize = dim == 0 ? numRows : numColumns ;

        HArray<ValueType> computedRes( nSize, initVal, testContext );

        // be careful, reduce updates its target array

        CSRUtils::reduce( computedRes, numRows, numColumns, csrIA, csrJA, csrValues,
                          dim, common::BinaryOp::ADD, common::UnaryOp::COPY, testContext );

        HArray<ValueType> expectedRes;

        data1::getReduceResult( expectedRes, dim );  // assumes ADD, SQR

        HArrayUtils::compute( expectedRes, expectedRes, common::BinaryOp::ADD, initVal );
    
        BOOST_TEST( hostReadAccess( computedRes ) == hostReadAccess( expectedRes ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    SCAI_ASSERT_EQ_ERROR( csrIA.size(), numRows + 1, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrJA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrValues.size(), numValues, "size mismatch" )

    HArray<ValueType> x( { 3, -3, 2, -2 }, testContext );
    HArray<ValueType> y( { 1, -1, 2, -2, 1, 1, -1 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numColumns, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, y.size(), "size mismatch" );

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

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * CSR * x + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", CSR: ia = " << csrIA << ", ja = " << csrJA << ", values = " << csrValues )

        common::MatrixOp op = common::MatrixOp::NORMAL;

        CSRUtils::gemv( res, alpha, x, beta, y, numRows, numColumns, csrIA, csrJA, csrValues, 
                        op, false, testContext );

        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTransTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    SCAI_ASSERT_EQ_ERROR( csrIA.size(), numRows + 1, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrJA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrValues.size(), numValues, "size mismatch" )

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

        SCAI_LOG_DEBUG( logger, "compute res = " << alpha << " * x * CSR + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", CSR: ia = " << csrIA << ", ja = " << csrJA << ", values = " << csrValues )

        auto op = common::MatrixOp::TRANSPOSE;

        CSRUtils::gemv( res, alpha, x, beta, y, numRows, numColumns, csrIA, csrJA, csrValues, 
                        op, false, testContext );

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEMVTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    data1::getRowIndexes( rowIndexes );

    SCAI_ASSERT_EQ_ERROR( csrIA.size(), numRows + 1, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrJA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrValues.size(), numValues, "size mismatch" )

    HArray<ValueType> x( { 3, -3, 2, -2 }, testContext );
    HArray<ValueType> y( { 1, -1, 2, -2, 1, 1, -1 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numColumns, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, y.size(), "size mismatch" );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];

        HArray<ValueType> res( y );

        SCAI_LOG_INFO( logger, "compute res += " << alpha << " * CSR * x "
                       << ", with x = " << x
                       << ", CSR: ia = " << csrIA << ", ja = " << csrJA << ", values = " << csrValues )

        auto op = common::MatrixOp::NORMAL;
        bool async = false;

        CSRUtils::gemvSp(  res, alpha, x, numRows, numColumns, csrIA, csrJA, csrValues, 
                           op, rowIndexes, async, testContext );

        ValueType beta = 1;  // res = alpha * A * x + 1 * res <-> res += alpha * A * x

        auto expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEVMTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    data1::getRowIndexes( rowIndexes );

    SCAI_ASSERT_EQ_ERROR( csrIA.size(), numRows + 1, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrJA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrValues.size(), numValues, "size mismatch" )

    const ValueType res_values[] = { 1, -1, 2, -2 };
    const ValueType x_values[]   = { 3, -2, -2, 3, 1, 0, 1 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_res = sizeof( res_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numRows, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, n_res, "size mismatch" );

    HArray<ValueType> x( numRows, x_values, testContext );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];

        HArray<ValueType> res( numColumns, res_values, testContext );

        SCAI_LOG_INFO( logger, "compute res += " << alpha << " * CSR * x "
                       << ", with x = " << x
                       << ", CSR: ia = " << csrIA << ", ja = " << csrJA << ", values = " << csrValues )

        auto op = common::MatrixOp::TRANSPOSE;

        CSRUtils::gemvSp(  res, alpha, x, numRows, numColumns, csrIA, csrJA, csrValues, 
                           op, rowIndexes, false, testContext );

        HArray<ValueType> expectedRes( numColumns, res_values );

        ValueType beta = 1;  // res = alpha * x * A + 1 * res <-> res += alpha * x * A

        expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, expectedRes );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemmSDTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    SCAI_ASSERT_EQ_ERROR( csrIA.size(), numRows + 1, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrJA.size(), numValues, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( csrValues.size(), numValues, "size mismatch" )

    // we compare gemmSD with two individual gemv operations

    const IndexType nVectors = 2;  // gemm is optimized for nVectors * gemv

    HArray<ValueType> x( { 3, -3, 2, -2,
                           2, -1, 1,  2 }, testContext );

    HArray<ValueType> y( { 1, -1, 2, -2, 1, 1, -1 ,
                           2, -3, 4,  1, 5, 3,  2  }, testContext );

    SCAI_ASSERT_EQ_ERROR( numColumns * nVectors, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows * nVectors, y.size(), "size mismatch" );

    // x and y are two vectors, but for verification we need the single ones

    HArray<ValueType> x1( numColumns );  // becomes x[ :, 0 ]
    HArray<ValueType> x2( numColumns );  // becomes x[ :, 1
    HArray<ValueType> y1( numRows ) ;    // becomes y[ :, 0 ]
    HArray<ValueType> y2( numRows );     // becomes y[ :, 1 ]

    HArrayUtils::setArraySection( y1, 0, 1, y, 0, nVectors, numRows );
    HArrayUtils::setArraySection( y2, 0, 1, y, 1, nVectors, numRows );
    HArrayUtils::setArraySection( x1, 0, 1, x, 0, nVectors, numColumns );
    HArrayUtils::setArraySection( x2, 0, 1, x, 1, nVectors, numColumns );

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

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * CSR * x + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", CSR: ia = " << csrIA << ", ja = " << csrJA << ", values = " << csrValues )

        bool async = false;

        auto op = common::MatrixOp::NORMAL;

        CSRUtils::gemmSD( res,  alpha, x, beta, y,
                          numRows, numColumns, nVectors, csrIA, csrJA, csrValues, op, async, testContext );

        HArray<ValueType> expectedRes1 = data1::getGEMVNormalResult( alpha, x1, beta, y1 );
        HArray<ValueType> expectedRes2 = data1::getGEMVNormalResult( alpha, x2, beta, y2 );

        {
            auto rComputed = hostReadAccess( res );
            auto rExpected1 = hostReadAccess( expectedRes1 );
            auto rExpected2 = hostReadAccess( expectedRes2 );

            for ( IndexType i = 0; i < numRows; ++i )
            {
                BOOST_CHECK_EQUAL( rExpected1[i], rComputed[ nVectors * i + 0 ] );
                BOOST_CHECK_EQUAL( rExpected2[i], rComputed[ nVectors * i + 1 ] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemmDSTest, ValueType, scai_numeric_test_types )
{
    /*   DenseMatrix = alpha * DenseMatrix * SparseMatrix
   
           1  -1  2  -2  1  1  -1       6  0  0  4           -3  -11  18   7
           2  -3  4   1  5  3   2   *   7  0  0  0     =      3    7  36  36 
           1   1  0  -1  2  1   3       0  0  9  4           15   -2   0   9
                                        2  5  0  3 
                                        2  0  0  1  
                                        0  0  0  0  
                                        0  1  0  2  
    */

    ContextPtr testContext = ContextFix::testContext;

    const IndexType numRows = 7;
    const IndexType numColumns = 4;

    HArray<IndexType> csrIA(     { 0,    2, 3,    5,       8,    10, 10,   12 }, testContext );
    HArray<IndexType> csrJA(     { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3,     1, 3 }, testContext );
    HArray<ValueType> csrValues( { 6, 4, 7, 9, 4, 2, 5, 3, 2, 1,     1, 2 }, testContext );

    const IndexType nVector = 3;
    
    HArray<ValueType> result( { 0,  0, 0,  0,
                                0,  0, 0,  0,
                                0,  0, 0,  0  }, testContext );

    HArray<ValueType> expResult( { -3, -11, 18, 5,
                                    3, 7, 36, 36,
                                   15, -2, 0,  9  }, testContext );

    HArray<ValueType> x( { 1, -1, 2, -2, 1, 1, -1 ,
                           2, -3, 4,  1, 5, 3,  2 ,
                           1,  1, 0, -1, 2, 1,  3  }, testContext );

    SCAI_ASSERT_EQ_ERROR( nVector * numColumns, result.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( nVector * numRows, x.size(), "size mismatch" );

    bool async = false;

    auto op = common::MatrixOp::NORMAL;

    ValueType alpha = 1;
    ValueType beta  = 1;

    CSRUtils::gemmDS( result,  alpha, x, beta,
                      numRows, numColumns, nVector, csrIA, csrJA, csrValues, op, async, testContext );

    SCAI_CHECK_EQUAL_ARRAY( result, expResult )

    alpha = 2;
    beta  = -2;

    std::shared_ptr<tasking::SyncToken> token;

    async = true;
    token.reset( CSRUtils::gemmDS( result,  alpha, x, beta,
                                   numRows, numColumns, nVector, 
                                   csrIA, csrJA, csrValues, op, async, testContext ) );

    token->wait();

    // result should be a zero vector

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    BOOST_CHECK( HArrayUtils::maxNorm( result ) < eps );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiTest, ValueType, scai_numeric_test_types )
{
    // Run one jacobi iteration step with a CSR storage

    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext )

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data2::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    HArray<ValueType> rhs( { 1, -1, 2, -2 }, testContext );
    HArray<ValueType> oldSolution( { 3, -2, -2, 3 }, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> newSolution( testContext );

        CSRUtils::jacobi( newSolution, omega, oldSolution, rhs, csrIA, csrJA, csrValues, testContext );

        HArray<ValueType> expSolution;

        data2::getJacobiResult( expSolution, oldSolution, omega, rhs );

        auto eps = common::TypeTraits<ValueType>::small();

        auto maxDiff = HArrayUtils::maxDiffNorm( expSolution, newSolution );

        BOOST_CHECK( maxDiff < eps );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiHaloTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "jacobiHalo test @ " << *testContext )

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    data1::getRowIndexes( rowIndexes );

    HArray<ValueType> oldSolution( { 3, -2, -2, 3, 1, 0, 2 }, testContext );
    HArray<ValueType> diag( { 9,  8,  7, 6, 7, 8, 9 }, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    // create a dummy array for the offsets of the local part

    HArray<IndexType> iaDummy;
    HArrayUtils::setOrder( iaDummy, numRows + 1, testContext );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> solution( numRows, ValueType( 0 ), testContext );

        CSRUtils::jacobiHalo( solution, omega, diag, oldSolution,
                              csrIA, csrJA, csrValues, rowIndexes, testContext );

        HArray<ValueType> expectedSol( numRows, ValueType( 0 ), testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        auto maxDiff = HArrayUtils::maxDiffNorm( expectedSol, solution );

        if ( maxDiff >= common::TypeTraits<ValueType>::small() )
        {
            BOOST_TEST( hostReadAccess( expectedSol ) == hostReadAccess( oldSolution ), per_element() );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()
