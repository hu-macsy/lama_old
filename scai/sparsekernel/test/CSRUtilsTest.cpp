/**
 * @file CSRUtilsTest.cpp
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
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/test/TestMacros.hpp>

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

BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxDiffValTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;

    ContextPtr loc = testContext;
    absMaxDiffVal.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "absMaxDiffVal< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    // input arrays
    //    Array1             Array2
    //
    //    1 2 3 0 0          1 2 0 0 0
    //    0 0 1 1 2          1 0 2 2 1
    const IndexType ia1[] =
    { 0, 3, 6 };
    const IndexType ja1[] =
    { 0, 1, 2, 2, 3, 4 };
    const IndexType ia2[] =
    { 0, 2, 6 };
    const IndexType ja2[] =
    { 0, 1, 0, 2, 3, 4 };
    const ValueType values1[] =
    { 1, 2, 3, 1, 1, 2 };
    const ValueType values2[] =
    { 1, 2, 1, 2, 2, 1 };
    const IndexType numRows = 2;
    // const IndexType numColumns = 5;
    const IndexType numValues1 = sizeof( ja1 ) / sizeof( IndexType );
    const IndexType numValues2 = sizeof( ja2 ) / sizeof( IndexType );
    HArray<IndexType> csrIA1( numRows + 1, ia1, testContext );
    HArray<IndexType> csrJA1( numValues1, ja1, testContext );
    HArray<ValueType> csrValues1( numValues1, values1, testContext );
    HArray<IndexType> csrIA2( numRows + 1, ia2, testContext );
    HArray<IndexType> csrJA2( numValues2, ja2, testContext );
    HArray<ValueType> csrValues2( numValues2, values2, testContext );
    ReadAccess<IndexType> rCSRIA1( csrIA1, loc );
    ReadAccess<IndexType> rCSRJA1( csrJA1, loc );
    ReadAccess<ValueType> rCSRValues1( csrValues1, loc );
    ReadAccess<IndexType> rCSRIA2( csrIA2, loc );
    ReadAccess<IndexType> rCSRJA2( csrJA2, loc );
    ReadAccess<ValueType> rCSRValues2( csrValues2, loc );
    SCAI_CONTEXT_ACCESS( loc );
    ValueType maxVal = absMaxDiffVal[loc]( numRows, false, rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), rCSRIA2.get(),
                                           rCSRJA2.get(), rCSRValues2.get() );
    BOOST_CHECK_EQUAL( 3, maxVal );
    // rows are sorted, so we can also apply sortFlag = true
    maxVal = absMaxDiffVal[loc]( numRows, true, rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), rCSRIA2.get(),
                                 rCSRJA2.get(), rCSRValues2.get() );
    BOOST_CHECK_EQUAL( 3, maxVal );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::scaleRows<ValueType> > scaleRows;

    ContextPtr loc = testContext;

    scaleRows.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "scaleRows test for " << *testContext << " on " << *loc )

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    HArray<ValueType> savedValues( csrValues );  // keep a copy for comparison later

    const ValueType row_factors[]   = { 2, 3, 4, 5, 1, 3, 2 };

    const IndexType n_factors = sizeof( row_factors ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( numRows, n_factors );

    HArray<ValueType> rows( n_factors, row_factors, testContext );

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<ValueType> wValues( csrValues, loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        ReadAccess<ValueType> rRows( rows, loc );

        scaleRows[loc]( wValues.get(), rIA.get(), numRows, rRows.get() );
    }

    // prove by hand

    {
        ReadAccess<IndexType> rIA( csrIA );
        ReadAccess<IndexType> rJA( csrJA );

        ReadAccess<ValueType> rRows( rows );
        ReadAccess<ValueType> rSavedValues( savedValues );
        ReadAccess<ValueType> rValues( csrValues );

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

    LAMAKernel<CSRKernelTrait::compress<ValueType> > compress;

    ContextPtr loc = testContext;
    compress.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

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

    LAMAKernel<CSRKernelTrait::decomposition<ValueType> > decomposition;

    ContextPtr loc = testContext;
    decomposition.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected
    SCAI_LOG_INFO( logger, "decomposition< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

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

    const IndexType ia[]     = { 0, 4, 8, 12, 15 };
    const IndexType ja[]     = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3 };
    const ValueType values[] = { 3.0,  4.0, -5.0,  6.0,
                                 6.0,  5.0, -6.0, 5.0,
                                 9.0, -4.0,  2.0, 3.0,
                                 2.0, -3.0, 1.0
                               };
    const ValueType rhsValues[] = { 39.0, 43.0, 6.0, 13.0 };
    const ValueType solValues[] = { 1.0, 2.0, -2.0, 3.0 };
    const IndexType numRows = 4;
    const IndexType nnz = 15;

    HArray<IndexType> csrIA( numRows + 1, ia, testContext );
    HArray<IndexType> csrJA( nnz, ja, testContext );
    HArray<ValueType> csrValues( nnz, values, testContext );

    HArray<ValueType> rhs( numRows, rhsValues, testContext );
    HArray<ValueType> solution;

    {
        ReadAccess<IndexType> rCSRIA( csrIA, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        ReadAccess<ValueType> rCSRValues( csrValues, loc );
        ReadAccess<ValueType> rRHS( rhs, loc );
        WriteOnlyAccess<ValueType> wSol( solution, loc, numRows );
        SCAI_CONTEXT_ACCESS( loc );
        decomposition[loc]( wSol.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(),
                            rRHS.get(), numRows, nnz, false );
    }

    {
        ContextPtr host = Context::getHostPtr();
        ReadAccess<ValueType> rSol( solution, host );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType x = rSol[i] - solValues[i];
            BOOST_CHECK_SMALL( common::Math::real( x ), common::TypeTraits<ValueType>::small() );
            BOOST_CHECK_SMALL( common::Math::imag( x ), common::TypeTraits<ValueType>::small() );
        }
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
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::reduce<ValueType> > reduce;

    ContextPtr loc = testContext;

    reduce.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    for ( IndexType dim = 0; dim < 2; ++dim )
    {
        HArray<ValueType> computedRes( loc );

        {
            SCAI_CONTEXT_ACCESS( loc );

            IndexType resSize = dim == 0 ? numRows : numColumns;

            computedRes.setSameValue( resSize, ValueType( 0 ) );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );
    
            WriteAccess<ValueType> wResult( computedRes, loc );

            reduce[loc]( wResult.get(),
                         rIA.get(), rJA.get(), rValues.get(), numRows, dim,
                         common::BinaryOp::ADD, common::UnaryOp::COPY );
        }

        HArray<ValueType> expectedRes;

        data1::getReduceResult( expectedRes, dim );  // assumes ADD, SQR
    
        BOOST_CHECK_EQUAL( computedRes.size(), expectedRes.size() );

        {
            ReadAccess<ValueType> rComputed( computedRes, hostContext );
            ReadAccess<ValueType> rExpected( expectedRes, hostContext );
    
            for ( IndexType i = 0; i < computedRes.size(); ++i )
            {
                BOOST_CHECK_EQUAL( rExpected[i], rComputed[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

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
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );
            WriteOnlyAccess<ValueType> wResult( res, loc, numRows );

            common::MatrixOp op = common::MatrixOp::NORMAL;

            if ( beta == 0 ) 
            {
                normalGEMV[loc]( wResult.get(),
                                 alpha, rX.get(), beta, NULL,
                                 numRows, numColumns, numValues, rIA.get(), rJA.get(), rValues.get(), op );
            }
            else
            {
                normalGEMV[loc]( wResult.get(),
                                 alpha, rX.get(), beta, rY.get(),
                                 numRows, numColumns, numValues, rIA.get(), rJA.get(), rValues.get(), op );
            }
        }

        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        {
            ReadAccess<ValueType> rComputed( res, hostContext );
            ReadAccess<ValueType> rExpected( expectedRes, hostContext );

            for ( IndexType i = 0; i < numRows; ++i )
            {
                BOOST_CHECK_EQUAL( rExpected[i], rComputed[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTransTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

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

    const ValueType y_values[]   = { 1, -1, 2, -2 };
    const ValueType x_values[]   = { 3, -2, -2, 3, 1, 0, 1 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numRows, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, n_y, "size mismatch" );

    HArray<ValueType> x( numRows, x_values, testContext );
    HArray<ValueType> y( numColumns, y_values, testContext );

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

        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );
            WriteOnlyAccess<ValueType> wResult( res, loc, numColumns );

            common::MatrixOp op = common::MatrixOp::TRANSPOSE;

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, rY.get(),
                             numRows, numColumns, numValues, rIA.get(), rJA.get(), rValues.get(), op );
        }

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEMVTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::sparseGEMV<ValueType> > sparseGEMV;

    ContextPtr loc = testContext;

    sparseGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    data1::getRowIndexes( rowIndexes );
    IndexType numNonEmptyRows = rowIndexes.size();

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
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );
            ReadAccess<IndexType> rIndexes( rowIndexes, loc );

            ReadAccess<ValueType> rX( x, loc );
            WriteAccess<ValueType> wResult( res, loc );

            auto op = common::MatrixOp::NORMAL;

            sparseGEMV[loc]( wResult.get(),
                             alpha, rX.get(),
                             numNonEmptyRows, rIndexes.get(), rIA.get(), rJA.get(), rValues.get(), op );

        }

        ValueType beta = 1;  // res = alpha * A * x + 1 * res <-> res += alpha * A * x

        auto expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEVMTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::sparseGEMV<ValueType> > sparseGEMV;

    ContextPtr loc = testContext;

    sparseGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;
    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    data1::getRowIndexes( rowIndexes );
    IndexType numNonEmptyRows = rowIndexes.size();

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
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );
            ReadAccess<IndexType> rIndexes( rowIndexes, loc );

            ReadAccess<ValueType> rX( x, loc );
            WriteAccess<ValueType> wResult( res, loc );

            sparseGEMV[loc]( wResult.get(),
                             alpha, rX.get(),
                             numNonEmptyRows, rIndexes.get(), rIA.get(), rJA.get(), rValues.get(), common::MatrixOp::TRANSPOSE );

        }

        HArray<ValueType> expectedRes( numColumns, res_values );

        ValueType beta = 1;  // res = alpha * x * A + 1 * res <-> res += alpha * x * A

        expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, expectedRes );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemmTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::gemm<ValueType> > gemm;

    ContextPtr loc = testContext;

    gemm.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

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

    // we compare gemm with two individual gemv operations

    const IndexType nVectors = 2;  // gemm is optimized for nVectors * gemv

    const ValueType y_values[]   = { 1, -1, 2, -2, 1, 1, -1 ,
                                     2, -3, 4,  1, 5, 3,  2
                                   };
    const ValueType x_values[]   = { 3, -3, 2, -2,
                                     2, -1, 1,  2
                                   };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numColumns * nVectors, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows * nVectors, n_y, "size mismatch" );

    HArray<ValueType> x( n_x, x_values, testContext );
    HArray<ValueType> y( n_y, y_values, testContext );

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
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );

            WriteOnlyAccess<ValueType> wResult( res, loc, nVectors * numRows );

            gemm[loc]( wResult.get(),
                       alpha, rX.get(), beta, rY.get(),
                       numRows, nVectors, numColumns, rIA.get(), rJA.get(), rValues.get() );
        }

        HArray<ValueType> expectedRes1 = data1::getGEMVNormalResult( alpha, x1, beta, y1 );
        HArray<ValueType> expectedRes2 = data1::getGEMVNormalResult( alpha, x2, beta, y2 );

        {
            ReadAccess<ValueType> rComputed( res, hostContext );
            ReadAccess<ValueType> rExpected1( expectedRes1, hostContext );
            ReadAccess<ValueType> rExpected2( expectedRes2, hostContext );

            for ( IndexType i = 0; i < numRows; ++i )
            {
                BOOST_CHECK_EQUAL( rExpected1[i], rComputed[ nVectors * i + 0 ] );
                BOOST_CHECK_EQUAL( rExpected2[i], rComputed[ nVectors * i + 1 ] );
            }
        }
    }
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
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = testContext;

    jacobiHalo.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    SCAI_LOG_INFO( logger, "jacobiHalo test for " << *testContext << " on " << *loc )

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    data1::getRowIndexes( rowIndexes );
    IndexType numNonEmptyRows = rowIndexes.size();

    const ValueType old_values[]   = { 3, -2, -2, 3, 1, 0, 2 };
    const ValueType diag_values[]  = { 9,  8,  7, 6, 7, 8, 9 };

    const IndexType n_old_values = sizeof( old_values ) / sizeof( ValueType );
    SCAI_ASSERT_EQ_ERROR( n_old_values, numRows, "size mismatch" );

    HArray<ValueType> oldSolution( numRows, old_values, testContext );
    HArray<ValueType> diag( numRows, diag_values, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    // create a dummy array for the offsets of the local part

    HArray<IndexType> iaDummy;
    HArrayUtils::setOrder( iaDummy, numRows + 1, testContext );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> solution( numRows, ValueType( 0 ), testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );
            ReadAccess<IndexType> rIndexes( rowIndexes, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<IndexType> rIADummy( iaDummy, loc );
            ReadAccess<ValueType> rDiag( diag, loc );
            WriteAccess<ValueType> wSolution( solution, loc, numColumns );

            jacobiHalo[loc]( wSolution.get(), rDiag.get(),
                             rIA.get(), rJA.get(), rValues.get(), 
                             rIndexes.get(), rOld.get(), omega, numNonEmptyRows );
        }

        HArray<ValueType> expectedSol( numRows, ValueType( 0 ), testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        auto maxDiff = HArrayUtils::maxDiffNorm( expectedSol, solution );

        BOOST_CHECK( maxDiff < 0.1 );

        bool mustBeIdentical = false;

        if ( mustBeIdentical )
        {
            ReadAccess<ValueType> rExpected( expectedSol );
            ReadAccess<ValueType> rComputed( solution );

            for ( IndexType i = 0; i < numRows; ++i )
            {
                BOOST_CHECK_EQUAL( rExpected[i], rComputed[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiHaloDiagTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = testContext;

    jacobiHalo.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    SCAI_LOG_INFO( logger, "jacobiHalo test for " << *testContext << " on " << *loc )

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;
    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    data1::getRowIndexes( rowIndexes );
    IndexType numNonEmptyRows = rowIndexes.size();

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

        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );
            ReadAccess<IndexType> rIndexes( rowIndexes, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<ValueType> rDiag( diag, loc );
            WriteAccess<ValueType> wSolution( solution, loc, numColumns );

            jacobiHalo[loc]( wSolution.get(),
                             rDiag.get(),
                             rIA.get(), rJA.get(), rValues.get(), rIndexes.get(),
                             rOld.get(), omega, numNonEmptyRows );
        }

        HArray<ValueType> expectedSol( numRows, ValueType( 0 ), testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        auto maxDiff = HArrayUtils::maxDiffNorm( expectedSol, solution );

        BOOST_CHECK( maxDiff < 0.1 );

        bool mustBeIdentical = false;

        if ( mustBeIdentical )
        {
            ReadAccess<ValueType> rExpected( expectedSol );
            ReadAccess<ValueType> rComputed( solution );

            for ( IndexType i = 0; i < numRows; ++i )
            {
                BOOST_CHECK_EQUAL( rExpected[i], rComputed[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()
