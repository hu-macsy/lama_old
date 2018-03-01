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
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
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

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CSRUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( validOffsetsTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::validOffsets> validOffsets;

    ContextPtr loc = testContext;

    validOffsets.getSupportedContext( loc );

    const IndexType ia[]   = { 0, 2, 5, 7 };
    const IndexType ia_f[] = { 0, 5, 2, 7 };

    IndexType numRows = sizeof( ia ) / sizeof( IndexType ) - 1;
    IndexType total   = ia[numRows];

    HArray<IndexType> csrIA( numRows + 1, ia, testContext );

    bool okay = false;

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        okay = validOffsets[loc]( rIA.get(), numRows, total );
    }

    BOOST_CHECK( okay );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        okay = validOffsets[loc]( rIA.get(), numRows - 1, total );
    }

    BOOST_CHECK( !okay );

    csrIA.setRawData( numRows + 1, ia_f );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        okay = validOffsets[loc]( rIA.get(), numRows, total );
    }

    BOOST_CHECK( !okay );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( offsets2sizesTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::offsets2sizes> offsets2sizes;

    ContextPtr loc = testContext;
    offsets2sizes.getSupportedContext( loc );

    const IndexType offsets[] = { 0, 5, 11, 16, 19, 19, 21, 23 };
    const IndexType sizes[]   = { 5, 6,  5,  3,  0,  2,  2 };

    IndexType numRows = sizeof( offsets ) / sizeof( IndexType ) - 1;

    HArray<IndexType> csrIA( numRows + 1, offsets, testContext );

    HArray<IndexType> csrSizes( testContext );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        WriteOnlyAccess<IndexType> wSizes( csrSizes, loc, numRows );
        offsets2sizes[loc]( wSizes.get(), rIA.get(), numRows );
    }

    {
        ReadAccess<IndexType> rSizes( csrSizes );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( sizes[i], rSizes[i] );
        }
    }

    LAMAKernel<CSRKernelTrait::offsets2sizesGather> offsets2sizesGather;

    loc = testContext;
    offsets2sizesGather.getSupportedContext( loc );

    const IndexType row_indexes[] = { 0, 2, 4, 6 };

    IndexType n = sizeof( row_indexes ) / sizeof( IndexType );

    HArray<IndexType> rows( n, row_indexes, testContext );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        ReadAccess<IndexType> rRows( rows, loc );
        WriteOnlyAccess<IndexType> wSizes( csrSizes, loc, n );
        offsets2sizesGather[loc]( wSizes.get(), rIA.get(), rRows.get(), n );
    }

    {
        ReadAccess<IndexType> rSizes( csrSizes );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( sizes[row_indexes[i]], rSizes[i] );
        }
    }
}


/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( nonEmptyRowsTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::countNonEmptyRowsByOffsets> countNonEmptyRowsByOffsets;

    ContextPtr loc = testContext;
    countNonEmptyRowsByOffsets.getSupportedContext( loc );

    const IndexType ia[]   = { 0, 2, 2, 5, 5, 7, 7, 9 };
    const IndexType rows[] = { 0,    2,    4,    6  };

    IndexType numRows = sizeof( ia ) / sizeof( IndexType ) - 1;

    HArray<IndexType> csrIA( numRows + 1, ia, testContext );

    IndexType nonEmptyRows = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        nonEmptyRows = countNonEmptyRowsByOffsets[loc]( rIA.get(), numRows );
    }

    const IndexType expectedNonEmptyRows = sizeof( rows ) / sizeof( IndexType );

    BOOST_REQUIRE_EQUAL( expectedNonEmptyRows, nonEmptyRows );

    HArray<IndexType> rowIndexes( testContext );

    LAMAKernel<CSRKernelTrait::setNonEmptyRowsByOffsets> setNonEmptyRowsByOffsets;

    loc = testContext;
    setNonEmptyRowsByOffsets.getSupportedContext( loc );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        WriteOnlyAccess<IndexType> wIndexes( rowIndexes, loc, nonEmptyRows );
        setNonEmptyRowsByOffsets[loc]( wIndexes.get(), nonEmptyRows, rIA.get(), numRows );
    }

    {
        ReadAccess<IndexType> rIndexes( rowIndexes );

        for ( IndexType i = 0; i < nonEmptyRows; ++i )
        {
            BOOST_CHECK_EQUAL( rows[i], rIndexes[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( hasDiagonalPropertyTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<CSRKernelTrait::hasDiagonalProperty> hasDiagonalProperty;

    ContextPtr loc = testContext;

    hasDiagonalProperty.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    LArray<IndexType> csrIA( testContext );
    LArray<IndexType> csrJA( testContext );
    LArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    bool okay;

    IndexType numDiagonals = common::Math::min( numRows, numColumns );

    {
        ReadAccess<IndexType> rIA( csrIA, loc );
        ReadAccess<IndexType> rJA( csrJA, loc );

        SCAI_CONTEXT_ACCESS( loc );

        okay = hasDiagonalProperty[loc]( numDiagonals, rIA.get(), rJA.get() );
    }

    BOOST_CHECK( !okay );

    // data set 2 has a square matrix with diagonal entries first

    data2::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    BOOST_REQUIRE_EQUAL( numRows, numColumns );

    {
        ReadAccess<IndexType> rIA( csrIA, loc );
        ReadAccess<IndexType> rJA( csrJA, loc );

        SCAI_CONTEXT_ACCESS( loc );

        okay = hasDiagonalProperty[loc]( numRows, rIA.get(), rJA.get() );
    }

    BOOST_CHECK( okay );
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

    LAMAKernel<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;

    ContextPtr loc = testContext;
    convertCSR2CSC.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected
    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    //    1.0   -   2.0       1.0  0.5   -
    //    0.5  0.3   -         -   0.3   -
    //     -    -   3.0       2.0   -   3.0
    const IndexType ia1[] =
    { 0, 2, 4, 5 };
    const IndexType ja1[] =
    { 0, 2, 0, 1, 2 };
    const IndexType ia2[] =
    { 0, 2, 3, 5 };
    const IndexType ja2[] =
    { 0, 1, 1, 0, 2 };
    const ValueType values1[] =
    { 1.0, 2.0, 0.5, 0.3, 3.0 };
    const ValueType values2[] =
    { 1.0, 0.5, 0.3, 2.0, 3.0 };
    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType numValues = 5;
    HArray<IndexType> csrIA( numRows + 1, ia1, testContext );
    HArray<IndexType> csrJA( numValues, ja1, testContext );
    HArray<ValueType> csrValues( numValues, values1, testContext );
    HArray<IndexType> cscIA;
    HArray<IndexType> cscJA;
    HArray<ValueType> cscValues;
    {
        ReadAccess<IndexType> rCSRIA( csrIA, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        ReadAccess<ValueType> rCSRValues( csrValues, loc );
        WriteOnlyAccess<IndexType> wCSCIA( cscIA, loc, numColumns + 1 );
        WriteOnlyAccess<IndexType> wCSCJA( cscJA, loc, numValues );
        WriteOnlyAccess<ValueType> wCSCValues( cscValues, loc, numValues );
        SCAI_CONTEXT_ACCESS( loc );
        convertCSR2CSC[loc]( wCSCIA.get(), wCSCJA.get(), wCSCValues.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(), numRows,
                             numColumns, numValues );
    }
    {
        ReadAccess<IndexType> rCSCIA( cscIA );
        WriteAccess<IndexType> wCSCJA( cscJA );
        WriteAccess<ValueType> wCSCValues( cscValues );

        for ( IndexType j = 0; j <= numColumns; ++j )
        {
            BOOST_CHECK_EQUAL( rCSCIA[j], ia2[j] );
        }

        // For comparison of cscJA and cscValue we need to sort it
        bool diagonalFlag = false;
        OpenMPCSRUtils::sortRowElements( wCSCJA.get(), wCSCValues.get(), rCSCIA.get(), numColumns, diagonalFlag );

        for ( IndexType j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCJA[j], ja2[j] );
        }

        for ( IndexType j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( wCSCValues[j], values2[j] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( sortRowTest )
{
    typedef double ValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;

    ContextPtr loc = testContext;
    sortRowElements.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    //    1.0   -   2.0  1.1
    //    0.5  0.3   -    -
    //     -    -   3.0   -
    //     -    -   4.0   1.5

    const IndexType ia[] = { 0, 3, 5, 6, 8 };

    const IndexType ja_unsorted[] = { 3, 2, 0, 0, 1, 2, 3, 2 };
    const ValueType values_unsorted[] = { 1.1, 2.0, 1.0, 0.5, 0.3 , 3.0, 1.5, 4.0 };

    const IndexType ja_sorted[] = { 0, 2, 3, 0, 1, 2, 2, 3 };
    const ValueType values_sorted[] = { 1.0, 2.0, 1.1, 0.5, 0.3, 3.0, 4.0, 1.5 };

    const IndexType ja_dia_sorted[] = { 0, 2, 3, 1, 0, 2, 3, 2 };
    const ValueType values_dia_sorted[] = { 1.0, 2.0, 1.1, 0.3, 0.5, 3.0, 1.5, 4.0 };

    const IndexType numRows = 4;
    const IndexType numValues = 8;

    HArray<IndexType> csrIA( numRows + 1, ia, testContext );
    HArray<IndexType> csrJA( numValues, ja_unsorted, testContext );
    HArray<double> csrValues( numValues, values_unsorted, testContext );

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( csrIA, loc );
        WriteAccess<IndexType> wJA( csrJA, loc );
        WriteAccess<ValueType> wValues( csrValues, loc );

        bool diagonalFlag = false;

        sortRowElements[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, diagonalFlag );
    }

    // now check that values are sorted
    {
        ReadAccess<IndexType> rJA( csrJA );
        ReadAccess<ValueType> rValues( csrValues );

        for ( IndexType k = 0; k < numValues; ++k )
        {
            BOOST_CHECK_EQUAL( rJA[k], ja_sorted[k] );
            BOOST_CHECK_EQUAL( rValues[k], values_sorted[k] );
        }
    }


    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( csrIA, loc );
        WriteAccess<IndexType> wJA( csrJA, loc );
        WriteAccess<ValueType> wValues( csrValues, loc );

        bool diagonalFlag = true;

        sortRowElements[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, diagonalFlag );
    }

    // now check that values are sorted
    {
        ReadAccess<IndexType> rJA( csrJA );
        ReadAccess<ValueType> rValues( csrValues );

        for ( IndexType k = 0; k < numValues; ++k )
        {
            BOOST_CHECK_EQUAL( rJA[k], ja_dia_sorted[k] );
            BOOST_CHECK_EQUAL( rValues[k], values_dia_sorted[k] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( countNonZeroTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::countNonZeros<ValueType> > countNonZeros;

    ContextPtr loc = testContext;
    countNonZeros.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    //    1.0   -   2.0  1.1
    //    0.5  0.0   -    -
    //     -    -   3.0   -
    //     -    -   4.0  0.0

    const IndexType ia_raw[]     = {   0,             3,        5,   6,       8 };
    const IndexType ja_raw[]     = {   3,   2,   0,   0,   1,   2,   3,   2 };
    const ValueType values_raw[] = { 1.1, 2.0, 1.0, 0.5, 0.0, 3.0, 0.0, 4.0 };

    const IndexType compress_sizes[]  = {  3, 1, 1, 1 };

    const IndexType numRows = 4;
    const IndexType numValues = 8;

    HArray<IndexType> ia( numRows + 1, ia_raw, testContext );
    HArray<IndexType> ja( numValues, ja_raw, testContext );
    HArray<ValueType> values( numValues, values_raw, testContext );

    HArray<IndexType> sizes;

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( ia, loc );
        ReadAccess<IndexType> rJA( ja, loc );
        ReadAccess<ValueType> rValues( values, loc );

        WriteOnlyAccess<IndexType> wSizes( sizes, loc, numRows );

        bool diagonalFlag = false;
        ValueType eps = 0.000001;

        countNonZeros[loc]( wSizes.get(), rIA.get(), rJA.get(), rValues.get(),
                            numRows, eps, diagonalFlag );
    }

    // now check the values of sizes
    {
        ReadAccess<IndexType> rSizes( sizes );

        for ( IndexType k = 0; k < numRows; ++k )
        {
            BOOST_CHECK_EQUAL( rSizes[k], compress_sizes[k] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( compressTest )
{
    typedef double ValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::compress<ValueType> > compress;

    ContextPtr loc = testContext;
    compress.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    //    1.0   -   2.0  1.1
    //    0.5  0.0   -    -
    //     -    -   3.0   -
    //     -    -   4.0  0.0

    const IndexType ia_original[]     = {   0,             3,        5,   6,       8 };
    const IndexType ja_original[]     = {   3,   2,   0,   0,   1,   2,   3,   2 };
    const ValueType values_original[] = { 1.1, 2.0, 1.0, 0.5, 0.0, 3.0, 0.0, 4.0 };

    const IndexType ia_compress[]     = {   0,             3,   4,   5,   6 };
    const IndexType ja_compress[]     = {   3,   2,   0,   0,   2,   2 };
    const ValueType values_compress[] = { 1.1, 2.0, 1.0, 0.5, 3.0, 4.0 };

    const IndexType numRows = 4;
    const IndexType oldNumValues = 8;
    const IndexType newNumValues = 6;

    HArray<IndexType> oldIA( numRows + 1, ia_original, testContext );
    HArray<IndexType> oldJA( oldNumValues, ja_original, testContext );
    HArray<double>    oldValues( oldNumValues, values_original, testContext );

    HArray<IndexType> newIA( numRows + 1, ia_compress, testContext );
    HArray<IndexType> newJA;
    HArray<double>    newValues;

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( newIA, loc );

        WriteOnlyAccess<IndexType> wJA( newJA, loc, newNumValues );
        WriteOnlyAccess<ValueType> wValues( newValues, loc, newNumValues );

        ReadAccess<IndexType> roIA( oldIA, loc );
        ReadAccess<IndexType> roJA( oldJA, loc );
        ReadAccess<ValueType> roValues( oldValues, loc );

        bool diagonalFlag = false;
        double eps = 0.000001;

        compress[loc]( wJA.get(), wValues.get(), rIA.get(), roIA.get(), roJA.get(), roValues.get(), numRows, eps, diagonalFlag );
    }

    // now check that non-zero values are at the right place
    {
        ReadAccess<IndexType> rJA( newJA );
        ReadAccess<ValueType> rValues( newValues );

        for ( IndexType k = 0; k < newNumValues; ++k )
        {
            SCAI_LOG_TRACE( logger, "value " << k << ": " << rJA[k] << ":" << rValues[k]
                            << ", expected " << ja_compress[k] << ":" << values_compress[k] )

            BOOST_CHECK_EQUAL( rJA[k], ja_compress[k] );
            BOOST_CHECK_EQUAL( rValues[k], values_compress[k] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<CSRKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = testContext;

    getValuePos.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    LArray<IndexType> csrIA( testContext );
    LArray<IndexType> csrJA( testContext );
    LArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

    HArray<ValueType> denseValues( testContext );

    data1::getDenseTestData( numRows, numColumns, denseValues );

    ValueType zero = 0;

    {
        ReadAccess<IndexType> rIa( csrIA, loc );
        ReadAccess<IndexType> rJa( csrJA, loc );

        // comparison is done via accesses on the host

        ReadAccess<ValueType> rValues( csrValues, hostContext );
        ReadAccess<ValueType> rDense( denseValues, hostContext );

        SCAI_CONTEXT_ACCESS( loc );

        for ( IndexType i = 0; i < numRows; i++ )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                IndexType pos = getValuePos[loc]( i, j, rIa.get(), rJa.get() );

                IndexType k   = i * numColumns + j;

                if ( pos == invalidIndex )
                {
                    BOOST_CHECK_EQUAL( rDense[ k ], zero );
                }
                else
                {
                    BOOST_CHECK_EQUAL( rDense[ k], rValues[pos] );
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getValuePosColTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = testContext;
    getValuePosCol.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    //    1.0   -   2.0
    //    0.5  0.3   -
    //     -    -   3.0

    const IndexType ia[] = { 0, 2, 4, 5 };
    const IndexType ja[] = { 0, 2, 0, 1, 2 };

    const IndexType numRows = 3;
    const IndexType numValues = 5;

    HArray<IndexType> csrIA( numRows + 1, ia, testContext );
    HArray<IndexType> csrJA( numValues, ja, testContext );

    HArray<IndexType> row;   // result for rowIndexes
    HArray<IndexType> pos;   // result for positions

    IndexType cnt;

    IndexType columnIndex = 1;   // has 1 entry

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rCSRIA( csrIA, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        WriteOnlyAccess<IndexType> wRow( row, loc, numRows );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numRows );
        cnt = getValuePosCol[loc]( wRow.get(), wPos.get(), columnIndex, rCSRIA.get(), numRows, rCSRJA.get(), numValues );
    }

    BOOST_REQUIRE_EQUAL( cnt, IndexType( 1 ) );   //  only one entry for column 1

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rRow( row );

        BOOST_CHECK_EQUAL( IndexType( 1 ), rRow[0] );   // is in entry row
        BOOST_CHECK_EQUAL( IndexType( 3 ), rPos[0] );   // value of for (1,1) is at pos 3
    }

    columnIndex = 2;
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rCSRIA( csrIA, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        WriteOnlyAccess<IndexType> wRow( row, loc, numRows );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numRows );
        cnt = getValuePosCol[loc]( wRow.get(), wPos.get(), columnIndex, rCSRIA.get(), numRows, rCSRJA.get(), numValues );
    }

    BOOST_REQUIRE_EQUAL( cnt, IndexType( 2 ) );   //  two entries for column 2, order might be arbitrary

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rRow( row );
        ReadAccess<IndexType> rJA( csrJA );
        ReadAccess<IndexType> rIA( csrIA );

        for ( IndexType k = 0; k < cnt; ++k )
        {
            IndexType p = rPos[k];
            IndexType i = rRow[k];
            BOOST_CHECK_EQUAL( rJA[ p ], columnIndex );
            BOOST_CHECK( rIA[i] <= p );
            BOOST_CHECK( p < rIA[i + 1] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeNonSquareTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;

    ContextPtr loc = testContext;
    convertCSR2CSC.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> non-square test for " << *testContext << " on " << *loc )

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

    {
        ReadAccess<IndexType> rCSRIA( csrIA, loc );
        ReadAccess<IndexType> rCSRJA( csrJA, loc );
        ReadAccess<ValueType> rCSRValues( csrValues, loc );
        WriteOnlyAccess<IndexType> wCSCIA( cscIA, loc, numColumns + 1 );
        WriteOnlyAccess<IndexType> wCSCJA( cscJA, loc, numValues );
        WriteOnlyAccess<ValueType> wCSCValues( cscValues, loc, numValues );
        SCAI_CONTEXT_ACCESS( loc );
        convertCSR2CSC[loc]( wCSCIA.get(), wCSCJA.get(), wCSCValues.get(),
                             rCSRIA.get(), rCSRJA.get(), rCSRValues.get(), numRows,
                             numColumns, numValues );
    }

    //  For comparison later we sort cscJA and cscValue

    LAMAKernel<CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;

    loc = testContext;
    sortRowElements.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    SCAI_LOG_INFO( logger, "sortRowElements< " << TypeTraits<ValueType>::id() << "> for " << *testContext << " on " << *loc )

    {
        ReadAccess<IndexType> rCSCIA( cscIA, loc );
        WriteAccess<IndexType> wCSCJA( cscJA, loc );
        WriteAccess<ValueType> wCSCValues( cscValues, loc );

        SCAI_CONTEXT_ACCESS( loc );

        bool diagonalFlag = false;

        // For comparison of cscJA and cscValue we need to sort it

        sortRowElements[loc]( wCSCJA.get(), wCSCValues.get(), rCSCIA.get(), numColumns, diagonalFlag );
    }

    //  compare with the CSC test data

    IndexType numRows1;
    IndexType numColumns1;
    IndexType numValues1;

    LArray<IndexType> expIA( testContext );
    LArray<IndexType> expJA( testContext );
    LArray<ValueType> expValues( testContext );

    data1::getCSCTestData( numRows1, numColumns1, numValues1, expIA, expJA, expValues );

    // test cscIA = expJA, cscJA = expIA, cscValues = expValues

    BOOST_REQUIRE_EQUAL( expJA.size(), cscIA.size() );
    BOOST_CHECK_EQUAL( IndexType( 0 ), expJA.maxDiffNorm( cscIA ) );

    BOOST_REQUIRE_EQUAL( expIA.size(), cscJA.size() );
    BOOST_CHECK_EQUAL( IndexType( 0 ), expIA.maxDiffNorm( cscJA ) );

    BOOST_REQUIRE_EQUAL( expValues.size(), cscValues.size() );
    BOOST_CHECK_EQUAL( ValueType( 0 ), expValues.maxDiffNorm( cscValues ) );
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

    LAMAKernel<CSRKernelTrait::matrixMultiplySizes> matrixMultiplySizes;

    ContextPtr loc = testContext;
    matrixMultiplySizes.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    SCAI_LOG_INFO( logger, "matmul< " << TypeTraits<ValueType>::id() << "> non-square test" )

    // REMARK: test with explicit zeros because cusparse often has problem with it
    //       array 1             array 2
    //    1.0   -   2.0       1.0  0.5  0.0  4.0
    //    0.5  0.3   -         -   0.3  0.0  1.5
    //     -    -   3.0       2.0   -   3.0   -
    //    4.0  1.5   -

    const IndexType ia1[] = { 0, 2, 4, 5, 7 };
    const IndexType ja1[] = { 0, 2, 0, 1, 2, 0, 1 };
    const ValueType values1[] = { 1.0, 2.0, 0.5, 0.3, 3.0, 4.0, 1.5 };

    const IndexType ia2[] = { 0, 4, 7, 9 };
    const IndexType ja2[] = { 0, 1, 2, 3, 1, 2, 3, 0, 2 };
    const ValueType values2[] = { 1.0, 0.5, 0.0, 4.0, 0.3, 0.0, 1.5, 2.0, 3.0 };

    // REMARK: explicit zeros are also in result, when no compress is called
    //       array3 ( 4 x 4 )
    //
    //     5.0  0.5   6.0  4.0
    //     0.5  0.34  0.0  2.45
    //     6.0   -    9.0   -
    //     4.0  2.45  0.0  18.25

    const IndexType ia3[] = { 0, 4, 8, 10, 14 };
    const IndexType ja3[] = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 2, 0, 1, 2, 3 };
    const ValueType values3[] = { 5.0, 0.5, 6.0, 4.0, 0.5, 0.34, 0.0, 2.45, 6.0, 9.0, 4.0, 2.45, 0.0, 18.25 };

    const IndexType n1 = 4;
    const IndexType n2 = 3;
    const IndexType n3 = 4;

    const IndexType nnz1 = sizeof( ja1 ) / sizeof( IndexType );
    const IndexType nnz2 = sizeof( ja2 ) / sizeof( IndexType );

    // csr arrays for matrix a

    HArray<IndexType> aIA( n1 + 1, ia1, testContext );
    HArray<IndexType> aJA( nnz1, ja1, testContext );
    HArray<ValueType> aValues( nnz1, values1, testContext );

    // csr arrays for matrix b

    HArray<IndexType> bIA( n2 + 1, ia2, testContext );
    HArray<IndexType> bJA( nnz2, ja2, testContext );
    HArray<ValueType> bValues( nnz2, values2, testContext );

    bool diagonalProperty = false;

    SCAI_LOG_INFO( logger, "matrixMultiplySizes< " << TypeTraits<ValueType>::id() << "> non-square test for " << *testContext << " on " << *loc )

    IndexType nnz3;
    HArray<IndexType> cSizes;

    {
        ReadAccess<IndexType> rAIA( aIA, loc );
        ReadAccess<IndexType> rAJA( aJA, loc );
        ReadAccess<IndexType> rBIA( bIA, loc );
        ReadAccess<IndexType> rBJA( bJA, loc );

        SCAI_CONTEXT_ACCESS( loc );

        WriteOnlyAccess<IndexType> wSizes( cSizes, loc, n1 + 1 );

        nnz3 = matrixMultiplySizes[loc]( wSizes.get(), n1, n3, n2, diagonalProperty,
                                         rAIA.get(), rAJA.get(), rBIA.get(), rBJA.get() );
    }

    BOOST_CHECK_EQUAL( ia3[ n1 ], nnz3 );

    // check sizes
    {
        ContextPtr host = Context::getHostPtr();

        ReadAccess<IndexType> rSizes( cSizes, host );

        for ( IndexType j = 0; j <= n1; ++j )
        {
            BOOST_CHECK_EQUAL( rSizes[j], ia3[j] );
        }
    }

    SCAI_LOG_INFO( logger, "matrixMultiply< " << TypeTraits<ValueType>::id() << "> non-square test for " << *testContext << " on " << *loc )

    LAMAKernel<CSRKernelTrait::matrixMultiply<ValueType> > matrixMultiply;

    loc = testContext;
    matrixMultiply.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    ValueType alpha = 1.0;

    HArray<IndexType> cJa;
    HArray<ValueType> cValues;
    {
        ReadAccess<IndexType> rAIA( aIA, loc );
        ReadAccess<IndexType> rAJA( aJA, loc );
        ReadAccess<ValueType> rAValues( aValues, loc );

        ReadAccess<IndexType> rBIA( bIA, loc );
        ReadAccess<IndexType> rBJA( bJA, loc );
        ReadAccess<ValueType> rBValues( bValues, loc );

        ReadAccess<IndexType> rIa( cSizes, loc );

        SCAI_CONTEXT_ACCESS( loc );

        WriteOnlyAccess<IndexType> wJa( cJa, loc, nnz3 );
        WriteOnlyAccess<ValueType> wValues( cValues, loc, nnz3 );

        matrixMultiply[loc]( rIa.get(), wJa.get(), wValues.get(), n1, n3, n2, alpha, diagonalProperty,
                             rAIA.get(), rAJA.get(), rAValues.get(), rBIA.get(), rBJA.get(), rBValues.get() );
    }

    // sort the columns, otherwise comparison might fail

    {
        ReadAccess<IndexType> rCIA( cSizes );
        WriteAccess<IndexType> wCJA( cJa );
        WriteAccess<ValueType> wCValues( cValues );

        OpenMPCSRUtils::sortRowElements( wCJA.get(), wCValues.get(), rCIA.get(), n1, false );
    }

    // check ja
    {
        ContextPtr host = Context::getHostPtr();

        ReadAccess<IndexType> rJA( cJa, host );

        for ( IndexType j = 0; j < nnz3; ++j )
        {
            BOOST_CHECK_EQUAL( rJA[j], ja3[j] );
        }
    }

    // check values
    {
        ContextPtr host = Context::getHostPtr();

        ReadAccess<ValueType> rValues( cValues, host );

        for ( IndexType j = 0; j < nnz3; ++j )
        {
            ValueType x = rValues[j] - values3[j];
            BOOST_CHECK_SMALL( common::Math::real( x ), common::TypeTraits<ValueType>::small() );
            BOOST_CHECK_SMALL( common::Math::imag( x ), common::TypeTraits<ValueType>::small() );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( matAddTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<CSRKernelTrait::matrixAddSizes> matrixAddSizes;
    LAMAKernel<CSRKernelTrait::matrixAdd<ValueType> > matrixAdd;

    ContextPtr loc = testContext;
    matrixAdd.getSupportedContext( loc, matrixAddSizes );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    SCAI_LOG_INFO( logger, "matmul< " << TypeTraits<ValueType>::id() << "> non-square test for " << *testContext << " on " << *loc )

    //       array 1                array 2
    //    1.0   -   2.0  -         1.0  0.5   -   -
    //    0.5  0.3   -   0.5        -    -    -   0.5
    //     -    -   3.0  -        2.0   -   1.0   0.5

    const IndexType ia1[]     = { 0,        2,             5, 6 };
    const IndexType ja1[]     = { 0,   2,   0,   1,   3,   2 };
    const ValueType values1[] = { 1.0, 2.0, 0.5, 0.3, 0.5, 3.0 };

    const IndexType ia2[]     = { 0,        2,   3,         6 };
    const IndexType ja2[]     = { 0,   1,   3,   0,   2,  3 };
    const ValueType values2[] = { 1.0, 0.5, 0.5, 2.0, 1.0, 0.5 };

    //       array3 = array 1 + array 2
    //
    //     2.0  0.5  2.0  -
    //     0.5  0.3   -   1.0
    //     2.0   -   4.0  0.5

    const IndexType ia3[]     = { 0,             3,             6,            9 };
    const IndexType ja3[]     = { 0,   1,   2,   0,   1,   3,   0,   2,   3 };
    const ValueType values3[] = { 2.0, 0.5, 2.0, 0.5, 0.3, 1.0, 2.0, 4.0, 0.5 };

    const IndexType n1 = 3;
    const IndexType n2 = 4;

    const IndexType nnz1 = sizeof( ja1 ) / sizeof( IndexType );
    const IndexType nnz2 = sizeof( ja2 ) / sizeof( IndexType );

    // csr arrays for matrix a

    HArray<IndexType> aIA( n1 + 1, ia1, testContext );
    HArray<IndexType> aJA( nnz1, ja1, testContext );
    HArray<ValueType> aValues( nnz1, values1, testContext );

    // csr arrays for matrix b

    HArray<IndexType> bIA( n2 + 1, ia2, testContext );
    HArray<IndexType> bJA( nnz2, ja2, testContext );
    HArray<ValueType> bValues( nnz2, values2, testContext );

    bool diagonalProperty = false;

    IndexType nnz3;
    HArray<IndexType> cSizes;

    SCAI_LOG_INFO( logger, "matrixAddSizes @ " << *loc )

    {
        ReadAccess<IndexType> rAIA( aIA, loc );
        ReadAccess<IndexType> rAJA( aJA, loc );
        ReadAccess<IndexType> rBIA( bIA, loc );
        ReadAccess<IndexType> rBJA( bJA, loc );

        SCAI_CONTEXT_ACCESS( loc );

        WriteOnlyAccess<IndexType> wSizes( cSizes, loc, n1 + 1 );

        nnz3 = matrixAddSizes[loc]( wSizes.get(), n1, n2, diagonalProperty,
                                    rAIA.get(), rAJA.get(), rBIA.get(), rBJA.get() );

    }

    BOOST_CHECK_EQUAL( ia3[ n1 ], nnz3 );

    // check sizes
    {
        ContextPtr host = Context::getHostPtr();

        ReadAccess<IndexType> rSizes( cSizes, host );

        for ( IndexType j = 0; j <= n1; ++j )
        {
            SCAI_LOG_TRACE( logger, "addSizes, size[" << j << "] = " << rSizes[j] )
            BOOST_CHECK_EQUAL( rSizes[j], ia3[j] );
        }
    }

    HArray<IndexType> cJA;
    HArray<ValueType> cValues;

    {
        ReadAccess<IndexType> rAIA( aIA, loc );
        ReadAccess<IndexType> rAJA( aJA, loc );
        ReadAccess<ValueType> rAValues( aValues, loc );

        ReadAccess<IndexType> rBIA( bIA, loc );
        ReadAccess<IndexType> rBJA( bJA, loc );
        ReadAccess<ValueType> rBValues( bValues, loc );

        ReadAccess<IndexType> rCIA( cSizes, loc );

        WriteOnlyAccess<IndexType> wCJA( cJA, loc, nnz3 );
        WriteOnlyAccess<ValueType> wCValues( cValues, loc, nnz3 );

        ValueType one = 1;

        SCAI_CONTEXT_ACCESS( loc );

        matrixAdd[loc]( wCJA.get(), wCValues.get(),
                        rCIA.get(), n1, n2, diagonalProperty,
                        one, rAIA.get(), rAJA.get(), rAValues.get(),
                        one, rBIA.get(), rBJA.get(), rBValues.get() );

    }

    // sort the columns, otherwise comparison might fail

    {
        ReadAccess<IndexType> rCIA( cSizes );
        WriteAccess<IndexType> wCJA( cJA );
        WriteAccess<ValueType> wCValues( cValues );

        for ( IndexType k = 0; k < nnz3; ++k )
        {
            SCAI_LOG_TRACE( logger, "matrixAdd, ja[" << k << "] = " << wCJA[k] )
            SCAI_LOG_TRACE( logger, "matrixAdd, values[" << k << "] = " << wCValues[k] )
        }

        OpenMPCSRUtils::sortRowElements( wCJA.get(), wCValues.get(), rCIA.get(), n1, false );
    }

    {
        ReadAccess<IndexType> rCJA( cJA );
        ReadAccess<ValueType> rCValues( cValues );

        for ( IndexType k = 0; k < nnz3; ++k )
        {
            BOOST_CHECK_EQUAL( rCJA[k], ja3[k] );
            BOOST_CHECK_EQUAL( rCValues[k], values3[k] );
        }
    }
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
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<CSRKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = testContext;

    jacobi.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext << " on " << *loc )

    HArray<IndexType> csrIA( testContext );
    HArray<IndexType> csrJA( testContext );
    HArray<ValueType> csrValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data2::getCSRTestData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );

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

        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( csrIA, loc );
            ReadAccess<IndexType> rJA( csrJA, loc );
            ReadAccess<ValueType> rValues( csrValues, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<ValueType> rRhs( rhs, loc );
            WriteOnlyAccess<ValueType> wSolution( res, loc, numColumns );

            jacobi[loc]( wSolution.get(),
                         rIA.get(), rJA.get(), rValues.get(),
                         rOld.get(), rRhs.get(), omega, numRows );
        }

        LArray<ValueType> expectedRes( testContext );

        data2::getJacobiResult( expectedRes, oldSolution, omega, rhs );

        ValueType maxDiff = expectedRes.maxDiffNorm( res );

        BOOST_CHECK( common::Math::real( maxDiff ) < 0.1 );

        bool mustBeIdentical = false;

        if ( mustBeIdentical )
        {
            ReadAccess<ValueType> rExpected( expectedRes );
            ReadAccess<ValueType> rComputed( res );

            for ( IndexType i = 0; i < numRows; ++i )
            {
                BOOST_CHECK_EQUAL( rExpected[i], rComputed[i] );
            }
        }
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

            jacobiHalo[loc]( wSolution.get(),
                             rIADummy.get(), rDiag.get(),
                             rIA.get(), rJA.get(), rValues.get(), rIndexes.get(),
                             rOld.get(), omega, numNonEmptyRows );
        }

        LArray<ValueType> expectedSol( numRows, ValueType( 0 ), testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        ValueType maxDiff = expectedSol.maxDiffNorm( solution );

        BOOST_CHECK( common::Math::real( maxDiff ) < 0.1 );

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

    static LAMAKernel<CSRKernelTrait::jacobiHaloWithDiag<ValueType> > jacobiHaloWithDiag;

    ContextPtr loc = testContext;

    jacobiHaloWithDiag.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    SCAI_LOG_INFO( logger, "jacobiHaloWithDiag test for " << *testContext << " on " << *loc )

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

            jacobiHaloWithDiag[loc]( wSolution.get(),
                                     rDiag.get(),
                                     rIA.get(), rJA.get(), rValues.get(), rIndexes.get(),
                                     rOld.get(), omega, numNonEmptyRows );
        }

        LArray<ValueType> expectedSol( numRows, ValueType( 0 ), testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        ValueType maxDiff = expectedSol.maxDiffNorm( solution );

        BOOST_CHECK( common::Math::real( maxDiff ) < 0.1 );

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
