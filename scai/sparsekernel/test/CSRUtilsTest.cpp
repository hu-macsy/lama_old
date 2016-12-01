/**
 * @file CSRUtilsTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/test/TestMacros.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;
using namespace hmemo;
using namespace sparsekernel;
using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CSRUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxDiffValTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    kregistry::KernelTraitContextFunction<CSRKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;
    ContextPtr loc = Context::getContextPtr( absMaxDiffVal.validContext( testContext->getType() ) );
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
    ValueType maxVal = absMaxDiffVal[loc->getType()]( numRows, false, rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), rCSRIA2.get(),
                       rCSRJA2.get(), rCSRValues2.get() );
    BOOST_CHECK_EQUAL( 3, maxVal );
    // rows are sorted, so we can also apply sortFlag = true
    maxVal = absMaxDiffVal[loc->getType()]( numRows, true, rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), rCSRIA2.get(),
                                            rCSRJA2.get(), rCSRValues2.get() );
    BOOST_CHECK_EQUAL( 3, maxVal );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeSquareTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    kregistry::KernelTraitContextFunction<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;
    ContextPtr loc = Context::getContextPtr( convertCSR2CSC.validContext( testContext->getType() ) );
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
        convertCSR2CSC[loc->getType()]( wCSCIA.get(), wCSCJA.get(), wCSCValues.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(), numRows,
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

    ContextPtr testContext = Context::getContextPtr();

    kregistry::KernelTraitContextFunction<CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;

    ContextPtr loc = Context::getContextPtr( sortRowElements.validContext( testContext->getType() ) );

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

        sortRowElements[loc->getType()]( wJA.get(), wValues.get(), rIA.get(), numRows, diagonalFlag );
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

        sortRowElements[loc->getType()]( wJA.get(), wValues.get(), rIA.get(), numRows, diagonalFlag );
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

    ContextPtr testContext = Context::getContextPtr();

    kregistry::KernelTraitContextFunction<CSRKernelTrait::countNonZeros<ValueType> > countNonZeros;

    ContextPtr loc = Context::getContextPtr( countNonZeros.validContext( testContext->getType() ) );

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

        countNonZeros[loc->getType()]( wSizes.get(), rIA.get(), rJA.get(), rValues.get(),
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

    ContextPtr testContext = Context::getContextPtr();

    kregistry::KernelTraitContextFunction<CSRKernelTrait::compress<ValueType> > compress;

    ContextPtr loc = Context::getContextPtr( compress.validContext( testContext->getType() ) );

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

        compress[loc->getType()]( wJA.get(), wValues.get(), rIA.get(), roIA.get(), roJA.get(), roValues.get(),
                                  numRows, eps, diagonalFlag );
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

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getValuePosColTest )
{
    ContextPtr testContext = Context::getContextPtr();

    kregistry::KernelTraitContextFunction<CSRKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = Context::getContextPtr( getValuePosCol.validContext( testContext->getType() ) );

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
        cnt = getValuePosCol[loc->getType()]( wRow.get(), wPos.get(), columnIndex, rCSRIA.get(), numRows, rCSRJA.get(), numValues );
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
        cnt = getValuePosCol[loc->getType()]( wRow.get(), wPos.get(), columnIndex, rCSRIA.get(), numRows, rCSRJA.get(), numValues );
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
            BOOST_CHECK( p < rIA[i+1] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeNonSquareTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    kregistry::KernelTraitContextFunction<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;
    ContextPtr loc = Context::getContextPtr( convertCSR2CSC.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected
    SCAI_LOG_INFO( logger, "transpose< " << TypeTraits<ValueType>::id() << "> non-square test for " << *testContext << " on " << *loc )

    //      input array           transpose
    //  0     1.0   -   2.0       1.0  0.5   -    4.0
    //  1     0.5  0.3   -         -   0.3   -    1.5
    //  2      -    -   3.0       2.0   -   3.0    -
    //  3     4.0  1.5   -

    const IndexType ia1[] = { 0, 2, 4, 5, 7 };
    const IndexType ja1[] = { 0, 2, 0, 1, 2, 0, 1 };

    const IndexType ia2[] = { 0, 3, 5, 7 };
    const IndexType ja2[] = { 0, 1, 3, 1, 3, 0, 2 };

    const ValueType values1[] =  { 1.0, 2.0, 0.5, 0.3, 3.0, 4.0, 1.5 };
    const ValueType values2[] =  { 1.0, 0.5, 4.0, 0.3, 1.5, 2.0, 3.0 };

    const IndexType numRows = 4;
    const IndexType numColumns = 3;
    const IndexType numValues = 7;

    const IndexType ia1_size = sizeof( ia1 ) / sizeof( IndexType );
    const IndexType ia2_size = sizeof( ia2 ) / sizeof( IndexType );

    BOOST_CHECK_EQUAL( numRows + 1, ia1_size );
    BOOST_CHECK_EQUAL( numColumns + 1, ia2_size );

    const size_t sizeArray = numValues;

    BOOST_CHECK_EQUAL( sizeArray, sizeof( ja1 ) / sizeof( IndexType ) );
    BOOST_CHECK_EQUAL( sizeArray, sizeof( ja2 ) / sizeof( IndexType ) );

    BOOST_CHECK_EQUAL( sizeArray, sizeof( values1 ) / sizeof( ValueType ) );
    BOOST_CHECK_EQUAL( sizeArray, sizeof( values2 ) / sizeof( ValueType ) );

    HArray<IndexType> csrIA( numRows + 1, ia1, testContext );
    HArray<IndexType> csrJA( numValues, ja1, testContext );
    HArray<ValueType> csrValues( numValues, values1, testContext );

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
        convertCSR2CSC[loc->getType()]( wCSCIA.get(), wCSCJA.get(), wCSCValues.get(), 
                                        rCSRIA.get(), rCSRJA.get(), rCSRValues.get(), numRows,
                                        numColumns, numValues );
    }

    BOOST_REQUIRE_EQUAL( numColumns + 1, cscIA.size() );
    {
        ContextPtr host = Context::getHostPtr();

        ReadAccess<IndexType> rCSCIA( cscIA, host );

        for ( IndexType j = 0; j <= numColumns; ++j )
        {
            BOOST_REQUIRE_EQUAL( rCSCIA[j], ia2[j] );
        }
    }

    BOOST_REQUIRE_EQUAL( numValues, cscJA.size() );
    BOOST_REQUIRE_EQUAL( numValues, cscValues.size() );

    //  For comparison later we sort cscJA and cscValue

    kregistry::KernelTraitContextFunction<CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;

    loc = Context::getContextPtr( sortRowElements.validContext( testContext->getType() ) );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    SCAI_LOG_INFO( logger, "sortRowElements< " << TypeTraits<ValueType>::id() << "> for " << *testContext << " on " << *loc )

    {
        ReadAccess<IndexType> rCSCIA( cscIA, loc );
        WriteAccess<IndexType> wCSCJA( cscJA, loc );
        WriteAccess<ValueType> wCSCValues( cscValues, loc );

        SCAI_CONTEXT_ACCESS( loc );

        bool diagonalFlag = false;

        // For comparison of cscJA and cscValue we need to sort it

        sortRowElements[loc->getType()]( wCSCJA.get(), wCSCValues.get(), rCSCIA.get(), numColumns, diagonalFlag );
    }

    // check CSC for correctness, done on host
    {
        ContextPtr host = Context::getHostPtr();

        ReadAccess<IndexType> rCSCJA( cscJA, host );
        ReadAccess<ValueType> rCSCValues( cscValues, host );

        for ( IndexType j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( rCSCJA[j], ja2[j] );
        }

        for ( IndexType j = 0; j < numValues; ++j )
        {
            BOOST_CHECK_EQUAL( rCSCValues[j], values2[j] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

typedef boost::mpl::list<SCAI_NUMERIC_TYPES_EXT_HOST> scai_ext_test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( decompositionTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();
    kregistry::KernelTraitContextFunction<CSRKernelTrait::decomposition<ValueType> > decomposition;
    ContextPtr loc = Context::getContextPtr( decomposition.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected
    SCAI_LOG_INFO( logger, "decomposition< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    if ( common::TypeTraits<IndexType>::stype != common::scalar::INT )
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
                                       2.0, -3.0, 1.0 };
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
        decomposition[loc->getType()]( wSol.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(),
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
    ContextPtr testContext = Context::getContextPtr();

    kregistry::KernelTraitContextFunction<CSRKernelTrait::matrixMultiplySizes> matrixMultiplySizes;

    ContextPtr loc = Context::getContextPtr( matrixMultiplySizes.validContext( testContext->getType() ) );

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

        nnz3 = matrixMultiplySizes[loc->getType()]( wSizes.get(), n1, n3, n2, diagonalProperty, 
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

    kregistry::KernelTraitContextFunction<CSRKernelTrait::matrixMultiply<ValueType> > matrixMultiply;

    loc = Context::getContextPtr( matrixMultiply.validContext( testContext->getType() ) );

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

        matrixMultiply[loc->getType()]( rIa.get(), wJa.get(), wValues.get(), n1, n3, n2, alpha, diagonalProperty,
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
    ContextPtr testContext = Context::getContextPtr();

    kregistry::KernelTraitContextFunction<CSRKernelTrait::matrixAddSizes> matrixAddSizes;
    kregistry::KernelTraitContextFunction<CSRKernelTrait::matrixAdd<ValueType> > matrixAdd;

    ContextPtr loc = Context::getContextPtr( matrixAdd.validContext( testContext->getType() ) );

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

        nnz3 = matrixAddSizes[loc->getType()]( wSizes.get(), n1, n2, diagonalProperty, 
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

        matrixAdd[loc->getType()]( wCJA.get(), wCValues.get(), 
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

BOOST_AUTO_TEST_SUITE_END()
