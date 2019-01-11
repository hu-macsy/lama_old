/**
 * @file ELLUtilsTest.cpp
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
 * @brief Contains tests for kernel implementations of ELLKernelTrait routines.
 * @author Thomas Brandes
 * @date 15.10.2012
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/sparsekernel/ELLKernelTrait.hpp>
#include <scai/sparsekernel/ELLUtils.hpp>
#include <scai/sparsekernel/openmp/OpenMPELLUtils.hpp>
#include <scai/utilskernel.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/kregistry.hpp>
#include <scai/hmemo.hpp>
#include <scai/common/Math.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>
#include <scai/sparsekernel/test/TestData1.hpp>
#include <scai/sparsekernel/test/TestData2.hpp>

#include <scai/hmemo/test/ContextFix.hpp>

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

BOOST_AUTO_TEST_SUITE( ELLUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.ELLUtilsTest" )

/* ------------------------------------------------------------------------------------- */

// BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxDiffValTest, ValueType, scai_numeric_test_types )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( diagonalPositionsTest )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data2::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    // all available test

    HArray<IndexType> diagonalPositions;

    const IndexType numDiagonals = common::Math::min( numRows, numColumns );

    ELLUtils::getDiagonalPositions( diagonalPositions, numRows, numColumns, ellIA, ellJA, testContext );

    BOOST_CHECK_EQUAL( diagonalPositions.size(), numDiagonals );

    HArray<IndexType> expDiagonalPositions( { 0, 5, 6, 11 } );

    BOOST_TEST( hostReadAccess( diagonalPositions ) == hostReadAccess( expDiagonalPositions ), per_element() );

    // negative test

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    // test empty array

    HArrayUtils::setSameValue<IndexType>( ellIA, numRows, 0 );

    ellJA.clear();

    IndexType zeroDiags = ELLUtils::getDiagonalPositions( diagonalPositions, numRows, numColumns, ellIA, ellJA, testContext );

    BOOST_CHECK_EQUAL( zeroDiags, 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( diagonalTest )
{   
    typedef DefaultReal ValueType;
    
    ContextPtr testContext = ContextFix::testContext;

    //   Input storage:
    //   -------------
    //    1.5   -   2.0  1.1  - 
    //    0.5  0.0   -    -   - 
    //     -    -   3.0   -   - 
    //     -   0.0  4.0  0.0 1.0
    
    const IndexType x = 0;
    const ValueType z = 0;

    // ELL data contains entries for all diagonal elements even if they are zero
    // For accessing the diagonal there is no need for sorted entries per row

    HArray<IndexType> ia(     {   3,             2,        1,   4  }, testContext );
    HArray<IndexType> ja(     {   2,   0,   2,   1,   0,   1, x,   2,   3, x, x,   3, x, x, x, 4  },  testContext );
    HArray<ValueType> values( { 2.0, 0.5, 3.0, 0.0, 1.5, 0.0, z, 4.0, 1.1, z, z, 0.0, z, z, z, 1.0 },  testContext );
    
    HArray<ValueType> expDiag( { 1.5, 0.0, 3.0, 0.0 } );

    const IndexType m = 4;
    const IndexType n = 5;
    
    HArray<ValueType> diag;

    ELLUtils::getDiagonal( diag, m, n, ia, ja, values, testContext );

    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( expDiag ), per_element() );
    
    HArray<ValueType> newDiag( { 1.2, 2.0, 3.3, 0.5 } ); 

    ELLUtils::setDiagonalV( values, newDiag, m, n, ia, ja, testContext );

    ELLUtils::getDiagonal( diag, m, n, ia, ja, values, testContext );

    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( newDiag ), per_element() );

    ValueType diagVal = 1;

    ELLUtils::setDiagonal( values, diagVal, m, n, ia, ja, testContext );
    ELLUtils::getDiagonal( diag, m, n, ia, ja, values, testContext );

    BOOST_TEST( hostReadAccess( diag ) == hostReadAccess( HArray<ValueType>( m, diagVal) ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( fillELlValuesTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<ELLKernelTrait::fillELLValues<ValueType> > fillELLValues;

    ContextPtr loc = testContext;
    fillELLValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    {
        WriteAccess<IndexType> wJA( ellJA, loc );
        WriteAccess<ValueType> wValues( ellValues, loc );
        ReadAccess<IndexType> rIA( ellIA, loc );

        SCAI_CONTEXT_ACCESS( loc );

        fillELLValues[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, numValuesPerRow );
    }

    HArray<ValueType> xDummy( numColumns, ValueType( 1 ) );

    {
        ReadAccess<IndexType> rIA( ellIA, hostContext );
        ReadAccess<IndexType> rJA( ellJA, hostContext );
        ReadAccess<ValueType> rValues( ellValues, hostContext );
        ReadAccess<ValueType> rX( xDummy, hostContext );

        ValueType testVal = 0;

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = rIA[i]; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = jj * numRows + i;
                IndexType j   = rJA[ pos ];
                SCAI_ASSERT_VALID_INDEX( j, numColumns, "illegal col pos" )
                testVal += rValues[ pos ] * rX[ j ];
            }
        }

        BOOST_CHECK_EQUAL( ValueType( 0 ), testVal );
    }
}


/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( checkTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::check> check;

    ContextPtr loc = testContext;
    check.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    // check with correct values
    {
        const IndexType valuesIa[] = { 4, 3, 5, 2 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof( IndexType );
        const IndexType minusOne = static_cast<IndexType>( -1 );
        const IndexType valuesJa[] = { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, minusOne, 4, minusOne, 4, 0, 0, 0, 5, 0 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = 5;
        const IndexType numColumns = 6;
        HArray<IndexType> ia( nIa, valuesIa, testContext );
        HArray<IndexType> ja( nJa, valuesJa, testContext );
        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );
        SCAI_CONTEXT_ACCESS( loc );
        BOOST_CHECK_NO_THROW( check[loc]( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ) );
    }
    // check with invalid ia
    {
        HArray<IndexType> ia( { 4, 3, 7, 2 }, testContext );   // contains at least one illegal size
        HArray<IndexType> ja( { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 0, 4, 0, 4, 0, 0, 0, 5, 0 }, testContext );
        
        const IndexType numRows = ia.size();
        const IndexType numValuesPerRow = ja.size() / numRows;
        SCAI_ASSERT_EQ_ERROR( numRows * numValuesPerRow, ja.size(), "illegal size of ja" )
        const IndexType numColumns = 5;

        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );
        SCAI_CONTEXT_ACCESS( loc );
        BOOST_CHECK_THROW( check[loc]( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ),
                           Exception );
    }
    // check with invalid ja
    {
        const IndexType valuesIa[] = { 4, 3, 5, 2 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof( IndexType );
        const IndexType valuesJa[] = { 1, 1, 1, 1, 2, 2, 2, static_cast<IndexType>( -1 ), 3, 3, 3, 0, 4, 0, 4, 0, 0, 0, 5, 0 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = 5;
        const IndexType numColumns = 5;
        HArray<IndexType> ia( nIa, valuesIa, testContext );
        HArray<IndexType> ja( nJa, valuesJa, testContext );
        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );
        SCAI_CONTEXT_ACCESS( loc );
        BOOST_CHECK_THROW( check[loc]( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ),
                           Exception );
    }
    // check with valid empty values
    {
        const IndexType numRows = 0;
        const IndexType numValuesPerRow = 0;
        const IndexType numColumns = 0;
        HArray<IndexType> ia;
        HArray<IndexType> ja;
        ReadAccess<IndexType> rIa( ia, loc );
        ReadAccess<IndexType> rJa( ja, loc );
        SCAI_CONTEXT_ACCESS( loc );
        BOOST_CHECK_NO_THROW( check[loc]( numRows, numValuesPerRow, numColumns, rIa.get(), rJa.get(), "checkTest" ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getRowTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::getRow<ValueType> > getRow;

    ContextPtr loc = testContext;
    getRow.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    // check with valid dense values
    {
        ValueType valuesValues[] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
        IndexType valuesIa[] = { 5, 5, 5 };
        const IndexType nIa = sizeof( valuesIa ) / sizeof( IndexType );
        IndexType valuesJa[] = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        const IndexType i = 1;
        const IndexType numRows = nIa;
        const IndexType numValuesPerRow = nJa / nIa;
        const IndexType numColumns = 5;
        HArray<ValueType> values( nValues, valuesValues, testContext );
        HArray<IndexType> ia( nIa, valuesIa, testContext );
        HArray<IndexType> ja( nJa, valuesJa, testContext );
        HArray<ValueType> row( numColumns, ValueType( 0 ) );
        {
            ReadAccess<ValueType> rValues( values, loc );
            ReadAccess<IndexType> rIa( ia, loc );
            ReadAccess<IndexType> rJa( ja, loc );
            WriteOnlyAccess<ValueType> wRow( row, loc, numColumns );
            SCAI_CONTEXT_ACCESS( loc );
            getRow[loc]( wRow.get(), i, numRows, numColumns, numValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
        }

        std::vector<ValueType> expectedValues( { 0, 1, 2, 3, 4 } );
        BOOST_TEST( hostReadAccess( row ) == expectedValues, per_element() );
    }
    // check with valid sparse values
    {
        HArray<IndexType> ia( { 5, 5, 5 }, testContext );
        HArray<ValueType> values( { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 }, testContext );
        HArray<IndexType> ja( { 0, 0, 0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 10, 10, 10 }, testContext );

        const IndexType numRows = ia.size();
        const IndexType numColumns = 11;
        const IndexType numValuesPerRow = ja.size() / numRows;

        BOOST_CHECK_EQUAL( numValuesPerRow * numRows, ja.size() );
        BOOST_CHECK_EQUAL( numValuesPerRow * numRows, values.size() );

        const IndexType i = 1;  // row to get

        HArray<ValueType> row( numColumns, ValueType( 0 ) );
        {
            ReadAccess<ValueType> rValues( values, loc );
            ReadAccess<IndexType> rIa( ia, loc );
            ReadAccess<IndexType> rJa( ja, loc );
            WriteOnlyAccess<ValueType> wRow( row, loc, numColumns );
            SCAI_CONTEXT_ACCESS( loc );
            getRow[loc]( wRow.get(), i, numRows, numColumns, numValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
        }

        HArray<ValueType> expectedValues( { 0, 0, 1, 0, 2, 0, 3, 0, 0, 0, 4 } );

        BOOST_TEST( hostReadAccess( row ) == hostReadAccess( expectedValues ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    IndexType ia_values[] = { 5, 5, 5 };

    HArray<ValueType> values(  { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 }, testContext );
    HArray<IndexType> ja(      { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 }, testContext );
    HArray<IndexType> ia( { 5, 5, 5 },  testContext );

    ValueType expectedValues[] =  { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };

    const IndexType numRows = ia.size();
    const IndexType numValuesPerRow = ja.size() / numRows;

    BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, ja.size() );
    BOOST_REQUIRE_EQUAL( numRows * numValuesPerRow, values.size() );

    auto rValues = hostReadAccess( values );  // for check on host

    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType j = 0; j < ia_values[i]; j++ )
        {
            IndexType pos = ELLUtils::getValuePos( i, j, ia, ja, testContext );
            BOOST_CHECK_EQUAL( expectedValues[j * numRows + i], rValues[ pos ] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    /*   Matrix:       6  0  0  4         0  3  -    6  4  -
                       7  0  0  0         0  -  -    7  -  -
                       0  0  9  4         2  3  -    9  4  -
                       2  5  0  3         0  1  3    2  5  3
                       2  0  0  1         0  3  -    2  1  -
                       0  0  0  0         -  -  -    -  -  -
                       0  1  0  2         1  3  -    1  2  -
     */

    const IndexType x = 0;
    const ValueType v = 0;

    HArray<IndexType> ellIA(     { 2, 1, 2, 3, 2, 0, 2 }, testContext );
    HArray<IndexType> ellJA(     { 0, 0, 2, 0, 0, x, 1, 3, x, 3, 1, 3, x, 3, x, x, x, 3, x, x, x }, testContext );
    HArray<ValueType> ellValues( { 6, 7, 9, 2, 2, v, 1, 4, v, 4, 5, 1, v, 2, v, v, v, 3, v, v, v }, testContext );
    
    const IndexType numRows         = 7;
    const IndexType numColumns      = 4;

    HArray<ValueType> diagonal( { -1, 2, 3, -4, -5, -6, 4 }, testContext );

    HArray<ValueType> expectedValues( { 6 * -1, 7 * 2, 9 * 3, 2 * -4, 2 * -5, v, 1 * 4, 
                                        4 * -1, v, 4 * 3, 5 * -4, 1 * -5, v, 2 * 4, 
                                        v, v, v, 3 * -4, v, v, v }, testContext );

    auto op = common::BinaryOp::MULT;

    ELLUtils::setRows( ellValues, numRows, numColumns, ellIA, ellJA, diagonal, op, testContext );

    BOOST_TEST( hostReadAccess( expectedValues ) == hostReadAccess( ellValues ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setColumnsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    /*   Matrix:       6  0  0  4         0  3  -    6  4  -
                       7  0  0  0         0  -  -    7  -  -
                       0  0  9  4         2  3  -    9  4  -
                       2  5  0  3         0  1  3    2  5  3
                       2  0  0  1         0  3  -    2  1  -
                       0  0  0  0         -  -  -    -  -  -
                       0  1  0  2         1  3  -    1  2  -
     */

    const IndexType x = 0;
    const ValueType v = 0;

    HArray<IndexType> ellIA(     { 2, 1, 2, 3, 2, 0, 2 }, testContext );
    HArray<IndexType> ellJA(     { 0, 0, 2, 0, 0, x, 1, 3, x, 3, 1, 3, x, 3, x, x, x, 3, x, x, x }, testContext );
    HArray<ValueType> ellValues( { 6, 7, 9, 2, 2, v, 1, 4, v, 4, 5, 1, v, 2, v, v, v, 3, v, v, v }, testContext );

    const IndexType numRows         = 7;
    const IndexType numColumns      = 4;

    HArray<ValueType> diagonal( { -1, 2, 3, -4 }, testContext );

    HArray<ValueType> expectedValues( { 6 * -1, 7 * -1, 9 * 3, 2 * -1, 2 * -1, v, 1 * 2, 
                                        4 * -4, v, 4 * -4, 5 * 2, 1 * -4, v, 2 * -4, 
                                        v, v, v, 3 * -4, v, v, v }, testContext );

    auto op = common::BinaryOp::MULT;

    ELLUtils::setColumns( ellValues, numRows, numColumns, ellIA, ellJA, diagonal, op, testContext );

    BOOST_TEST( hostReadAccess( expectedValues ) == hostReadAccess( ellValues ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( convertCSR2ELLTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    const IndexType numRows = 3;
    const IndexType numColumns = 5;

    const HArray<IndexType> csrIA( { 0, 5, 10, 15 }, testContext );
    const HArray<IndexType> csrJA( { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 }, testContext );
    const HArray<ValueType> csrValues( { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 }, testContext );

    HArray<ValueType> ellValues;
    HArray<IndexType> ellJA;
    HArray<IndexType> ellIA;

    ELLUtils::convertCSR2ELL( ellIA, ellJA, ellValues, 
                              numRows, numColumns, csrIA, csrJA, csrValues, testContext );

    HArray<IndexType> expEllIA( { 5, 5, 5 } );
    HArray<ValueType> expEllValues( { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 } );
    HArray<IndexType> expEllJA(     { 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 } );

    BOOST_TEST( hostReadAccess( ellIA ) == hostReadAccess( expEllIA ), per_element() );
    BOOST_TEST( hostReadAccess( ellJA ) == hostReadAccess( expEllJA ), per_element() );
    BOOST_TEST( hostReadAccess( ellValues ) == hostReadAccess( expEllValues ), per_element() );

    HArray<ValueType> newCsrValues;
    HArray<IndexType> newCsrJA;
    HArray<IndexType> newCsrIA;

    ELLUtils::convertELL2CSR( newCsrIA, newCsrJA, newCsrValues,
                              numRows, numColumns, ellIA, ellJA, ellValues, testContext );

    BOOST_TEST( hostReadAccess( csrIA ) == hostReadAccess( newCsrIA ), per_element() );
    BOOST_TEST( hostReadAccess( csrJA ) == hostReadAccess( newCsrJA ), per_element() );
    BOOST_TEST( hostReadAccess( csrValues ) == hostReadAccess( newCsrValues ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( nonEmptyRowsTest )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> ellIA{ 2, 0, 3, 0, 2, 0, 2 };
    HArray<IndexType> expRowIndexes( { 0,    2,    4,    6  } );

    HArray<IndexType> rowIndexes;  // will be set if there are less than threshold non-zero rows

    ELLUtils::nonEmptyRows( rowIndexes, ellIA, 1.0f, testContext );

    BOOST_TEST( hostReadAccess( rowIndexes ) == hostReadAccess( expRowIndexes ), per_element() );

    // set a threshold, as there are more than 20% non-empty row, no indexes are built.

    ELLUtils::nonEmptyRows( rowIndexes, ellIA, 0.2f, testContext );

    BOOST_CHECK_EQUAL( rowIndexes.size(), IndexType( 0 ) );

    ellIA = { 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
    expRowIndexes = { 0, 6 };

    // here we have less than 20% non-empty rows

    ELLUtils::nonEmptyRows( rowIndexes, ellIA, 0.2f, testContext );

    BOOST_TEST( hostReadAccess( rowIndexes ) == hostReadAccess( expRowIndexes ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( compressTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::compressValues<ValueType> > compressValues;

    ContextPtr loc = testContext;
    compressValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // Check without epsilon

    {
        /* Input Matrix:     1  4  0  0  0     0  3  4  5  6
                             2  5  0  0  8     1  3  4  5  6
                             3  6  0  7  9     2  3  4  5  6

           Output Matrix:    1  4  0  0        0  3  0  0
                             2  5  8  0        1  3  6  0
                             3  6  7  9        2  3  5  6
        */
        HArray<IndexType> ellIa( { 5, 5, 5 }, testContext );

        HArray<IndexType> expEllIa( { 2, 3, 4 } );

        HArray<ValueType> ellValues( { 1, 2, 3,
                                       4, 5, 6,
                                       0, 0, 0,
                                       0, 0, 7,
                                       0, 8, 9  }, testContext );

        HArray<ValueType> expEllValues( { 1, 2, 3,
                                          4, 5, 6,
                                          0, 8, 7,
                                          0, 0, 9  } );

        HArray<IndexType> ellJa(    { 0, 1, 2,
                                      3, 3, 3,
                                      4, 4, 4,
                                      5, 5, 5,
                                      6, 6, 6  }, testContext );

        HArray<IndexType> expEllJa( { 0, 1, 2,
                                      3, 3, 3,
                                      0, 6, 5,
                                      0, 0, 6  } );

  
        const IndexType numRows = ellIa.size();

        IndexType numValuesPerRow = ellValues.size() / numRows;

        RealType<ValueType> eps = 0.0;

        ELLUtils::compress( ellIa, ellJa, ellValues, numValuesPerRow, eps, testContext );
        
        BOOST_TEST( hostReadAccess( ellIa ) == hostReadAccess( expEllIa ), per_element() );
        BOOST_TEST( hostReadAccess( ellJa ) == hostReadAccess( expEllJa ), per_element() );
        BOOST_TEST( hostReadAccess( ellValues ) == hostReadAccess( expEllValues ), per_element() );

        eps = 1;   // entry( 0, 0 ) is now also considered as zero

        ELLUtils::compress( ellIa, ellJa, ellValues, numValuesPerRow, eps, testContext );
        
        HArray<IndexType> expEllIa1( { 1, 3, 4 } );
        HArray<ValueType> expEllValues1( { 4, 2, 3,
                                           0, 5, 6,
                                           0, 8, 7,
                                           0, 0, 9  } );
        HArray<IndexType> expEllJa1( { 3, 1, 2,
                                       0, 3, 3,
                                       0, 6, 5,
                                       0, 0, 6  } );
  
        BOOST_TEST( hostReadAccess( ellIa ) == hostReadAccess( expEllIa1 ), per_element() );
        BOOST_TEST( hostReadAccess( ellJa ) == hostReadAccess( expEllJa1 ), per_element() );
        BOOST_TEST( hostReadAccess( ellValues ) == hostReadAccess( expEllValues1 ), per_element() );
    }
}

BOOST_AUTO_TEST_CASE( matrixMultiplySizesTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::matrixMultiplySizes> matrixMultiplySizes;

    ContextPtr loc         = testContext;
    matrixMultiplySizes.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // Check with symmetric matrix
    {
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 3 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType expectedCIa[] =
        { 4, 4, 3, 4, 4 };
        IndexType numValues = 5; // all matrices have shape 5 x 5
        IndexType aNumValuesPerRow = aNumValues / numValues;
        IndexType bNumValuesPerRow = bNumValues / numValues;
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( testContext );
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            WriteOnlyAccess<IndexType> wCIa( CIa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            matrixMultiplySizes[loc]( wCIa.get(), numValues, numValues, numValues, false, rAIa.get(), rAJa.get(),
                                      aNumValuesPerRow, rBIa.get(), rBJa.get(), bNumValuesPerRow );
        }
        BOOST_CHECK_EQUAL( numValues, CIa.size() );
        ReadAccess<IndexType> rCIa( CIa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
        }
    }
    // Check with asymmetric matrix
    {
        //   A       B
        //   x x 0   x  x
        //   0 x x   x
        //   0 0 x   x
        IndexType valuesAIa[] =
        { 2, 2, 2 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 0, 1, 0, 2, 3, 3 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        IndexType valuesBIa[] =
        { 2, 1, 2, 3 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 2 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType expectedCIa[] =
        { 3, 3, 3 };
        IndexType cNumRows = sizeof( expectedCIa ) / sizeof( IndexType );
        BOOST_REQUIRE_EQUAL( aNumRows, cNumRows );
        // a and a * b have same number rows
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( testContext );
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            WriteOnlyAccess<IndexType> wCIa( CIa, loc, cNumRows );
            SCAI_CONTEXT_ACCESS( loc );
            IndexType numColumns = 5; // does not really matter
            matrixMultiplySizes[loc]( wCIa.get(), aNumRows, numColumns, bNumRows, false, rAIa.get(), rAJa.get(),
                                      aNumValuesPerRow, rBIa.get(), rBJa.get(), bNumValuesPerRow );
        }
        BOOST_CHECK_EQUAL( cNumRows, CIa.size() );
        ReadAccess<IndexType> rCIa( CIa );

        for ( IndexType i = 0; i < cNumRows; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( matrixMultiplyTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::matrixMultiply<ValueType> > matrixMultiply;

    ContextPtr loc         = testContext;
    matrixMultiply.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    // Check with symmetric matrix
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType valuesCIa[] =
        { 4, 4, 3, 4, 4 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 28, 71, 8, 84, 73, 27, 63, 27, 144, 153, 15, 99, 37, 77, 85, 8, 56, 0, 49, 112 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 2, 2, 3, 3, 0, 3, 3 };
        IndexType numValues = 20;
        ValueType alpha = 1;
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays
        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = numValues / cNumRows;
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, numValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            IndexType numColumns = 5; // not really needed here but internally used
            bool diagonalProperty = false; // do not care about it here
            matrixMultiply[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, bNumRows,
                                 diagonalProperty, alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow,
                                 rBIa.get(), rBJa.get(), rBValues.get(), bNumValuesPerRow );
        }
        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
    // Check with set alpha
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType valuesCIa[] =
        { 4, 4, 3, 4, 4 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 28, 71, 8, 84, 73, 27, 63, 27, 144, 153, 15, 99, 37, 77, 85, 8, 56, 0, 49, 112 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 2, 2, 3, 3, 0, 3, 3 };
        IndexType numValues = 20;
        ValueType alpha = 2.5;
        // input arrays, directly initialized on context device
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays
        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = numValues / cNumRows;
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, numValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            bool diagonalProperty = false; // do not care about it
            IndexType numColumns = 15; // does not matter here but internally used for optimizations
            matrixMultiply[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, bNumRows,
                                 diagonalProperty, alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow,
                                 rBIa.get(), rBJa.get(), rBValues.get(), bNumValuesPerRow );
        }
        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i]*alpha, rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
    // Check with asymmetric matrix
    {
        ValueType valuesAValues[] =
        { 2, 4, 4, 3, 1, 5 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 2, 2 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 0, 1, 0, 2, 3, 3 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        BOOST_REQUIRE_EQUAL( aNumValues, nAValues );
        ValueType valuesBValues[] =
        { 4, 3, 7, 5, 9, 0, 6, 8, 0, 0, 0, 9 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 1, 2, 3 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 2 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        BOOST_REQUIRE_EQUAL( bNumValues, nBValues );
        IndexType valuesCIa[] =
        { 3, 3, 3 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 29, 5, 41, 18, 20, 40, 18, 9, 81 };
        IndexType expectedCJa[] =
        { 0, 0, 0, 1, 1, 1, 2, 2, 2 };
        IndexType cNumValues = 9;
        ValueType alpha = 1;
        // input arrays, directly initialized on context device
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays
        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = cNumValues / cNumRows;
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, cNumValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, cNumValues );
            SCAI_CONTEXT_ACCESS( loc );
            bool diagonalProperty = false; // do not care about it
            IndexType numColumns = 15; // does not matter here but internally used for optimizations
            matrixMultiply[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, bNumRows,
                                 diagonalProperty, alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow,
                                 rBIa.get(), rBJa.get(), rBValues.get(), bNumValuesPerRow );
        }
        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < cNumValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
}

BOOST_AUTO_TEST_CASE( matrixAddSizesTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<ELLKernelTrait::matrixAddSizes > matrixAddSizes;

    ContextPtr loc         = testContext;
    matrixAddSizes.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    IndexType valuesAIa[] =
    { 2, 3, 2, 3, 4 };
    const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
    IndexType valuesAJa[] =
    { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
    const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
    IndexType valuesBIa[] =
    { 2, 2, 2, 3, 3 };
    const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
    IndexType valuesBJa[] =
    { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
    const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
    IndexType expectedCIa[] =
    { 4, 3, 3, 4, 4 };
    const IndexType expectedCNumRows = sizeof( expectedCIa ) / sizeof( IndexType );
    // for matrix add A and B must have same shape
    BOOST_REQUIRE_EQUAL( aNumRows, bNumRows );
    IndexType cNumRows = aNumRows; // C gets same shape as A and B
    BOOST_REQUIRE_EQUAL( cNumRows, expectedCNumRows );
    // values per row needed, verify that numValues is multiple of numRows
    IndexType aNumValuesPerRow = aNumValues / aNumRows;
    BOOST_REQUIRE_EQUAL( aNumRows * aNumValuesPerRow, aNumValues );
    IndexType bNumValuesPerRow = bNumValues / bNumRows;
    BOOST_REQUIRE_EQUAL( bNumRows * bNumValuesPerRow, bNumValues );
    HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
    HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
    HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
    HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
    HArray<IndexType> CIa( testContext );
    {
        ReadAccess<IndexType> rAIa( AIa, loc );
        ReadAccess<IndexType> rAJa( AJa, loc );
        ReadAccess<IndexType> rBIa( BIa, loc );
        ReadAccess<IndexType> rBJa( BJa, loc );
        WriteOnlyAccess<IndexType> wCIa( CIa, loc, cNumRows );
        SCAI_CONTEXT_ACCESS( loc );
        bool diagonalProperty = false;
        IndexType numColumns = aNumRows; // square matrices here
        matrixAddSizes[loc]( wCIa.get(), aNumRows, numColumns, diagonalProperty, rAIa.get(), rAJa.get(), aNumValuesPerRow,
                             rBIa.get(), rBJa.get(), bNumValuesPerRow );
    }
    ReadAccess<IndexType> rCIa( CIa );

    for ( IndexType i = 0; i < cNumRows; i++ )
    {
        BOOST_CHECK_EQUAL( expectedCIa[i], rCIa[i] );
    }
}

BOOST_AUTO_TEST_CASE( matrixAddTest )
{
    typedef DefaultReal ValueType;
    ContextPtr testContext = ContextFix::testContext;
    LAMAKernel<ELLKernelTrait::matrixAdd<ValueType> > matrixAdd;

    ContextPtr loc         = testContext;
    matrixAdd.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // Check with neutral beta
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType valuesCIa[] =
        { 4, 3, 3, 4, 4 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 3, 9, 11, 12, 3, 1, 15, 3, 9, 12, 8, 8, 7, 12, 15, 3, 0, 0, 9, 7 };
        IndexType expectedCJa[] =
        { 0, 0, 1, 0, 0, 1, 3, 2, 1, 2, 2, 4, 3, 2, 3, 3, 0, 0, 3, 4 };
        IndexType cNumValues = 20;
        IndexType numColumns = 5; // for convenience
        ValueType alpha = 1;
        ValueType beta = 1;
        // input arrays, directly initialized on context device
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays, CValues, CJa
        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );
        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = cNumValues / cNumRows;
        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, cNumValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, cNumValues );
            SCAI_CONTEXT_ACCESS( loc );
            bool diagonalProperty = false; // does not matter here
            matrixAdd[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, diagonalProperty,
                            alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow, beta, rBIa.get(), rBJa.get(),
                            rBValues.get(), bNumValuesPerRow );

        }

        // sort the columns, otherwise comparison might fail

        {
            ReadAccess<IndexType> rCIA( CIa );
            WriteAccess<IndexType> wCJA( CJa );
            WriteAccess<ValueType> wCValues( CValues );

            OpenMPELLUtils::sortRowElements( wCJA.get(), wCValues.get(), rCIA.get(), cNumRows, cNumValuesPerRow, false );
        }

        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < cNumValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
        }
    }
    // Check with set beta
    {
        ValueType valuesAValues[] =
        { 1, 5, 2, 4, 3, 3, 7, 3, 7, 9, 0, 8, 0, 9, 8, 0, 0, 0, 0, 7 };
        const IndexType nAValues = sizeof( valuesAValues ) / sizeof( ValueType );
        IndexType valuesAIa[] =
        { 2, 3, 2, 3, 4 };
        const IndexType aNumRows = sizeof( valuesAIa ) / sizeof( IndexType );
        IndexType valuesAJa[] =
        { 1, 0, 1, 0, 0, 3, 3, 2, 2, 2, 0, 4, 0, 3, 3, 0, 0, 0, 0, 4 };
        const IndexType aNumValues = sizeof( valuesAJa ) / sizeof( IndexType );
        ValueType valuesBValues[] =
        { 3, 4, 9, 8, 3, 8, 8, 7, 9, 7, 0, 0, 0, 5, 0 };
        const IndexType nBValues = sizeof( valuesBValues ) / sizeof( ValueType );
        IndexType valuesBIa[] =
        { 2, 2, 2, 3, 2 };
        const IndexType bNumRows = sizeof( valuesBIa ) / sizeof( IndexType );
        IndexType valuesBJa[] =
        { 0, 0, 1, 0, 2, 2, 3, 3, 1, 3, 0, 0, 0, 2, 0 };
        const IndexType bNumValues = sizeof( valuesBJa ) / sizeof( IndexType );
        IndexType valuesCIa[] =
        { 4, 3, 3, 4, 4 };
        const IndexType cNumRows = sizeof( valuesCIa ) / sizeof( IndexType );
        ValueType expectedCValues[] =
        { 6, 13, 20, 20, 3, 1, 23, 3, 18, 15, 16, 8, 14, 17, 22, 3, 0, 0, 9, 7 };
        IndexType expectedCJa[] =
        { 0, 0, 1, 0, 0, 1, 3, 2, 1, 2, 2, 4, 3, 2, 3, 3, 0, 0, 3, 4 };
        IndexType cNumValues = 20;
        IndexType numColumns = 5; // for convenience
        ValueType alpha = 1;
        ValueType beta = 2;
        // input arrays, directly initialized on context device
        HArray<ValueType> AValues( nAValues, valuesAValues, testContext );
        HArray<IndexType> AIa( aNumRows, valuesAIa, testContext );
        HArray<IndexType> AJa( aNumValues, valuesAJa, testContext );
        HArray<ValueType> BValues( nBValues, valuesBValues, testContext );
        HArray<IndexType> BIa( bNumRows, valuesBIa, testContext );
        HArray<IndexType> BJa( bNumValues, valuesBJa, testContext );
        HArray<IndexType> CIa( cNumRows, valuesCIa, testContext );
        // output arrays

        HArray<ValueType> CValues( testContext );
        HArray<IndexType> CJa( testContext );

        IndexType aNumValuesPerRow = aNumValues / aNumRows;
        IndexType bNumValuesPerRow = bNumValues / bNumRows;
        IndexType cNumValuesPerRow = cNumValues / cNumRows;

        {
            ReadAccess<IndexType> rAIa( AIa, loc );
            ReadAccess<IndexType> rAJa( AJa, loc );
            ReadAccess<ValueType> rAValues( AValues, loc );
            ReadAccess<IndexType> rBIa( BIa, loc );
            ReadAccess<IndexType> rBJa( BJa, loc );
            ReadAccess<ValueType> rBValues( BValues, loc );
            ReadAccess<IndexType> rCIa( CIa, loc );
            WriteOnlyAccess<ValueType> wCValues( CValues, loc, cNumValues );
            WriteOnlyAccess<IndexType> wCJa( CJa, loc, cNumValues );
            SCAI_CONTEXT_ACCESS( loc );
            bool diagonalProperty = false; // does not matter here
            matrixAdd[loc]( wCJa.get(), wCValues.get(), rCIa.get(), cNumValuesPerRow, aNumRows, numColumns, diagonalProperty,
                            alpha, rAIa.get(), rAJa.get(), rAValues.get(), aNumValuesPerRow, beta, rBIa.get(), rBJa.get(),
                            rBValues.get(), bNumValuesPerRow );

        }
        // sort the columns, otherwise comparison might fail

        {
            ReadAccess<IndexType> rCIA( CIa );
            WriteAccess<IndexType> wCJA( CJa );
            WriteAccess<ValueType> wCValues( CValues );

            OpenMPELLUtils::sortRowElements( wCJA.get(), wCValues.get(), rCIA.get(), cNumRows, cNumValuesPerRow, false );
        }

        ReadAccess<ValueType> rCValues( CValues );
        ReadAccess<IndexType> rCJa( CJa );

        for ( IndexType i = 0; i < cNumValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedCValues[i], rCValues[i] );
            BOOST_CHECK_EQUAL( expectedCJa[i], rCJa[i] );
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

    const IndexType numRows = 3;

    HArray<IndexType> ellIA( { 2, 2, 1 },  testContext );
    HArray<IndexType> ellJA( { 0, 0, 2, 2, 1, invalidIndex }, testContext );

    HArray<IndexType> row;   // result for rowIndexes
    HArray<IndexType> pos;   // result for positions

    IndexType columnIndex = 1;   // has 1 entry

    ELLUtils::getColumnPositions( row, pos, ellIA, ellJA, columnIndex, testContext );

    BOOST_TEST( hostReadAccess( row ) == std::vector<IndexType>( { 1 } ), per_element() );
    BOOST_TEST( hostReadAccess( pos ) == std::vector<IndexType>( { 4 } ), per_element() );

    columnIndex = 2;

    ELLUtils::getColumnPositions( row, pos, ellIA, ellJA, columnIndex, testContext );

    BOOST_REQUIRE_EQUAL( row.size(), IndexType( 2 ) );   //  two entries for column 2, order might be arbitrary

    // Verifiy: ellJA[ pos[ i ] ] == j, row[i] == pos[i] % numRows

    HArray<IndexType> ja;
    HArrayUtils::gather( ja, ellJA, pos, common::BinaryOp::COPY, testContext );
    BOOST_TEST( hostReadAccess( ja ) == std::vector<IndexType>( ja.size(), columnIndex ), per_element() );

    // Verifiy: row[i] == pos[i] % numRows

    HArrayUtils::compute( pos, pos, common::BinaryOp::MODULO, numRows, testContext );
    BOOST_TEST( hostReadAccess( pos ) == hostReadAccess( row ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    SCAI_ASSERT_EQ_ERROR( ellIA.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), numRows * numValuesPerRow, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellValues.size(), numRows * numValuesPerRow, "size mismatch" )

    HArray<ValueType> x( { 3, -3, 2, -2 },  testContext );
    HArray<ValueType> y( { 1, -1, 2, -2, 1, 1, -1 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numColumns, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, y.size(), "size mismatch" );

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        HArray<ValueType> res( testContext );

        SCAI_LOG_DEBUG( logger, "compute res = " << alpha << " * ELL * x + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", ELL: ia = " << ellIA << ", ja = " << ellJA << ", values = " << ellValues )

        common::MatrixOp op = common::MatrixOp::NORMAL;

        ELLUtils::gemv( res, alpha, x, beta, y, 
                        numRows, numColumns, numValuesPerRow, 
                        ellIA, ellJA, ellValues,
                        op, false, testContext );

        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gevmTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();
    ContextPtr loc         = testContext;

    static LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    SCAI_ASSERT_EQ_ERROR( ellIA.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), numValuesPerRow * numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellValues.size(), numValuesPerRow * numRows, "size mismatch" )

    HArray<ValueType> x( { 3, -2, -2, 3, 1, 0, 1 },  testContext );
    HArray<ValueType> y( { 1, -1, 2, -2 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numRows, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, y.size(), "size mismatch" );

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        HArray<ValueType> res( testContext );

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x * ELL + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", ELL: ia = " << ellIA << ", ja = " << ellJA << ", values = " << ellValues )

        common::MatrixOp op = common::MatrixOp::TRANSPOSE;

        ELLUtils::gemv( res, alpha, x, beta, y, 
                        numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues,
                        op, false, testContext );

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEMVTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;
    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    data1::getRowIndexes( rowIndexes );

    SCAI_ASSERT_EQ_ERROR( numRows, ellIA.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( numRows * numValuesPerRow, ellJA.size(), "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), ellValues.size(), "size mismatch" )

    HArray<ValueType> x( { 3, -3, 2, -2 },  testContext );
    HArray<ValueType> y( { 1, -1, 2, -2, 1, 1, -1 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numColumns, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, y.size(), "size mismatch" );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha; ++icase )
    {
        ValueType alpha = alpha_values[icase];

        HArray<ValueType> res( y );

        SCAI_LOG_INFO( logger, "compute res += " << alpha << " * CSR * x "
                       << ", with x = " << x
                       << ", CSR: ia = " << ellIA << ", ja = " << ellJA << ", values = " << ellValues )


        auto op = common::MatrixOp::NORMAL;

        bool async = false;

        ELLUtils::gemvSp(  res, alpha, x, 
                           numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues,
                           op, rowIndexes, async, testContext );

        HArray<ValueType> expectedRes( y );

        ValueType beta = 1;  // res = alpha * A * x + 1 * res <-> res += alpha * A * x

        expectedRes = data1::getGEMVNormalResult( alpha, x, beta, expectedRes );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( spGEVMTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;
    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    data1::getRowIndexes( rowIndexes );

    SCAI_ASSERT_EQ_ERROR( ellIA.size(), numRows, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), numRows * numValuesPerRow, "size mismatch" )
    SCAI_ASSERT_EQ_ERROR( ellValues.size(), ellJA.size(), "size mismatch" )

    HArray<ValueType> x( { 3, -2, -2, 3, 1, 0, 1 }, testContext );
    HArray<ValueType> y( { 1, -1, 2, -2 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numRows, x.size(), "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, y.size(), "size mismatch" );

    // use different alpha and beta values as kernels might be optimized for it

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha; ++icase )
    {
        ValueType alpha = alpha_values[icase];

        HArray<ValueType> res( y );

        SCAI_LOG_INFO( logger, "compute res += " << alpha << " * CSR * x "
                       << ", with x = " << x
                       << ", CSR: ia = " << ellIA << ", ja = " << ellJA << ", values = " << ellValues )

        auto op = common::MatrixOp::TRANSPOSE;
        bool async = false;

        ELLUtils::gemvSp(  res, alpha, x, numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues,
                           op, rowIndexes, async, testContext );

        ValueType beta = 1;  // res = alpha * x * A + 1 * res <-> res += alpha * x * A

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( res ) == hostReadAccess( expectedRes ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext )

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data2::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    HArray<ValueType> rhs( { 1, -1, 2, -2 }, testContext );
    HArray<ValueType> oldSolution( { 3, -2, -2, 3 }, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> newSolution( testContext );
        
        bool async = true;

        std::unique_ptr<tasking::SyncToken> token;

        token.reset( ELLUtils::jacobi( newSolution, omega, oldSolution, rhs, 
                                       ellIA, ellJA, ellValues, async, testContext ) );
        
        if ( token.get() )
        {
            token->wait();
        }

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

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiHaloTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "jacobiHalo test @ " << *testContext )

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );
    HArray<IndexType> rowIndexes( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;
    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    data1::getRowIndexes( rowIndexes );

    HArray<ValueType> oldSolution( { 3, -2, -2, 3, 1, 0, 2 }, testContext );
    HArray<ValueType> diag( { 9,  8,  7, 6, 7, 8, 9 }, testContext );

    SCAI_ASSERT_EQ_ERROR( oldSolution.size(), numRows, "test size mismatch" )
    SCAI_ASSERT_EQ_ERROR( diag.size(), numRows, "test size mismatch" )

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> solution( numRows, ValueType( 0 ), testContext );
       
        ELLUtils::jacobiHalo( solution, omega, diag, oldSolution,
                             ellIA, ellJA, ellValues, rowIndexes, testContext );

        HArray<ValueType> expectedSol( numRows, ValueType( 0 ), testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        auto maxDiff = HArrayUtils::maxDiffNorm( expectedSol, solution );

        if ( maxDiff >= common::TypeTraits<ValueType>::small() )
        {
            BOOST_TEST( hostReadAccess( expectedSol ) == hostReadAccess( oldSolution ), per_element() );
        }
        else 
        {
            BOOST_CHECK( maxDiff < common::TypeTraits<ValueType>::small() );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxValTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    LAMAKernel<ELLKernelTrait::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = testContext;
    absMaxVal.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<IndexType> ellIA( testContext );
    HArray<IndexType> ellJA( testContext );
    HArray<ValueType> ellValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValuesPerRow;

    data1::getELLTestData( numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues );

    ValueType maxVal = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rIA( ellIA, loc );
        ReadAccess<ValueType> rValues( ellValues, loc );

        maxVal = absMaxVal[loc]( numRows, numValuesPerRow, rIA.get(), rValues.get() );
    }

    ValueType expectedMaxVal = data1::getMaxVal<ValueType>();

    BOOST_CHECK_EQUAL( expectedMaxVal, maxVal );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()
