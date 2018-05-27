/**
 * @file JDSUtilsTest.cpp
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
 * @brief Contains tests for the class CUDAJDSUtils and OpenMPJDSUtils
 * @author Thomas Brandes
 * @date 16.10.2012
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/kregistry.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel.hpp>
#include <scai/sparsekernel/JDSKernelTrait.hpp>
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
    ContextPtr loc         = testContext;

    LAMAKernel<JDSKernelTrait::getRow<ValueType> > getRow;
    getRow.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    ValueType valuesValues[] =
    { 1, 7, 12, 2, 8, 13, 3, 9, 14, 4, 10, 15, 5, 11, 6 };
    const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
    IndexType valuesJa[] =
    { 0, 1, 5, 2, 3, 7, 4, 5, 12, 6, 7, 15, 8, 9, 10 };
    const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
    IndexType valuesDlg[] =
    { 3, 3, 3, 3, 2, 1 };
    const IndexType nDlg = sizeof( valuesDlg ) / sizeof( IndexType );
    IndexType valuesIlg[] =
    { 6, 5, 4 };
    const IndexType nIlg = sizeof( valuesIlg ) / sizeof( IndexType );
    IndexType valuesPerm[] =
    { 1, 2, 0 };
    const IndexType nPerm = sizeof( valuesPerm ) / sizeof( IndexType );
    const IndexType i = 2;
    const IndexType numColumns = 16;
    const IndexType numRows = 3;
    HArray<ValueType> values( nValues, valuesValues, testContext );
    HArray<IndexType> ja( nJa, valuesJa, testContext );
    HArray<IndexType> dlg( nDlg, valuesDlg, testContext );
    HArray<IndexType> ilg( nIlg, valuesIlg, testContext );
    HArray<IndexType> perm( nPerm, valuesPerm, testContext );
    HArray<ValueType> row;
    ReadAccess<ValueType> rValues( values, loc );
    ReadAccess<IndexType> rJa( ja, loc );
    ReadAccess<IndexType> rDlg( dlg, loc );
    ReadAccess<IndexType> rIlg( ilg, loc );
    ReadAccess<IndexType> rPerm( perm, loc );
    {
        WriteOnlyAccess<ValueType> wRow( row, loc, numColumns );
        SCAI_CONTEXT_ACCESS( loc );
        getRow[loc]( wRow.get(), i, numColumns, numRows, rPerm.get(), rIlg.get(), rDlg.get(), rJa.get(), rValues.get() );
    }

    std::vector<ValueType> expectedValues( { 0, 7, 0, 8, 0, 9, 0, 10, 0, 11, 0, 0, 0, 0, 0, 0 } );
    BOOST_TEST( hostReadAccess( row ) == expectedValues, per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( getRowPositionsTest )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<JDSKernelTrait::getRowPositions> getRowPositions;

    ContextPtr loc = testContext;

    getRowPositions.getSupportedContext( loc );

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

    ReadAccess<IndexType> rDlg( jdsDLG, loc );
    ReadAccess<IndexType> rIlg( jdsILG, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    WriteOnlyAccess<IndexType> wPos( pos, loc, numColumns );

    SCAI_CONTEXT_ACCESS( loc );

    IndexType n = 0;    // count number of entries

    for ( IndexType i = 0; i < numRows; ++i )
    {
        n += getRowPositions[loc]( wPos.get(), i, numRows, rIlg.get(), rDlg.get(), rPerm.get() );
    }

    BOOST_CHECK_EQUAL( jdsJA.size(), n );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setRowTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<JDSKernelTrait::setRow<ValueType> > setRow;
    LAMAKernel<JDSKernelTrait::getRow<ValueType> > getRow;

    ContextPtr loc = testContext;
    getRow.getSupportedContext( loc, setRow );

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

    for ( IndexType i = 0; i < numRows; ++i )
    {
        WriteAccess<ValueType> wValues( jdsValues, loc );
        ReadAccess<IndexType> rJa( jdsJA, loc );
        ReadAccess<IndexType> rDlg( jdsDLG, loc );
        ReadAccess<IndexType> rIlg( jdsILG, loc );
        ReadAccess<IndexType> rPerm( jdsPerm, loc );

        SCAI_CONTEXT_ACCESS( loc );

        WriteOnlyAccess<ValueType> wRow( row, loc, numColumns );

        common::BinaryOp op = common::BinaryOp::SUB;

        getRow[loc]( wRow.get(), i, numColumns, numRows, rPerm.get(), rIlg.get(), rDlg.get(), rJa.get(), wValues.get() );
        setRow[loc]( wValues.get(), i, numColumns, numRows, rPerm.get(), rIlg.get(), rDlg.get(), rJa.get(), wRow.get(), op );
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

    LAMAKernel<JDSKernelTrait::scaleRows<ValueType> > scaleRows;

    ContextPtr loc = testContext;
    scaleRows.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    ValueType valuesValues[] = { 1, 7, 12, 2, 8, 13, 3, 9, 14, 4, 10, 15, 5, 11, 6 };
    IndexType valuesDlg[]    = { 3, 3, 3, 3, 2, 1 };
    IndexType valuesIlg[]    = { 6, 5, 4 };
    IndexType valuesPerm[]   = { 1, 2, 0 };
    ValueType valuesDiagonal[] = { 3, 1, 2 };
    ValueType expectedValues[] = { 1, 14, 36, 2, 16, 39, 3, 18, 42, 4, 20, 45, 5, 22, 6 };
    const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
    const IndexType nDlg = sizeof( valuesDlg ) / sizeof( IndexType );
    const IndexType nIlg = sizeof( valuesIlg ) / sizeof( IndexType );
    const IndexType nPerm = sizeof( valuesPerm ) / sizeof( IndexType );
    const IndexType nDiagonal = sizeof( valuesDiagonal ) / sizeof( ValueType );
    const IndexType numRows = 3;
    HArray<ValueType> values( nValues, valuesValues, testContext );
    HArray<IndexType> dlg( nDlg, valuesDlg, testContext );
    HArray<IndexType> ilg( nIlg, valuesIlg, testContext );
    HArray<IndexType> perm( nPerm, valuesPerm, testContext );
    HArray<ValueType> diagonal( nDiagonal, valuesDiagonal, testContext );
    {
        ReadAccess<IndexType> rDlg( dlg, loc );
        ReadAccess<IndexType> rIlg( ilg, loc );
        ReadAccess<IndexType> rPerm( perm, loc );
        ReadAccess<ValueType> rDiagonal( diagonal, loc );
        WriteAccess<ValueType> wValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        scaleRows[loc]( wValues.get(), numRows, rPerm.get(), rIlg.get(), rDlg.get(), rDiagonal.get() );
    }
    {
        ReadAccess<ValueType> rValues( values );

        for ( IndexType i = 0; i < nValues; i++ )
        {
            BOOST_CHECK_EQUAL( expectedValues[i], rValues[i] );
        }
    }
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

BOOST_AUTO_TEST_CASE_TEMPLATE( setCSRValuesTest, ValueType, scai_numeric_test_types )
{
    typedef DefaultReal OtherValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<JDSKernelTrait::setCSRValues<ValueType, OtherValueType> > setCSRValues;

    ContextPtr loc = testContext;
    setCSRValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    /*
     * Testmatrix:
     * 0 0 5 3 0 0 4 0
     * 3 0 4 0 3 5 0 0
     * 0 2 0 8 7 9 0 5
     * 2 0 3 0 0 0 0 0
     * 5 0 0 7 0 0 0 9
     */
    IndexType valuesJDSDlg[] =
    { 5, 5, 4, 2, 1 };
    const IndexType nJDSDlg = sizeof( valuesJDSDlg ) / sizeof( IndexType );
    IndexType valuesJDSIlg[] =
    { 5, 4, 3, 3, 2 };
    const IndexType nJDSIlg = sizeof( valuesJDSIlg ) / sizeof( IndexType );
    IndexType valuesJDSPerm[] =
    { 2, 1, 0, 4, 3 };
    const IndexType nJDSPerm = sizeof( valuesJDSPerm ) / sizeof( IndexType );
    IndexType valuesCSRIa[] =
    { 0, 3, 7, 12, 14, 17 };
    const IndexType nCSRIa = sizeof( valuesCSRIa ) / sizeof( IndexType );
    IndexType valuesCSRJa[] =
    { 2, 3, 6, 0, 2, 4, 5, 1, 3, 4, 5, 7, 0, 2, 0, 3, 7 };
    const IndexType nCSRJa = sizeof( valuesCSRJa ) / sizeof( IndexType );
    OtherValueType valuesCSRValues[] =
    { 5, 3, 4, 3, 4, 3, 5, 2, 8, 7, 9, 5, 2, 3, 5, 7, 9 };
    const IndexType nCSRValues = sizeof( valuesCSRValues ) / sizeof( OtherValueType );
    IndexType expectedJDSJa[] =
    { 1, 0, 2, 0, 0, 3, 2, 3, 3, 2, 4, 4, 6, 7, 5, 5, 7 };
    ValueType expectedJDSValues[] =
    { 2, 3, 5, 5, 2, 8, 4, 3, 7, 3, 7, 3, 4, 9, 9, 5, 5 };

    const IndexType numRows = 5;
    const IndexType nJDS = nCSRValues;
    HArray<IndexType> JDSDlg( nJDSDlg, valuesJDSDlg );
    HArray<IndexType> JDSIlg( nJDSIlg, valuesJDSIlg );
    HArray<IndexType> JDSPerm( nJDSPerm, valuesJDSPerm );
    HArray<IndexType> CSRIa( nCSRIa, valuesCSRIa );
    HArray<IndexType> CSRJa( nCSRJa, valuesCSRJa );
    HArray<OtherValueType> CSRValues( nCSRValues, valuesCSRValues );
    // output arrays
    HArray<IndexType> JDSJa;
    HArray<ValueType> JDSValues;
    {
        WriteOnlyAccess<IndexType> wJDSJa( JDSJa, loc, nJDS );
        WriteOnlyAccess<ValueType> wJDSValues( JDSValues, loc, nJDS );
        ReadAccess<IndexType> rJDSDlg( JDSDlg, loc );
        ReadAccess<IndexType> rJDSIlg( JDSIlg, loc );
        ReadAccess<IndexType> rJDSPerm( JDSPerm, loc );
        ReadAccess<IndexType> rCSRIa( CSRIa, loc );
        ReadAccess<IndexType> rCSRJa( CSRJa, loc );
        ReadAccess<OtherValueType> rCSRValues( CSRValues, loc );
        SCAI_CONTEXT_ACCESS( loc );
        setCSRValues[loc]( wJDSJa.get(), wJDSValues.get(), numRows, rJDSPerm.get(), rJDSIlg.get(), nJDSDlg, rJDSDlg.get(),
                           rCSRIa.get(), rCSRJa.get(), rCSRValues.get() );
    }
    ReadAccess<IndexType> rJDSJa( JDSJa );
    ReadAccess<ValueType> rJDSValues( JDSValues );

    for ( IndexType i = 0; i < nJDS; i++ )
    {
        BOOST_CHECK_EQUAL( expectedJDSJa[i], rJDSJa.get()[i] );
        BOOST_CHECK_EQUAL( expectedJDSValues[i], rJDSValues.get()[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRValuesTest, ValueType, scai_numeric_test_types )
{
    typedef DefaultReal OtherValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<JDSKernelTrait::getCSRValues<ValueType, OtherValueType> > getCSRValues;

    ContextPtr loc = testContext;
    getCSRValues.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );
    /*
     * Testmatrix:
     * 0 0 5 3 0 0 4 0
     * 3 0 4 0 3 5 0 0
     * 0 2 0 8 7 9 0 5
     * 2 0 3 0 0 0 0 0
     * 5 0 0 7 0 0 0 9
     */
    IndexType valuesJDSDlg[] =
    { 5, 5, 4, 2, 1 };
    const IndexType nJDSDlg = sizeof( valuesJDSDlg ) / sizeof( IndexType );
    IndexType valuesJDSIlg[] =
    { 5, 4, 3, 3, 2 };
    const IndexType nJDSIlg = sizeof( valuesJDSIlg ) / sizeof( IndexType );
    IndexType valuesJDSPerm[] =
    { 2, 1, 0, 4, 3 };
    const IndexType nJDSPerm = sizeof( valuesJDSPerm ) / sizeof( IndexType );
    IndexType valuesCSRIa[] =
    { 0, 3, 7, 12, 14, 17 };
    const IndexType nCSRIa = sizeof( valuesCSRIa ) / sizeof( IndexType );
    IndexType valuesJDSJa[] =
    { 1, 0, 2, 0, 0, 3, 2, 3, 3, 2, 4, 4, 6, 7, 5, 5, 7 };
    const IndexType nJDSJa = sizeof( valuesJDSJa ) / sizeof( IndexType );
    ValueType valuesJDSValues[] =
    { 2, 3, 5, 5, 2, 8, 4, 3, 7, 3, 7, 3, 4, 9, 9, 5, 5 };
    const IndexType nJDSValues = sizeof( valuesJDSValues ) / sizeof( ValueType );
    IndexType expectedCSRJa[] =
    { 2, 3, 6, 0, 2, 4, 5, 1, 3, 4, 5, 7, 0, 2, 0, 3, 7 };
    OtherValueType expectedCSRValues[] =
    { 5, 3, 4, 3, 4, 3, 5, 2, 8, 7, 9, 5, 2, 3, 5, 7, 9 };
    const IndexType numRows = 5;
    const IndexType nJDS = nJDSValues;
    // Input arrays: initialized directly on testContext
    HArray<IndexType> JDSJa( nJDSJa, valuesJDSJa, testContext );
    HArray<ValueType> JDSValues( nJDSValues, valuesJDSValues, testContext );
    HArray<IndexType> JDSDlg( nJDSDlg, valuesJDSDlg, testContext );
    HArray<IndexType> JDSIlg( nJDSIlg, valuesJDSIlg, testContext );
    HArray<IndexType> JDSPerm( nJDSPerm, valuesJDSPerm, testContext );
    HArray<IndexType> CSRIa( nCSRIa, valuesCSRIa );
    // Output arrays: no initialization
    HArray<IndexType> CSRJa;
    HArray<OtherValueType> CSRValues;
    {
        ReadAccess<IndexType> rJDSJa( JDSJa, loc );
        ReadAccess<ValueType> rJDSValues( JDSValues, loc );
        ReadAccess<IndexType> rJDSDlg( JDSDlg, loc );
        ReadAccess<IndexType> rJDSIlg( JDSIlg, loc );
        ReadAccess<IndexType> rJDSPerm( JDSPerm, loc );
        ReadAccess<IndexType> rCSRIa( CSRIa, loc );
        WriteOnlyAccess<IndexType> wCSRJa( CSRJa, loc, nJDS );
        WriteOnlyAccess<OtherValueType> wCSRValues( CSRValues, loc, nJDS );
        SCAI_CONTEXT_ACCESS( loc );
        getCSRValues[loc]( wCSRJa.get(), wCSRValues.get(), rCSRIa.get(), numRows, rJDSPerm.get(), rJDSIlg.get(),
                           rJDSDlg.get(), rJDSJa.get(), rJDSValues.get() );
    }
    ReadAccess<IndexType> rCSRJa( CSRJa );
    ReadAccess<OtherValueType> rCSRValues( CSRValues );

    for ( IndexType i = 0; i < nJDS; i++ )
    {
        BOOST_CHECK_EQUAL( expectedCSRJa[i], rCSRJa.get()[i] );
        BOOST_CHECK_EQUAL( expectedCSRValues[i], rCSRValues.get()[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

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

    HArray<ValueType> res( testContext );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * JDS * x + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", JDS: ilg = " << jdsILG << ", ja = " << jdsJA << ", values = " << jdsValues )
        {
            std::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );

            SCAI_ASYNCHRONOUS( syncToken.get() );

            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rPerm( jdsPerm, loc );
            ReadAccess<IndexType> rDLG( jdsDLG, loc );
            ReadAccess<IndexType> rILG( jdsILG, loc );
            ReadAccess<IndexType> rJA( jdsJA, loc );
            ReadAccess<ValueType> rValues( jdsValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );
            WriteOnlyAccess<ValueType> wResult( res, loc, numRows );

            auto op = common::MatrixOp::NORMAL;

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, rY.get(),
                             numRows, numColumns, rPerm.get(), rILG.get(), numDiagonals, rDLG.get(), rJA.get(), rValues.get(), op );
        }

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
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

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
        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rPerm( jdsPerm, loc );
            ReadAccess<IndexType> rDLG( jdsDLG, loc );
            ReadAccess<IndexType> rILG( jdsILG, loc );
            ReadAccess<IndexType> rJA( jdsJA, loc );
            ReadAccess<ValueType> rValues( jdsValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            WriteAccess<ValueType> wResult( res, loc, numRows );

            auto op = common::MatrixOp::NORMAL;

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, wResult.get(),
                             numRows, numColumns, rPerm.get(), rILG.get(), numDiagonals, rDLG.get(), rJA.get(), rValues.get(), op );
        }

        // compare against mv product of dense matrix

        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_CHECK( HArrayUtils::maxDiffNorm( expectedRes, res ) < 0.001 );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTransposeTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

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
        {
            std::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
            SCAI_ASYNCHRONOUS( syncToken.get() );

            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rPerm( jdsPerm, loc );
            ReadAccess<IndexType> rDLG( jdsDLG, loc );
            ReadAccess<IndexType> rILG( jdsILG, loc );
            ReadAccess<IndexType> rJA( jdsJA, loc );
            ReadAccess<ValueType> rValues( jdsValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );
            WriteOnlyAccess<ValueType> wResult( res, loc, numColumns );

            auto op = common::MatrixOp::TRANSPOSE;

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, rY.get(),
                             numRows, numColumns, rPerm.get(), rILG.get(), numDiagonals, rDLG.get(), rJA.get(), rValues.get(), op );

            // implicit synchronizsation at end of this scope
        }

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
        
        JDSUtils::jacobi( newSolution, omega, oldSolution, rhs, 
                          jdsILG, jdsDLG, jdsPerm, jdsJA, jdsValues, testContext );
        
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
