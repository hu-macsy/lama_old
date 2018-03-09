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

BOOST_AUTO_TEST_CASE( getValuePosRowTest )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<JDSKernelTrait::getValuePosRow> getValuePosRow;

    ContextPtr loc = testContext;

    getValuePosRow.getSupportedContext( loc );

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
        n += getValuePosRow[loc]( wPos.get(), i, numRows, rIlg.get(), rDlg.get(), rPerm.get() );
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

BOOST_AUTO_TEST_CASE_TEMPLATE( getValuePosColTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<JDSKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = testContext;

    getValuePosCol.getSupportedContext( loc );

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

        IndexType n = 0;

        {
            ReadAccess<IndexType> rJa( jdsJA, loc );
            ReadAccess<IndexType> rDlg( jdsDLG, loc );
            ReadAccess<IndexType> rIlg( jdsILG, loc );
            ReadAccess<IndexType> rPerm( jdsPerm, loc );

            SCAI_CONTEXT_ACCESS( loc );

            WriteOnlyAccess<IndexType> wRow( row, loc, numColumns );
            WriteOnlyAccess<IndexType> wPos( pos, loc, numColumns );

            n = getValuePosCol[loc]( wRow.get(), wPos.get(), j, numRows, rIlg.get(), rDlg.get(), rPerm.get(), rJa.get() );

        }

        nTotal += n;
    }

    BOOST_CHECK_EQUAL( nTotal, jdsJA.size() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<JDSKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = testContext;
    getValuePos.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    ValueType valuesValues[] =
    { 1, 5, 4, 3, 1, 3, 2, 2, 2, 8, 4, 9, 9, 7, 8, 7, 2 };
    const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
    IndexType valuesJa[] =
    { 0, 0, 2, 2, 0, 1, 3, 5, 4, 1, 3, 7, 6, 6, 9, 9, 9 };
    const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
    IndexType valuesDlg[] =
    { 5, 5, 3, 3, 1 };
    const IndexType nDlg = sizeof( valuesDlg ) / sizeof( IndexType );
    IndexType valuesIlg[] =
    { 5, 4, 4, 2, 2 };
    const IndexType nIlg = sizeof( valuesIlg ) / sizeof( IndexType );
    IndexType valuesPerm[] =
    { 0, 2, 3, 1, 4 };
    const IndexType nPerm = sizeof( valuesPerm ) / sizeof( IndexType );
    ValueType expectedValues[5][10] =
    {
        { 1, 3, 0, 4, 0, 0, 7, 0, 0, 2 },
        { 0, 0, 3, 0, 2, 0, 0, 0, 0, 0 },
        { 5, 0, 0, 2, 0, 0, 0, 9, 0, 8 },
        { 0, 0, 4, 0, 0, 2, 9, 0, 0, 7 },
        { 1, 8, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    const IndexType numColumns = 10;
    const IndexType numRows = 5;
    HArray<ValueType> values( nValues, valuesValues, testContext );
    HArray<IndexType> ja( nJa, valuesJa, testContext );
    HArray<IndexType> dlg( nDlg, valuesDlg, testContext );
    HArray<IndexType> ilg( nIlg, valuesIlg, testContext );
    HArray<IndexType> perm( nPerm, valuesPerm, testContext );
    ReadAccess<ValueType> rValues( values, loc );
    ReadAccess<IndexType> rJa( ja, loc );
    ReadAccess<IndexType> rDlg( dlg, loc );
    ReadAccess<IndexType> rIlg( ilg, loc );
    ReadAccess<IndexType> rPerm( perm, loc );

    IndexType nnz = 0;

    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType j = 0; j < numColumns; j++ )
        {
            SCAI_CONTEXT_ACCESS( loc );
            IndexType pos = getValuePos[loc]( i, j, numRows, rDlg.get(), rIlg.get(), rPerm.get(), rJa.get() );

            if ( pos == invalidIndex )
            {
                BOOST_CHECK_EQUAL( expectedValues[i][j], 0 );
            }
            else
            {
                BOOST_CHECK_EQUAL( expectedValues[i][j], valuesValues[pos] );
                nnz++;
            }
        }
    }

    BOOST_CHECK_EQUAL( nValues, nnz );
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

BOOST_AUTO_TEST_CASE( checkDiagonalPropertyTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<JDSKernelTrait::checkDiagonalProperty> checkDiagonalProperty;

    ContextPtr loc = testContext;
    checkDiagonalProperty.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    // check with matrix without diagonal property
    {
        IndexType valuesJa[]   = { 0, 1, 5, 2, 3, 7, 4, 5, 12, 6, 7, 15, 8, 9, 10 };
        IndexType valuesDlg[]  = { 3, 3, 3, 3, 2, 1 };
        IndexType valuesIlg[]  = { 6, 5, 4 };
        IndexType valuesPerm[] = { 1, 2, 0 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        const IndexType nDlg = sizeof( valuesDlg ) / sizeof( IndexType );
        const IndexType nIlg = sizeof( valuesIlg ) / sizeof( IndexType );
        const IndexType nPerm = sizeof( valuesPerm ) / sizeof( IndexType );
        const IndexType numRows = 3;
        const IndexType numColumns = 16;
        const IndexType numDiagonals = 3;
        HArray<IndexType> ja( nJa, valuesJa, testContext );
        HArray<IndexType> dlg( nDlg, valuesDlg, testContext );
        HArray<IndexType> ilg( nIlg, valuesIlg, testContext );
        HArray<IndexType> perm( nPerm, valuesPerm, testContext );
        ReadAccess<IndexType> rJa( ja, loc );
        ReadAccess<IndexType> rDlg( dlg, loc );
        ReadAccess<IndexType> rIlg( ilg, loc );
        ReadAccess<IndexType> rPerm( perm, loc );
        SCAI_CONTEXT_ACCESS( loc );
        bool diagonalProperty =
            checkDiagonalProperty[loc]( numDiagonals, numRows, numColumns,
                                        rPerm.get(), rJa.get(), rDlg.get() );
        BOOST_CHECK_EQUAL( false, diagonalProperty );
    }
    // check with matrix with diagonal property
    {
        IndexType valuesJa[] =
        { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        IndexType valuesDlg[] =
        { 3, 3, 3 };
        const IndexType nDlg = sizeof( valuesDlg ) / sizeof( IndexType );
        IndexType valuesIlg[] =
        { 3, 3, 3 };
        const IndexType nIlg = sizeof( valuesIlg ) / sizeof( IndexType );
        IndexType valuesPerm[] =
        { 0, 1, 2 };
        const IndexType nPerm = sizeof( valuesPerm ) / sizeof( IndexType );
        const IndexType numRows = 3;
        const IndexType numColumns = 3;
        const IndexType numDiagonals = 3;
        HArray<IndexType> ja( nJa, valuesJa, testContext );
        HArray<IndexType> dlg( nDlg, valuesDlg, testContext );
        HArray<IndexType> ilg( nIlg, valuesIlg, testContext );
        HArray<IndexType> perm( nPerm, valuesPerm, testContext );
        ReadAccess<IndexType> rJa( ja, loc );
        ReadAccess<IndexType> rDlg( dlg, loc );
        ReadAccess<IndexType> rIlg( ilg, loc );
        ReadAccess<IndexType> rPerm( perm, loc );
        SCAI_CONTEXT_ACCESS( loc );
        bool diagonalProperty;
        diagonalProperty = checkDiagonalProperty[loc]( numDiagonals, numRows, numColumns,
                           rPerm.get(), rJa.get(), rDlg.get() );
        BOOST_CHECK_EQUAL( true, diagonalProperty );
    }
    // check with empty matrix
    {
        const IndexType numRows = 0;
        const IndexType numColumns = 0;
        const IndexType numDiagonals = 0;
        HArray<IndexType> ja;
        HArray<IndexType> dlg;
        HArray<IndexType> ilg;
        HArray<IndexType> perm;
        ReadAccess<IndexType> rJa( ja, loc );
        ReadAccess<IndexType> rDlg( dlg, loc );
        ReadAccess<IndexType> rIlg( ilg, loc );
        ReadAccess<IndexType> rPerm( perm, loc );
        SCAI_CONTEXT_ACCESS( loc );
        bool diagonalProperty = checkDiagonalProperty[loc]( numDiagonals, numRows, numColumns,
                                rPerm.get(), rJa.get(), rDlg.get() );
        BOOST_CHECK_EQUAL( false, diagonalProperty );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( ilg2dlgTest )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<JDSKernelTrait::ilg2dlg> ilg2dlg;

    ContextPtr loc = testContext;
    ilg2dlg.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    HArray<IndexType> ilg( { 7, 7, 5, 4, 4, 1 }, testContext );
    const IndexType numRows = ilg.size();

    const IndexType numDiagonals = 7;   // must also be maxval( valuesIlg )

    HArray<IndexType> dlg( testContext );

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteOnlyAccess<IndexType> wDlg( dlg, loc, numDiagonals );
        ReadAccess<IndexType> rIlg( ilg, loc );
        ilg2dlg[loc]( wDlg.get(), numDiagonals, rIlg.get(), numRows );
    }

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
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<JDSKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = testContext;

    jacobi.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext << " on " << *loc )

    HArray<IndexType> jdsPerm( testContext );
    HArray<IndexType> jdsILG( testContext );
    HArray<IndexType> jdsDLG( testContext );
    HArray<IndexType> jdsJA( testContext );
    HArray<ValueType> jdsValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data2::getJDSTestData( numRows, numColumns, numDiagonals, jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues );

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

            ReadAccess<IndexType> rPerm( jdsPerm, loc );
            ReadAccess<IndexType> rDLG( jdsDLG, loc );
            ReadAccess<IndexType> rILG( jdsILG, loc );
            ReadAccess<IndexType> rJA( jdsJA, loc );
            ReadAccess<ValueType> rValues( jdsValues, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<ValueType> rRhs( rhs, loc );
            WriteOnlyAccess<ValueType> wSolution( res, loc, numColumns );

            jacobi[loc]( wSolution.get(), numRows,
                         rPerm.get(), rILG.get(), numDiagonals, rDLG.get(), rJA.get(), rValues.get(),
                         rOld.get(), rRhs.get(), omega );

        }

        HArray<ValueType> expectedRes( testContext );

        data2::getJacobiResult( expectedRes, oldSolution, omega, rhs );

        BOOST_CHECK( HArrayUtils::maxDiffNorm( expectedRes, res ) < 0.1 );

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

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiHaloTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<JDSKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = testContext;

    jacobiHalo.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc, testContext );

    SCAI_LOG_INFO( logger, "jacobiHalo test for " << *testContext << " on " << *loc )

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

        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rPerm( jdsPerm, loc );
            ReadAccess<IndexType> rILG( jdsILG, loc );
            ReadAccess<IndexType> rDLG( jdsDLG, loc );
            ReadAccess<IndexType> rJA( jdsJA, loc );
            ReadAccess<ValueType> rValues( jdsValues, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<ValueType> rDiag( diag, loc );
            WriteAccess<ValueType> wSolution( solution, loc );

            jacobiHalo[loc]( wSolution.get(), numRows, rDiag.get(),
                             numDiagonals, rPerm.get(), rILG.get(), rDLG.get(),
                             rJA.get(), rValues.get(),
                             rOld.get(), omega );
        }

        auto expectedSol = utilskernel::fillHArray<ValueType>( numRows, 0, testContext );

        data1::getJacobiHaloResult( expectedSol, oldSolution, diag, omega );

        // be careful: expectedSol == oldSolution does not hold on GPUs

        BOOST_CHECK( HArrayUtils::maxDiffNorm( expectedSol, solution ) < 0.001 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
