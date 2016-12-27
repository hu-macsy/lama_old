/**
 * @file COOUtilsTest.cpp
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
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>

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

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( COOUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.COOUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( offsets2iaTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<COOKernelTrait::offsets2ia > offsets2ia;

    ContextPtr loc = testContext;
    offsets2ia.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "offsets2ia test for " << *testContext << " on " << *loc )
    // Test without diagonal property
    {
        const IndexType offsets_values[] =
        { 0, 2, 2, 3, 5 };
        const IndexType ia_values[] =
        { 0, 0, 2, 3, 3 };
        const IndexType numOffsets = sizeof( offsets_values ) / sizeof( IndexType );
        const IndexType numRows = numOffsets - 1;
        const IndexType numValues = sizeof( ia_values ) / sizeof( IndexType );
        // verify that offsets and ia fit
        BOOST_REQUIRE_EQUAL( numValues, offsets_values[numRows] );
        HArray<IndexType> offsets( numOffsets, offsets_values, testContext );
        HArray<IndexType> ia;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        const IndexType numDiagonals = 0;
        {
            WriteOnlyAccess<IndexType> wIA( ia, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            offsets2ia[loc]( wIA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<IndexType> rIA( ia );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( rIA[i], ia_values[i] );
        }
    }
    // Test with diagonal property
    {
        const IndexType offsets_values[] =
        { 0, 2, 5, 7, 9 };
        // Result with diagonals = 0: ia_values[] = { 0, 0, 1, 1, 1, 2, 2, 3, 3 };
        // But with 3 diagonals we get this:
        const IndexType ia_values[] =
        { 0, 1, 2, 0, 1, 1, 2, 3, 3 };
        const IndexType numOffsets = sizeof( offsets_values ) / sizeof( IndexType );
        const IndexType numRows = numOffsets - 1;
        const IndexType numValues = sizeof( ia_values ) / sizeof( IndexType );
        // verify that offsets and ia fit
        BOOST_REQUIRE_EQUAL( numValues, offsets_values[numRows] );
        HArray<IndexType> offsets( numOffsets, offsets_values, testContext );
        HArray<IndexType> ia;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        const IndexType numDiagonals = 3;
        {
            WriteOnlyAccess<IndexType> wIA( ia, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            offsets2ia[loc]( wIA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<IndexType> rIA( ia );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            // SCAI_LOG_TRACE( logger,  "rIA[" << i << "] = " << rIA[i] << ", expects " << ia_values[i] )
            BOOST_CHECK_EQUAL( rIA[i], ia_values[i] );
        }
    }
} // offsets2iaTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setCSRDataTest, ValueType, scai_numeric_test_types )
{
    typedef float CSRValueType;

    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<COOKernelTrait::setCSRData<ValueType, CSRValueType> > setCSRData;

    ContextPtr loc = testContext;
    setCSRData.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "setCSRData< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    // setCSRData is for conversion of CSR storage to COO storage
    // is usually just a copy but has some reordering if diagonal property is required
    // here we test only for csrJA
    {
        const IndexType offsets_values[] = { 0, 2, 5, 7, 9 };
        const CSRValueType csrja_values[]   = { 0, 5, 1, 4, 5, 2, 0, 4, 3 };
        const ValueType cooja_values[]   = { 0, 1, 2, 5, 4, 5, 0, 4, 3 };
        const IndexType numOffsets = sizeof( offsets_values ) / sizeof( IndexType );
        const IndexType numRows = numOffsets - 1;
        const IndexType numDiagonals = 3;
        const IndexType numValues = sizeof( csrja_values ) / sizeof( CSRValueType );
        // verify that offsets and ia fit
        BOOST_REQUIRE_EQUAL( numValues, offsets_values[numRows] );
        BOOST_REQUIRE( numDiagonals <= numRows );
        HArray<IndexType> offsets( numOffsets, offsets_values, testContext );
        HArray<CSRValueType> csrJA( numValues, csrja_values, testContext );
        HArray<ValueType> cooJA;
        ReadAccess<IndexType> rOffsets( offsets, loc );
        ReadAccess<CSRValueType> rCSRJA( csrJA, loc );
        {
            WriteOnlyAccess<ValueType> wCOOJA( cooJA, loc, numValues );
            SCAI_CONTEXT_ACCESS( loc );
            setCSRData[loc]( wCOOJA.get(), rCSRJA.get(), numValues, rOffsets.get(), numRows, numDiagonals );
        }
        ReadAccess<ValueType> rCOOJA( cooJA );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( rCOOJA[i], cooja_values[i] );
        }
    }
} // setCSRData

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<COOKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = testContext;

    getValuePos.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    LArray<IndexType> cooIA( testContext );
    LArray<IndexType> cooJA( testContext );
    LArray<ValueType> cooValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    HArray<ValueType> denseValues( testContext );

    data1::getDenseTestData( numRows, numColumns, denseValues );

    ValueType zero = 0;

    {
        ReadAccess<IndexType> rIa( cooIA, loc );
        ReadAccess<IndexType> rJa( cooJA, loc );

        // comparison is done via accesses on the host

        ReadAccess<ValueType> rValues( cooValues, hostContext );
        ReadAccess<ValueType> rDense( denseValues, hostContext );

        SCAI_CONTEXT_ACCESS( loc );

        for ( IndexType i = 0; i < numRows; i++ )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                IndexType pos = getValuePos[loc]( i, j, rIa.get(), rJa.get(), numValues );

                IndexType k   = i * numColumns + j;

                if ( pos == nIndex )
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

    LAMAKernel<COOKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = testContext;
    getValuePosCol.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    //    1.0   -   2.0      
    //    0.5  0.3   -       
    //     -    -   3.0      

    const IndexType ia[] = { 0, 0, 1, 1, 2 };
    const IndexType ja[] = { 0, 2, 0, 1, 2 };

    const IndexType numRows = 3;
    const IndexType numValues = 5;

    HArray<IndexType> cooIA( numValues, ia, testContext );
    HArray<IndexType> cooJA( numValues, ja, testContext );

    HArray<IndexType> row;   // result for rowIndexes 
    HArray<IndexType> pos;   // result for positions

    IndexType cnt;

    IndexType columnIndex = 1;   // has 1 entry

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );
        WriteOnlyAccess<IndexType> wRow( row, loc, numRows );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numRows );
        cnt = getValuePosCol[loc]( wRow.get(), wPos.get(), columnIndex, rIA.get(), numRows, rJA.get(), numValues );
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

        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );
        WriteOnlyAccess<IndexType> wRow( row, loc, numRows );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numRows );
        cnt = getValuePosCol[loc]( wRow.get(), wPos.get(), columnIndex, rIA.get(), numRows, rJA.get(), numValues );
    }

    BOOST_REQUIRE_EQUAL( cnt, IndexType( 2 ) );   //  two entries for column 2, order might be arbitrary

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rRow( row );
        ReadAccess<IndexType> rJA( cooJA );
        ReadAccess<IndexType> rIA( cooIA );

        for ( IndexType k = 0; k < cnt; ++k )
        {
            IndexType p = rPos[k];
            IndexType i = rRow[k];
            BOOST_CHECK_EQUAL( rJA[ p ], columnIndex );
            BOOST_CHECK( rIA[p] == i );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getValuePosRowTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<COOKernelTrait::getValuePosRow> getValuePosRow;

    ContextPtr loc = testContext;
    getValuePosRow.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );   // give warning if other context is selected

    //    1.0   -   2.0      
    //    0.5  0.3   -       
    //     -    -   3.0      

    const IndexType ia[] = { 0, 0, 1, 1, 2 };
    const IndexType ja[] = { 0, 2, 0, 1, 2 };

    const IndexType numCols = 3;
    const IndexType numValues = 5;

    HArray<IndexType> cooIA( numValues, ia, testContext );
    HArray<IndexType> cooJA( numValues, ja, testContext );

    HArray<IndexType> col;   // result for colIndexes 
    HArray<IndexType> pos;   // result for positions

    IndexType cnt;

    IndexType rowIndex = 2;   // has 1 entry

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );
        WriteOnlyAccess<IndexType> wCol( col, loc, numCols );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numCols );
        cnt = getValuePosRow[loc]( wCol.get(), wPos.get(), rowIndex, rIA.get(), numCols, rJA.get(), numValues );
    }

    BOOST_REQUIRE_EQUAL( cnt, IndexType( 1 ) );   //  only one entry for column 1

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rCol( col );

        BOOST_CHECK_EQUAL( IndexType( 2 ), rCol[0] );   // is in entry row
        BOOST_CHECK_EQUAL( IndexType( 4 ), rPos[0] );   // value of for (2,2) is at pos 4
    }

    rowIndex = 1;
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );
        WriteOnlyAccess<IndexType> wCol( col, loc, numCols );
        WriteOnlyAccess<IndexType> wPos( pos, loc, numCols );
        cnt = getValuePosRow[loc]( wCol.get(), wPos.get(), rowIndex, rIA.get(), numCols, rJA.get(), numValues );
    }

    BOOST_REQUIRE_EQUAL( cnt, IndexType( 2 ) );   //  two entries for row 2, order might be arbitrary

    {
        ReadAccess<IndexType> rPos( pos );
        ReadAccess<IndexType> rCol( col );
        ReadAccess<IndexType> rJA( cooJA );
        ReadAccess<IndexType> rIA( cooIA );

        for ( IndexType k = 0; k < cnt; ++k )
        {
            IndexType p = rPos[k];
            IndexType j = rCol[k];
            BOOST_CHECK_EQUAL( rIA[ p ], rowIndex );
            BOOST_CHECK( rJA[p] == j );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( hasDiagonalPropertyTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    LAMAKernel<COOKernelTrait::hasDiagonalProperty> hasDiagonalProperty;

    ContextPtr loc = testContext;

    hasDiagonalProperty.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    LArray<IndexType> cooIA( testContext );
    LArray<IndexType> cooJA( testContext );
    LArray<ValueType> cooValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data1::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    bool okay;

    IndexType numDiagonals = common::Math::min( numRows, numColumns );

    {
        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );

        SCAI_CONTEXT_ACCESS( loc );

        if ( numValues < numDiagonals )
        { 
            okay = false;
        }
        else
        {
            okay = hasDiagonalProperty[loc]( rIA.get(), rJA.get(), numDiagonals );
        }
    }

    BOOST_CHECK( !okay );

    // data set 2 has a square matrix with diagonal entries first

    data2::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    BOOST_REQUIRE_EQUAL( numRows, numColumns );

    {
        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );

        SCAI_CONTEXT_ACCESS( loc );

        if ( numValues < numRows )
        {
            okay = false;
        }
        else
        {
            okay = hasDiagonalProperty[loc]( rIA.get(), rJA.get(), numRows );
        }
    }

    BOOST_CHECK( okay );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<COOKernelTrait::scaleRows<ValueType, ValueType> > scaleRows;

    ContextPtr loc = testContext;

    scaleRows.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "scaleRows test for " << *testContext << " on " << *loc )

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

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<ValueType> wValues( cooValues, loc );
        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<ValueType> rRows( rows, loc );

        scaleRows[loc]( wValues.get(), rRows.get(), rIA.get(), numValues );
    }

    // prove by hand on host

    {
        ReadAccess<IndexType> rIA( cooIA, hostContext );
        ReadAccess<ValueType> rRows( rows, hostContext );
        ReadAccess<ValueType> rSavedValues( savedValues, hostContext );
        ReadAccess<ValueType> rValues( cooValues, hostContext );

        for ( IndexType k = 0; k < numValues; ++k )
        {
            ValueType f = rRows[ rIA[k] ];
            BOOST_CHECK_EQUAL( f * rSavedValues[k], rValues[k] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<COOKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

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

    ValueType alpha = 1;
    ValueType beta  = -1;

    const ValueType y_values[]   = { 1, -1, 2, -2, 1, 1, -1 };
    const ValueType x_values[]   = { 3, -3, 2, -2 };
    const ValueType res_values[] = { 9, 22, 8, -13, 3, -1, -6 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );
    const IndexType n_res = sizeof( res_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numColumns, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, n_y, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, n_res, "size mismatch" );

    HArray<ValueType> x( numColumns, x_values, testContext );
    HArray<ValueType> y( numRows, y_values, testContext );

    HArray<ValueType> res( testContext );

    SCAI_LOG_INFO( logger, "compute res = " << alpha << " * COO * x + " << beta << " * y "
                            << ", with x = " << x << ", y = " << y
                            << ", COO: ia = " << cooIA << ", ja = " << cooJA << ", values = " << cooValues )
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );
        ReadAccess<ValueType> rValues( cooValues, loc );

        ReadAccess<ValueType> rX( x, loc );
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( res, loc, numRows );

        normalGEMV[loc]( wResult.get(),
                         alpha, rX.get(), beta, rY.get(),
                         numRows, numValues, rIA.get(), rJA.get(), rValues.get() );
    }

    {
        ReadAccess<ValueType> rResult( res, hostContext );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( rResult[i], res_values[i] );
        }
    }
} 

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gevmTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<COOKernelTrait::normalGEVM<ValueType> > normalGEVM;

    ContextPtr loc = testContext;

    normalGEVM.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

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

    ValueType alpha = 1;
    ValueType beta  = 1;

    const ValueType y_values[]   = { 1, -1, 2, -2 };
    const ValueType x_values[]   = { 3, -2, -2, 3, 1, 0, 1 };
    const ValueType res_values[] = { 13, 15, -16, 14 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );
    const IndexType n_res = sizeof( res_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numRows, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, n_y, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, n_res, "size mismatch" );

    HArray<ValueType> x( numRows, x_values, testContext );
    HArray<ValueType> y( numColumns, y_values, testContext );
    HArray<ValueType> res( testContext );

    SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x * COO + " << beta << " * y "
                            << ", with x = " << x << ", y = " << y
                            << ", COO: ia = " << cooIA << ", ja = " << cooJA << ", values = " << cooValues )
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( cooIA, loc );
        ReadAccess<IndexType> rJA( cooJA, loc );
        ReadAccess<ValueType> rValues( cooValues, loc );

        ReadAccess<ValueType> rX( x, loc );
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( res, loc, numColumns );

        normalGEVM[loc]( wResult.get(),
                         alpha, rX.get(), beta, rY.get(),
                         numColumns, numValues, rIA.get(), rJA.get(), rValues.get() );
    }

    {
        ReadAccess<ValueType> rResult( res, hostContext );

        for ( IndexType i = 0; i < numColumns; ++i )
        {
            BOOST_CHECK_EQUAL( rResult[i], res_values[i] );
        }
    }
} 

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<COOKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = testContext;

    jacobi.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext << " on " << *loc )

    HArray<IndexType> cooIA( testContext );
    HArray<IndexType> cooJA( testContext );
    HArray<ValueType> cooValues( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numValues;

    data2::getCOOTestData( numRows, numColumns, numValues, cooIA, cooJA, cooValues );

    const ValueType rhs_values[]   = { 1, -1, 2, -2 };
    const ValueType old_values[]   = { 3, -2, -2, 3 };

    HArray<ValueType> rhs( numRows, rhs_values, testContext );
    HArray<ValueType> oldSolution( numRows, old_values, testContext );

    // const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const ValueType omega_values[] = { 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> res( testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rIA( cooIA, loc );
            ReadAccess<IndexType> rJA( cooJA, loc );
            ReadAccess<ValueType> rValues( cooValues, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<ValueType> rRhs( rhs, loc );
            WriteOnlyAccess<ValueType> wSolution( res, loc, numColumns );

            jacobi[loc]( wSolution.get(), 
                         numValues, rIA.get(), rJA.get(), rValues.get(),
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

BOOST_AUTO_TEST_SUITE_END()

