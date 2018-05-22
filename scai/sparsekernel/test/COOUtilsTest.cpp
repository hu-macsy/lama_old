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
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( offsets2iaTest )
{
    ContextPtr testContext = ContextFix::testContext;

    LAMAKernel<COOKernelTrait::offsets2ia > offsets2ia;

    ContextPtr loc = testContext;
    offsets2ia.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "offsets2ia test for " << *testContext << " on " << *loc )

    HArray<IndexType> csrIA( { 0, 2, 2, 3, 5 }, testContext );

    const IndexType numRows = csrIA.size() - 1;

    IndexType numValues = HArrayUtils::getVal( csrIA, numRows );

    HArray<IndexType> cooIA;  // result array

    COOUtils::convertCSR2COO( cooIA, csrIA, numValues, testContext );

    HArray<IndexType> expectedIA( { 0, 0, 2, 3, 3 } );

    BOOST_TEST( hostReadAccess( cooIA ) == hostReadAccess( expectedIA ), boost::test_tools::per_element() );

} // offsets2iaTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setCSRDataTest, ValueType, scai_numeric_test_types )
{
    typedef DefaultReal CSRValueType;

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

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvNormalTest, ValueType, scai_numeric_test_types )
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

            auto op = common::MatrixOp::NORMAL;

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, rY.get(),
                             numRows, numColumns, numValues, rIA.get(), rJA.get(), rValues.get(), op );
        }

        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( expectedRes ) == hostReadAccess( res ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTransposeTest, ValueType, scai_numeric_test_types )
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

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

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

BOOST_AUTO_TEST_SUITE_END()

