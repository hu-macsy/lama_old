/**
 * @file DIAUtilsTest.cpp
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
 * @brief Contains tests for the DIAUtils interface to be tested on different devices
 * @author Thomas Brandes
 * @date 15.12.2016 
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/DIAKernelTrait.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>
#include <scai/sparsekernel/test/TestData1.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;

using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;

using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DIAUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DIAUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DIAKernelTrait::getCSRSizes<ValueType> > getCSRSizes;
    static LAMAKernel<DIAKernelTrait::getCSRValues<ValueType, ValueType> > getCSRValues;
    static LAMAKernel<UtilKernelTrait::scan<IndexType> > scan;

    ContextPtr loc = testContext;

    getCSRSizes.getSupportedContext( loc, getCSRValues, scan );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "getCSRSizes/getCSRValues test for " << *testContext << " on " << *loc )

    const IndexType ia_values[]  = { 2,    1, 2,    3,       2,    0, 2 };
    const IndexType ja_values[]  = { 0, 3, 0, 2, 3, 3, 0, 1, 0, 3,    1, 3 };
    const ValueType csr_values[] = { 6, 4, 7, 9, 4, 3, 2, 5, 2, 1,    1, 2 };

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

    getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    ValueType eps = 0;
    bool diagonalProperty = false;

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );
        WriteOnlyAccess<IndexType> wIA( csrIA, loc, numRows );

        getCSRSizes[loc]( wIA.get(), diagonalProperty, numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get(), eps );
    }

    {
        ReadAccess<IndexType> rIA( csrIA, hostContext );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( rIA[i], ia_values[i] );
        }
    }

    IndexType numValues = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<IndexType> wIA( csrIA, loc );
        wIA.resize( numRows + 1 );
        numValues = scan[loc]( wIA.get(), numRows );
    }

    BOOST_REQUIRE_EQUAL( 12, numValues );

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rIA( csrIA, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );
        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        WriteOnlyAccess<IndexType> wJA( csrJA, loc, numValues );
        WriteOnlyAccess<ValueType> wValues( csrValues, loc, numValues );

        getCSRValues[loc]( wJA.get(), wValues.get(), rIA.get(), diagonalProperty, 
                           numRows, numColumns, numDiagonals,
                           rOffsets.get(), rValues.get(), eps );
    }

    {
        ReadAccess<IndexType> rJA( csrJA, hostContext );
        ReadAccess<ValueType> rValues( csrValues, hostContext );

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( rJA[i], ja_values[i] );
            BOOST_CHECK_EQUAL( rValues[i], csr_values[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "normalGEMV test for " << *testContext << " on " << *loc )

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

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

    SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x + " << beta << " * y " 
                            << ", with x = " << x << ", y = " << y )
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );

        ReadAccess<ValueType> rX( x, loc );
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( res, loc, numRows );

        normalGEMV[loc]( wResult.get(), 
                         alpha, rX.get(), beta, rY.get(), 
                         numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get() );

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

    static LAMAKernel<DIAKernelTrait::normalGEVM<ValueType> > normalGEVM;

    ContextPtr loc = testContext;

    normalGEVM.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "normalGEVM test for " << *testContext << " on " << *loc )

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

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

    SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x * DIA + " << beta << " * y " 
                            << ", with x = " << x << ", y = " << y 
                            << ", DIA: offsets = " << diaOffsets << ", values = " << diaValues )
    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );

        ReadAccess<ValueType> rX( x, loc );
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( res, loc, numColumns );

        normalGEVM[loc]( wResult.get(), 
                         alpha, rX.get(), beta, rY.get(), 
                         numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get() );

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

BOOST_AUTO_TEST_SUITE_END()

