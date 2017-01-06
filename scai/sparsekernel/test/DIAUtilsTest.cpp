/**
 * @file DIAUtilsTest.cpp
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
#include <scai/utilskernel/LArray.hpp>
#include <scai/tasking/SyncToken.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>

#include <scai/sparsekernel/test/TestData1.hpp>
#include <scai/sparsekernel/test/TestData2.hpp>

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

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRTest0, ValueType, scai_numeric_test_types )
{
    // check to get a correct CSR storage from empty DIA storage

    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DIAKernelTrait::getCSRSizes<ValueType> > getCSRSizes;
    static LAMAKernel<DIAKernelTrait::getCSRValues<ValueType, ValueType> > getCSRValues;
    static LAMAKernel<UtilKernelTrait::scan<IndexType> > scan;

    ContextPtr loc = testContext;

    getCSRSizes.getSupportedContext( loc, getCSRValues, scan );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "getCSRSizes/getCSRValues test for " << *testContext << " on " << *loc )

    IndexType numRows = 5;
    IndexType numColumns = 5;
    IndexType numDiagonals = 0;

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    ValueType eps = 0;
    bool diagonalProperty = true;

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
            BOOST_CHECK_EQUAL( rIA[i], 1 );
        }
    }

    IndexType numValues = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<IndexType> wIA( csrIA, loc );
        wIA.resize( numRows + 1 );
        numValues = scan[loc]( wIA.get(), numRows );
    }

    BOOST_REQUIRE_EQUAL( numRows, numValues );

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
            BOOST_CHECK_EQUAL( rJA[i], i );
            BOOST_CHECK_EQUAL( rValues[i], ValueType( 0 ) );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRTest1, ValueType, scai_numeric_test_types )
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

    data1::getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRTest2, ValueType, scai_numeric_test_types )
{
    // Test conversion DIA -> CSR  this time with diagonal property

    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DIAKernelTrait::getCSRSizes<ValueType> > getCSRSizes;
    static LAMAKernel<DIAKernelTrait::getCSRValues<ValueType, ValueType> > getCSRValues;
    static LAMAKernel<UtilKernelTrait::scan<IndexType> > scan;

    ContextPtr loc = testContext;

    getCSRSizes.getSupportedContext( loc, getCSRValues, scan );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "getCSRSizes/getCSRValues test for " << *testContext << " on " << *loc )

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

    data2::getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    ValueType eps = 0;
    bool diagonalProperty = true;

    {
        SCAI_CONTEXT_ACCESS( loc );

        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );
        WriteOnlyAccess<IndexType> wIA( csrIA, loc, numRows );

        getCSRSizes[loc]( wIA.get(), diagonalProperty, numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get(), eps );
    }

    IndexType numValues = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );

        WriteAccess<IndexType> wIA( csrIA, loc );
        wIA.resize( numRows + 1 );
        numValues = scan[loc]( wIA.get(), numRows );
    }

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

    // we compare the computed CSR data with expected CSR data from test data2

    IndexType numRows1;
    IndexType numColumns1;
    IndexType numValues1;

    LArray<IndexType> expIA;
    LArray<IndexType> expJA;
    LArray<ValueType> expValues;

    data2::getCSRTestData( numRows1, numColumns1, numValues1, expIA, expJA, expValues );

    BOOST_REQUIRE_EQUAL( numRows1, numRows );
    BOOST_REQUIRE_EQUAL( numColumns1, numColumns );
    BOOST_REQUIRE_EQUAL( numValues1, numValues );

    BOOST_CHECK_EQUAL( expIA.maxDiffNorm( csrIA ), 0 );
    BOOST_CHECK_EQUAL( expJA.maxDiffNorm( csrJA ), 0 );
    BOOST_CHECK_EQUAL( expValues.maxDiffNorm( csrValues ), 0 );
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

    data1::getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    const ValueType y_values[]   = { 1, -1, 2, -2, 1, 1, -1 };
    const ValueType x_values[]   = { 3, -3, 2, -2 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numColumns, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numRows, n_y, "size mismatch" );

    HArray<ValueType> x( numColumns, x_values, testContext );
    HArray<ValueType> y( numRows, y_values, testContext );

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        HArray<ValueType> res( testContext );

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y )
        {
            common::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
            SCAI_ASYNCHRONOUS( syncToken.get() );

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

        HArray<ValueType> expectedRes;

        data1::getGEMVResult( expectedRes, alpha, x, beta, y );

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

    data1::getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    const ValueType y_values[]   = { 1, -1, 2, -2 };
    const ValueType x_values[]   = { 3, -2, -2, 3, 1, 0, 1 };

    const IndexType n_x   = sizeof( x_values ) / sizeof( ValueType );
    const IndexType n_y   = sizeof( y_values ) / sizeof( ValueType );

    SCAI_ASSERT_EQ_ERROR( numRows, n_x, "size mismatch" );
    SCAI_ASSERT_EQ_ERROR( numColumns, n_y, "size mismatch" );

    HArray<ValueType> x( numRows, x_values, testContext );
    HArray<ValueType> y( numColumns, y_values, testContext );

    const ValueType alpha_values[] = { -3, 1, -1, 0, 2 };
    const ValueType beta_values[]  = { -2, 0, 1 };

    const IndexType n_alpha = sizeof( alpha_values ) / sizeof( ValueType );
    const IndexType n_beta  = sizeof( beta_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_alpha * n_beta; ++icase )
    {
        ValueType alpha = alpha_values[icase % n_alpha ];
        ValueType beta  = beta_values[icase / n_alpha ];

        HArray<ValueType> res( testContext );

        SCAI_LOG_INFO( logger, "compute res = " << alpha << " * x * DIA + " << beta << " * y "
                       << ", with x = " << x << ", y = " << y
                       << ", DIA: offsets = " << diaOffsets << ", values = " << diaValues )
        {
            common::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
            SCAI_ASYNCHRONOUS( syncToken.get() );

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

        HArray<ValueType> expectedRes( testContext );

        data1::getGEVMResult( expectedRes, alpha, x, beta, y );

        {
            ReadAccess<ValueType> rComputed( res, hostContext );
            ReadAccess<ValueType> rExpected( expectedRes, hostContext );

            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( rExpected[j], rComputed[j] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    ContextPtr hostContext = Context::getHostPtr();

    static LAMAKernel<DIAKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = testContext;

    jacobi.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext << " on " << *loc )

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data2::getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

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

            ReadAccess<IndexType> rOffsets( diaOffsets, loc );
            ReadAccess<ValueType> rValues( diaValues, loc );

            ReadAccess<ValueType> rOld( oldSolution, loc );
            ReadAccess<ValueType> rRhs( rhs, loc );
            WriteOnlyAccess<ValueType> wSolution( res, loc, numColumns );

            jacobi[loc]( wSolution.get(),
                         numColumns, numDiagonals,
                         rOffsets.get(), rValues.get(),
                         rOld.get(), rRhs.get(), omega, numRows );

        }

        LArray<ValueType> expectedRes( testContext );

        data2::getJacobiResult( expectedRes, oldSolution, omega, rhs );

        ValueType maxDiff = expectedRes.maxDiffNorm( res );

        BOOST_CHECK( common::Math::real( maxDiff ) < 0.001 );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    LAMAKernel<DIAKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = testContext;
    getValuePos.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    LArray<ValueType> diaValues( testContext );
    LArray<IndexType> diaOffsets( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    LArray<ValueType> dense( numRows * numColumns, 0, testContext );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rOffsets( diaOffsets, loc );

        for ( IndexType i = 0; i < numRows; i++ )
        {
            for ( IndexType j = 0; j < numColumns; j++ )
            {
                IndexType pos = getValuePos[loc]( i, j, numRows, rOffsets.get(), numDiagonals );

                ValueType x   = 0;

                if ( pos != nIndex )
                {
                    x = diaValues[ pos ];
                }

                dense[ i * numColumns + j ] = x;
            }
        }
    }

    LArray<ValueType> expectedDense( testContext );

    data1::getDenseTestData( numRows, numColumns, expectedDense );

    BOOST_CHECK_EQUAL( expectedDense.maxDiffNorm( dense ), ValueType( 0 ) );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxValTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    LAMAKernel<DIAKernelTrait::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = testContext;
    absMaxVal.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    LArray<ValueType> diaValues( testContext );
    LArray<IndexType> diaOffsets( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    ValueType maxVal = 0;

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rOffsets( diaOffsets, loc );
        ReadAccess<ValueType> rValues( diaValues, loc );

        maxVal = absMaxVal[loc]( numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get() );
    }

    ValueType expectedMaxVal = data1::getMaxVal<ValueType>();

    BOOST_CHECK_EQUAL( expectedMaxVal, maxVal );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()

