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
#include <scai/sparsekernel/DIAUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/tasking/SyncToken.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>

#include <scai/sparsekernel/test/TestData1.hpp>
#include <scai/sparsekernel/test/TestData2.hpp>

#include <numeric>

/*--------------------------------------------------------------------- */

using namespace scai;

using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;

using boost::test_tools::per_element;

using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DIAUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DIAUtilsTest" )

/* ------------------------------------------------------------------------------------- */

static IndexType op_convert( int i )
{ 
    return static_cast<IndexType>( i ); 
}

BOOST_AUTO_TEST_CASE( convertCSRTest )
{
    typedef DefaultReal ValueType;

    ContextPtr testContext = ContextFix::testContext;

    const IndexType numRows = 7;
    const IndexType numColumns = 4;

    HArray<IndexType> csrIA(     { 0, 2, 3, 5, 8, 10, 10, 12 } );
    HArray<IndexType> csrJA(     { 0, 3,   0, 2,   3,   0, 1, 3,   0, 3,    1, 3 } );
    HArray<ValueType> csrValues( { 6, 4,   7, 9,   4,   2, 5, 3,   2, 1,    1, 2 } );

    std::vector<int> offset_vals( { -5, -4, -3, -2, -1, 0, 1, 3 } );

    HArray<IndexType> expOffset;
    expOffset.resize( offset_vals.size() );

    transform( offset_vals.begin(), offset_vals.end(), 
               hostWriteOnlyAccess( expOffset, offset_vals.size() ).begin(), op_convert );

    // HArray<IndexType> expOffset( { -5, -4, -3, -2, -1, 0, 1, 3 } );

    HArray<ValueType> expValues( { 0, 0, 0, 0, 0, 0, 1,
                                   0, 0, 0, 0, 2, 0, 0,
                                   0, 0, 0, 2, 0, 0, 2,
                                   0, 0, 0, 5, 0, 0, 0,
                                   0, 7, 0, 0, 1, 0, 0,
                                   6, 0, 9, 3, 0, 0, 0,
                                   0, 0, 4, 0, 0, 0, 0,
                                   4, 0, 0, 0, 0, 0, 0 } );

    HArray<IndexType> diaOffset;
    HArray<ValueType> diaValues;

    DIAUtils::convertCSR2DIA( diaOffset, diaValues, 
                              numRows, numColumns, csrIA, csrJA, csrValues, testContext );

    BOOST_TEST( hostReadAccess(	diaOffset ) == hostReadAccess( expOffset ), per_element() );
    BOOST_TEST( hostReadAccess(	diaValues ) == hostReadAccess( expValues ), per_element() );

    HArray<IndexType> newIA;
    HArray<IndexType> newJA;
    HArray<ValueType> newValues;

    DIAUtils::convertDIA2CSR( newIA, newJA, newValues, 
                              numRows, numColumns, diaOffset, diaValues, testContext );

    // Note: conversion DIA -> CSR gives sorted CSR entries

    // BOOST_CHECK( CSRUtils::isSorted( .. )

    BOOST_TEST( hostReadAccess(	newIA ) == hostReadAccess( csrIA ), per_element() );
    BOOST_TEST( hostReadAccess(	newJA ) == hostReadAccess( csrJA ), per_element() );
    BOOST_TEST( hostReadAccess(	newValues ) == hostReadAccess( csrValues ), per_element() );
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

    HArray<ValueType> x( { 3, -3, 2, -2 }, testContext );
    HArray<ValueType> y( { 1, -1, 2, -2, 1, 1, -1 }, testContext );

    SCAI_ASSERT_EQ_ERROR( numColumns, x.size(), "illegally sized x" );
    SCAI_ASSERT_EQ_ERROR( numRows, y.size(), "illegally sized y" );

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
            std::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
            SCAI_ASYNCHRONOUS( syncToken.get() );

            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rOffsets( diaOffsets, loc );
            ReadAccess<ValueType> rValues( diaValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );
            WriteOnlyAccess<ValueType> wResult( res, loc, numRows );

            auto op = common::MatrixOp::NORMAL;

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, rY.get(),
                             numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get(), op );

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

    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = testContext;

    normalGEMV.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    SCAI_LOG_INFO( logger, "normalGEMV (transpose) test for " << *testContext << " on " << *loc )

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
            std::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
            SCAI_ASYNCHRONOUS( syncToken.get() );

            SCAI_CONTEXT_ACCESS( loc );

            ReadAccess<IndexType> rOffsets( diaOffsets, loc );
            ReadAccess<ValueType> rValues( diaValues, loc );

            ReadAccess<ValueType> rX( x, loc );
            ReadAccess<ValueType> rY( y, loc );
            WriteOnlyAccess<ValueType> wResult( res, loc, numColumns );

            auto op = common::MatrixOp::TRANSPOSE;

            normalGEMV[loc]( wResult.get(),
                             alpha, rX.get(), beta, rY.get(),
                             numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get(), op );

        }

        HArray<ValueType> expectedRes = data1::getGEMVTransposeResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( expectedRes ) == hostReadAccess( res ), boost::test_tools::per_element() );
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

    HArray<ValueType> rhs(  { 1, -1, 2, -2 }, testContext );
    HArray<ValueType> oldSolution( { 3, -2, -2, 3 }, testContext );

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

        HArray<ValueType> expectedRes( testContext );

        data2::getJacobiResult( expectedRes, oldSolution, omega, rhs );

        auto maxDiff = HArrayUtils::maxDiffNorm( expectedRes, res );

        BOOST_CHECK( maxDiff < 0.001 );
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

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

    IndexType numRows;
    IndexType numColumns;
    IndexType numDiagonals;

    data1::getDIATestData( numRows, numColumns, numDiagonals, diaOffsets, diaValues );

    HArray<ValueType> dense( numRows * numColumns, ValueType( 0 ), testContext );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<IndexType> rOffsets( diaOffsets, loc );

        for ( IndexType i = 0; i < numRows; i++ )
        {
            for ( IndexType j = 0; j < numColumns; j++ )
            {
                IndexType pos = getValuePos[loc]( i, j, numRows, rOffsets.get(), numDiagonals );

                ValueType x   = 0;

                if ( pos != invalidIndex )
                {
                    x = diaValues[ pos ];
                }

                dense[ i * numColumns + j ] = x;
            }
        }
    }

    HArray<ValueType> expectedDense( testContext );

    data1::getDenseTestData( numRows, numColumns, expectedDense );

    BOOST_TEST( hostReadAccess( expectedDense ) == hostReadAccess( dense ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( absMaxValTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = Context::getContextPtr();

    LAMAKernel<DIAKernelTrait::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = testContext;
    absMaxVal.getSupportedContext( loc );

    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );

    HArray<ValueType> diaValues( testContext );
    HArray<IndexType> diaOffsets( testContext );

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

