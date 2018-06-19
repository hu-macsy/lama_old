/**
 * @file DenseUtilsTest.cpp
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
 * @brief Contains tests for the class DenseUtils that run on different devices
 * @author Thomas Brandes
 * @date 19.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/sparsekernel/DenseUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/sparsekernel/test/TestMacros.hpp>

#include <scai/sparsekernel/test/TestData1.hpp>
#include <scai/sparsekernel/test/TestData2.hpp>

/*--------------------------------------------------------------------- */

using namespace scai;

using namespace hmemo;
using namespace sparsekernel;
using namespace utilskernel;

using common::TypeTraits;
using common::BinaryOp;
using common::UnaryOp;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DenseUtilsTest" )

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( nonZeroValuesTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "nonZeroValues test for " << *testContext )

    const ValueType dense_values[] = { 1, 0, 2, -3, 1, 0, 2, 4, -5, 6, 3, 1 };

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    const IndexType dense_n = sizeof( dense_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( dense_n, numRows * numColumns );

    HArray<ValueType> dense( { 1, 0, 2, -3, 1, 0, 2, 4, -5, 6, 3, 1 }, testContext );

    IndexType count = DenseUtils::getNumValues( dense, numRows, numColumns, testContext );

    BOOST_CHECK_EQUAL( count, dense.size() - 2 );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( getCSRTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    HArray<ValueType> dense( { 1, 0, 2,
                              -3, 1, 0,
                               0, 0, -5,
                               6, 0, 1   }, testContext );

    // Note: compress keeps order of values

    HArray<IndexType> expIA( { 2, 2, 1, 2 } );
    HArray<IndexType> expJA( { 0, 2, 0, 1, 2, 0, 2 } );
    HArray<ValueType> expValues( { 1, 2, -3, 1, -5, 6, 1 } );

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    BOOST_REQUIRE_EQUAL( dense.size(), numRows * numColumns );

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    DenseUtils::getSparseRowSizes( csrIA, numRows, numColumns, dense, testContext );
    BOOST_TEST( hostReadAccess( expIA ) == hostReadAccess( csrIA ), per_element() );

    DenseUtils::convertDense2CSR( csrIA, csrJA, csrValues, numRows, numColumns, dense, testContext );

    BOOST_TEST( hostReadAccess( expJA ) == hostReadAccess( csrJA ), per_element() );
    BOOST_TEST( hostReadAccess( expValues ) == hostReadAccess( csrValues ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setCSRTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "setCSRValues test for " << *testContext )

    HArray<IndexType> csrIA( { 0, 2, 4, 5, 7 }, testContext );
    HArray<IndexType> csrJA( { 0, 2, 0, 1, 2, 0, 2 }, testContext );
    HArray<ValueType> csrValues( { 1, 2, -3, 1, -5, 6, 1 }, testContext );

    // expected dense storage data, stored row-wise

    HArray<ValueType> expDense( { 1, 0, 2,
                                 -3, 1, 0,
                                  0, 0, -5,
                                  6, 0, 1   } );

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    BOOST_REQUIRE_EQUAL( csrValues.size(), csrJA.size() );
    BOOST_REQUIRE_EQUAL( csrIA.size(), numRows + 1 );
    BOOST_REQUIRE_EQUAL( expDense.size(), numRows * numColumns );

    HArray<ValueType> dense;

    DenseUtils::convertCSR2Dense( dense, numRows, numColumns, csrIA, csrJA, csrValues, testContext );

    BOOST_TEST( hostReadAccess( expDense ) == hostReadAccess( dense ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setValueTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "setValue test for " << *testContext )

    const ValueType dense_values[] = { 1, 1, 2,
                                       3, 1, 3,
                                       2, 4, 5,
                                       6, 9, 1
                                     };

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    const IndexType dense_n = sizeof( dense_values ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( dense_n, numRows * numColumns );

    ValueType val = 2;

    for ( int op = 0; op < static_cast<int>( BinaryOp::MAX_BINARY_OP ); ++op )
    {
        BinaryOp binOp = BinaryOp( op );

        HArray<ValueType> dense( numRows * numColumns, dense_values, testContext );

        DenseUtils::setScalar( dense, numRows, numColumns, val, binOp, testContext );

        auto rDense = hostReadAccess<ValueType>( dense );

        for ( IndexType i = 0; i < dense_n; ++i )
        {
            ValueType expectedVal = applyBinary( dense_values[i], binOp, val );
            BOOST_CHECK_EQUAL( expectedVal, rDense[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleRowsTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "scaleRows test for " << *testContext )

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    HArray<ValueType> rows( { 1, 2, 0, 1 }, testContext );

    HArray<ValueType> dense( { 1, 1, 2,
                               3, 1, 3,
                               2, 4, 5,
                               6, 9, 1 }, testContext );

    HArray<ValueType> expDense( { 1, 1, 2,
                                  6, 2, 6,
                                  0, 0, 0,
                                  6, 9, 1 } );

    DenseUtils::setRows( dense, numRows, numColumns, rows, common::BinaryOp::MULT, testContext );

    BOOST_TEST( hostReadAccess( expDense ) == hostReadAccess( dense ), per_element() );
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "jacobi test for " << *testContext )

    HArray<ValueType> denseValues( testContext );

    IndexType numRows;
    IndexType numColumns;

    data2::getDenseTestData( numRows, numColumns, denseValues );

    BOOST_CHECK_EQUAL( numRows, numColumns );  // jacobi only for square matrices

    HArray<ValueType> rhs(  { 1, -1, 2, -2 }, testContext );
    HArray<ValueType> oldSolution( { 3, -2, -2, 3 }, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> res( testContext );

        DenseUtils::jacobi( res, omega, oldSolution, rhs, numRows, denseValues, testContext );

        HArray<ValueType> expectedRes( testContext );

        data2::getJacobiResult( expectedRes, oldSolution, omega, rhs );

        auto maxDiff = utilskernel::HArrayUtils::maxDiffNorm( expectedRes, res );

        BOOST_CHECK( maxDiff < 0.001 );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( jacobiHaloTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "jacobi halo test @ " << *testContext )

    HArray<ValueType> denseValues( testContext );

    IndexType numRows;
    IndexType numColumns;

    data1::getDenseTestData( numRows, numColumns, denseValues );

    HArray<ValueType> oldSolution( { 3, -2, -2, 3, 1, 0, 2 }, testContext );
    HArray<ValueType> diag( { 9, 8, 7, 6, 7, 8, 9 }, testContext );

    const ValueType omega_values[] = { 0, 0.5, 0.7, 1 };

    const IndexType n_omega  = sizeof( omega_values ) / sizeof( ValueType );

    for ( IndexType icase = 0; icase < n_omega; ++icase )
    {
        ValueType omega  = omega_values[icase];

        HArray<ValueType> res( numRows, ValueType( 0 ), testContext );

        DenseUtils::jacobiHalo( res, omega,diag, oldSolution, numRows, numColumns, denseValues, testContext );

        HArray<ValueType> expectedRes( numRows, ValueType( 0 ), testContext );

        data1::getJacobiHaloResult( expectedRes, oldSolution, diag, omega );

        auto maxDiff = utilskernel::HArrayUtils::maxDiffNorm( expectedRes, res );
        auto eps = common::TypeTraits<ValueType>::small();

        if ( maxDiff >= eps )
        {
            // this test gives more detailled information where the result values differ from expected ones
            BOOST_TEST( hostReadAccess( expectedRes ) == hostReadAccess( res ), per_element() );
        }
        else
        {
            BOOST_CHECK( maxDiff < eps );
        }
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    SCAI_LOG_INFO( logger, "normalGEMV test @ " << *testContext )

    HArray<ValueType> denseValues( testContext );

    IndexType numRows;
    IndexType numColumns;

    data1::getDenseTestData( numRows, numColumns, denseValues );

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

        DenseUtils::gemv( res, alpha, x, beta, y,
                          numRows, numColumns, denseValues,
                          common::MatrixOp::NORMAL, false, testContext );

        HArray<ValueType> expectedRes = data1::getGEMVNormalResult( alpha, x, beta, y );

        BOOST_TEST( hostReadAccess( expectedRes ) == hostReadAccess( res ), boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END()

