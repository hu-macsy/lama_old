/**
 * @file BLAS2Test.cpp
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
 * @brief Contains tests for the blas2 methods.
 * @author Bea Hornef
 * @date 15.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/hmemo.hpp>
#include <scai/kregistry/KernelContextFunction.hpp>

#include <scai/common/TypeTraits.hpp>

#include <scai/blaskernel/test/TestMacros.hpp>

using namespace scai;
using namespace scai::hmemo;
using common::TypeTraits;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( BLAS2Test )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.BLAS2Test" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
    ContextPtr loc = Context::getContextPtr( gemv.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "gemv< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    const CBLAS_TRANSPOSE trans = CblasNoTrans;   // stands for exactly gemv

    ValueType matrix_row_vals[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };   // row-major
    ValueType matrix_col_vals[] = { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };   // col-major

    ValueType x_vals[]      = { 2.0, -1.0, 4.0 };
    ValueType y_vals[]      = { 10.0, -20.0, 30.0 };

    const IndexType n_matrix = sizeof( matrix_row_vals ) / sizeof( ValueType );
    const IndexType n_x      = sizeof( x_vals ) / sizeof( ValueType );
    const IndexType n_y      = sizeof( y_vals ) / sizeof( ValueType );

    const IndexType m = 2;
    const IndexType n = 3;
    const ValueType alpha = 17.0;
    const IndexType incX = 1;
    const ValueType beta = 13.0;
    const ValueType result_vals[] = { -74, -617, 33 };

    HArray<ValueType> x( n_x, x_vals, testContext );

    // for each of CblasRowMajor, CblasColMajor

    for ( int n_order = 0; n_order < 2; ++ n_order )
    {
        const CBLAS_ORDER order       = n_order == 0 ? CblasRowMajor : CblasColMajor;
        const ValueType*  matrix_vals = order == CblasRowMajor ? matrix_row_vals : matrix_col_vals;
        const IndexType   lda         = order == CblasRowMajor ? n : m ;

        HArray<ValueType> matrix( n_matrix , matrix_vals, testContext );

        for ( IndexType incY = 1; incY <= 2; ++incY )
        {
            HArray<ValueType> y( n_y, y_vals, testContext );

            {
                SCAI_CONTEXT_ACCESS( loc );
                ReadAccess<ValueType> rMatrix( matrix, loc );
                ReadAccess<ValueType> rX( x, loc );
                WriteAccess<ValueType> wY( y, loc );
                gemv[loc->getType()]( order, trans, m, n, alpha, rMatrix.get(), lda, rX.get(), incX, beta, wY.get(), incY );
            }
            {
                ReadAccess<ValueType> rY( y );
                BOOST_CHECK_EQUAL( result_vals[0], rY[0] );
                BOOST_CHECK_EQUAL( result_vals[incY], rY[incY] );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( gevmTest, ValueType, blas_test_types )
{
    const CBLAS_TRANSPOSE trans = CblasTrans;  // Note: stands for gevm

    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
    ContextPtr loc = Context::getContextPtr( gemv.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "gemv< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    ValueType matrix_row_vals[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
    ValueType matrix_col_vals[] = { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };

    ValueType x_vals[] = { 2.0, -1.0, 4.0 };
    ValueType y_vals[] = { 10.0, -20.0, 30.0 };
    const IndexType m = 2;
    const IndexType n = 3;
    const ValueType alpha = 17.0;
    const ValueType beta = 13.0;
    const IndexType incY = 1;
    const ValueType result_inc2[] = { 436.0, 148.0, -120.0 };
    const ValueType result_inc1[] = { 96, -277, 390 };

    HArray<ValueType> x( 3, x_vals, testContext );

    for ( int n_order = 0; n_order < 2; ++ n_order )
    {
        const CBLAS_ORDER order       = n_order == 0 ? CblasRowMajor : CblasColMajor;
        const ValueType*  matrix_vals = order == CblasRowMajor ? matrix_row_vals : matrix_col_vals;
        const IndexType   lda         = order == CblasRowMajor ? n : m ;

        HArray<ValueType> matrix( 6, matrix_vals, testContext );

        for ( IndexType incX = 1; incX <= 2; ++incX )
        {
            HArray<ValueType> y( 3, y_vals, testContext );
            {
                SCAI_CONTEXT_ACCESS( loc );
                ReadAccess<ValueType> rMatrix( matrix, loc );
                ReadAccess<ValueType> rX( x, loc );
                WriteAccess<ValueType> wY( y, loc );
                gemv[loc->getType()]( order, trans, m, n, alpha, rMatrix.get(), lda, rX.get(), incX, beta, wY.get(), incY );
            }

            const ValueType* result = incX == 1 ? result_inc1 : result_inc2;

            {
                ReadAccess<ValueType> rY( y );
                BOOST_CHECK_EQUAL( result[0], rY[0] );
                BOOST_CHECK_EQUAL( result[1], rY[1] );
                BOOST_CHECK_EQUAL( result[2], rY[2] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( cgevmTest, ValueType, blas_test_types )
{
    // ToDo: do this tes with complex values

    const CBLAS_TRANSPOSE trans = CblasConjTrans;  

    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
    ContextPtr loc = Context::getContextPtr( gemv.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "gemv< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    ValueType matrix_row_vals[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
    ValueType matrix_col_vals[] = { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };

    ValueType x_vals[] = { 2.0, -1.0, 4.0 };
    ValueType y_vals[] = { 10.0, -20.0, 30.0 };
    const IndexType m = 2;
    const IndexType n = 3;
    const ValueType alpha = 17.0;
    const ValueType beta = 13.0;
    const IndexType incY = 1;
    const ValueType result_inc2[] = { 436.0, 148.0, -120.0 };
    const ValueType result_inc1[] = { 96, -277, 390 };

    HArray<ValueType> x( 3, x_vals, testContext );

    for ( int n_order = 0; n_order < 2; ++ n_order )
    {
        const CBLAS_ORDER order       = n_order == 0 ? CblasRowMajor : CblasColMajor;
        const ValueType*  matrix_vals = order == CblasRowMajor ? matrix_row_vals : matrix_col_vals;
        const IndexType   lda         = order == CblasRowMajor ? n : m ;

        HArray<ValueType> matrix( 6, matrix_vals, testContext );

        for ( IndexType incX = 1; incX <= 2; ++incX )
        {
            HArray<ValueType> y( 3, y_vals, testContext );
            {
                SCAI_CONTEXT_ACCESS( loc );
                ReadAccess<ValueType> rMatrix( matrix, loc );
                ReadAccess<ValueType> rX( x, loc );
                WriteAccess<ValueType> wY( y, loc );
                gemv[loc->getType()]( order, trans, m, n, alpha, rMatrix.get(), lda, rX.get(), incX, beta, wY.get(), incY );
            }

            const ValueType* result = incX == 1 ? result_inc1 : result_inc2;

            {
                ReadAccess<ValueType> rY( y );
                BOOST_CHECK_EQUAL( result[0], rY[0] );
                BOOST_CHECK_EQUAL( result[1], rY[1] );
                BOOST_CHECK_EQUAL( result[2], rY[2] );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
