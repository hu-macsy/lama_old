/**
 * @file BLAS2Test.cpp
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

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( BLAS2Test )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.BLAS2Test" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( geamTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;

    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::geam<ValueType> > geam;

    ContextPtr loc = Context::getContextPtr( geam.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "geam< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    HArray<ValueType> matrix_a( { 1, 2, -3, 4, 5, -6 } );
    HArray<ValueType> matrix_b( { 1, 4, 2, 5, -3, -6 } );

    HArray<ValueType> matrix_c;

    const IndexType m = 2;
    const IndexType n = 3;
    const ValueType alpha = 2;
    const ValueType beta = -1;

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<ValueType> rA( matrix_a, loc );
        ReadAccess<ValueType> rB( matrix_b, loc );
        WriteOnlyAccess<ValueType> wC( matrix_c, loc, m * n );
        geam[loc->getType()]( wC.get(), n, m, n, 
                              alpha, rA.get(), n, common::MatrixOp::NORMAL,
                              beta, rB.get(), n, common::MatrixOp::NORMAL );
    }

    HArray<ValueType> exp_c( { 1, 0, -8, 3, 13, -6 } );

    BOOST_TEST( hostReadAccess( matrix_c ) == hostReadAccess( exp_c ), per_element() );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<ValueType> rA( matrix_a, loc );
        ReadAccess<ValueType> rB( matrix_b, loc );
        WriteOnlyAccess<ValueType> wC( matrix_c, loc, m * n );
        geam[loc->getType()]( wC.get(), n, m, n,
                              ValueType( 1 ), rA.get(), m, common::MatrixOp::TRANSPOSE,
                              ValueType( 0 ), rB.get(), n, common::MatrixOp::NORMAL );
    }
 
    HArray<ValueType> exp_cT( { 1, -3, 5, 2, 4, -6 } );

    BOOST_TEST( hostReadAccess( matrix_c ) == hostReadAccess( exp_cT ), per_element() );

    {
        SCAI_CONTEXT_ACCESS( loc );
        ReadAccess<ValueType> rA( matrix_a, loc );
        ReadAccess<ValueType> rB( matrix_b, loc );
        WriteOnlyAccess<ValueType> wC( matrix_c, loc, m * n );
        geam[loc->getType()]( wC.get(), n, m, n,
                              alpha, rA.get(), n, common::MatrixOp::NORMAL,
                              beta, rB.get(), m, common::MatrixOp::TRANSPOSE );
    }

    // 2 * [ 1 2 -3 ; 4 5 -6 ] - [ 1 2 -3; 4 5 -6 ] -> [1 2 -3; 4 5 -6 ] 

    HArray<ValueType> exp_cNT( { 1, 2, -3, 4, 5, -6 } );

    BOOST_TEST( hostReadAccess( matrix_c ) == hostReadAccess( exp_cNT ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeSquareInPlaceTest, ValueType, blas_test_types )
{
    ContextPtr testContext = Context::getHostPtr();   // in-place not supported by CUDA

    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::geam<ValueType> > geam;

    ContextPtr loc = Context::getContextPtr( geam.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "geam< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    HArray<ValueType> matrix_a( { 1, 2, 3, 4, 5, 6, 7, 8, 9 } );
    HArray<ValueType> matrix_t( { 1, 4, 7, 2, 5, 8, 3, 6, 9 } );

    const IndexType n = 3;

    SCAI_ASSERT_EQ_ERROR( n * n, matrix_a.size(), "serious mismatch" )

    // CUDA does not support in place 

    {
        SCAI_CONTEXT_ACCESS( loc );
        WriteAccess<ValueType> wA( matrix_a, loc );
        geam[loc->getType()]( wA.get(), n, n, n,
                              ValueType( 0 ), wA.get(), n, common::MatrixOp::TRANSPOSE,
                              ValueType( 1 ), wA.get(), n, common::MatrixOp::TRANSPOSE );
    }
 
    BOOST_TEST( hostReadAccess( matrix_a ) == hostReadAccess( matrix_t ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeRectInPlaceTest, ValueType, blas_test_types )
{
    ContextPtr testContext = Context::getHostPtr();   // in-place not supported by CUDA

    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::geam<ValueType> > geam;

    ContextPtr loc = Context::getContextPtr( geam.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "geam< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    /*  Matrix A =  [ 1 -2 3  -4         At =  [  1  5   9
                      5  6 7   8                 -2  6  10
                      9 10 11 12 ]                3  7  11
                                                 -4  8  12 ]
    */

    HArray<ValueType> matrix_a( {  1, -2,  3, -4,  5,   6,   7,  8,   9, 10, 11,  12 } );
    HArray<ValueType> matrix_t( { -1, -5, -9,  2, -6, -10,  -3, -7, -11,  4, -8, -12 } );

    const IndexType m = 3;
    const IndexType n = 4;

    SCAI_ASSERT_EQ_ERROR( m * n, matrix_a.size(), "serious mismatch" )

    {
        SCAI_CONTEXT_ACCESS( loc );
        WriteAccess<ValueType> wA( matrix_a, loc );
        geam[loc->getType()]( wA.get(), n, m, n,
                              ValueType( -1 ), wA.get(), m, common::MatrixOp::TRANSPOSE,
                              ValueType( 0 ), wA.get(), n, common::MatrixOp::NORMAL );
    }
 
    BOOST_TEST( hostReadAccess( matrix_a ) == hostReadAccess( matrix_t ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
    ContextPtr loc = Context::getContextPtr( gemv.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "gemv< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    const common::MatrixOp trans = common::MatrixOp::NORMAL;   // stands for exactly gemv

    ValueType matrix_vals[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };   // row-major

    ValueType x_vals[]      = { 2.0, -1.0, 4.0 };
    ValueType y_vals[]      = { 10.0, -20.0, 30.0 };

    const IndexType n_matrix = sizeof( matrix_vals ) / sizeof( ValueType );
    const IndexType n_x      = sizeof( x_vals ) / sizeof( ValueType );
    const IndexType n_y      = sizeof( y_vals ) / sizeof( ValueType );

    const IndexType m = 2;
    const IndexType n = 3;
    const ValueType alpha = 17.0;
    const IndexType incX = 1;
    const ValueType beta = 13.0;
    const ValueType result_vals[] = { -74, -617, 33 };

    HArray<ValueType> x( n_x, x_vals, testContext );

    HArray<ValueType> matrix( n_matrix , matrix_vals, testContext );

    for ( IndexType incY = 1; incY <= 2; ++incY )
    {
        HArray<ValueType> y( n_y, y_vals, testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rMatrix( matrix, loc );
            ReadAccess<ValueType> rX( x, loc );
            WriteAccess<ValueType> wY( y, loc );
            gemv[loc->getType()]( trans, m, n, alpha, rMatrix.get(), n, rX.get(), incX, beta, wY.get(), incY );
        }
        {
            ReadAccess<ValueType> rY( y );
            BOOST_CHECK_EQUAL( result_vals[0], rY[0] );
            BOOST_CHECK_EQUAL( result_vals[incY], rY[incY] );
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( gevmTest, ValueType, blas_test_types )
{
    const common::MatrixOp trans = common::MatrixOp::TRANSPOSE;   // stands for exactly gemv

    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
    ContextPtr loc = Context::getContextPtr( gemv.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "gemv< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    ValueType matrix_vals[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };

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

    const IndexType   lda  = n;

    HArray<ValueType> matrix( 6, matrix_vals, testContext );

    for ( IndexType incX = 1; incX <= 2; ++incX )
    {
        HArray<ValueType> y( 3, y_vals, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rMatrix( matrix, loc );
            ReadAccess<ValueType> rX( x, loc );
            WriteAccess<ValueType> wY( y, loc );
            gemv[loc->getType()]( trans, m, n, alpha, rMatrix.get(), lda, rX.get(), incX, beta, wY.get(), incY );
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

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( cgevmTest, ValueType, blas_test_types )
{
    // ToDo: do this tes with complex values

    const common::MatrixOp trans = common::MatrixOp::CONJ_TRANSPOSE;

    ContextPtr testContext = ContextFix::testContext;
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
    ContextPtr loc = Context::getContextPtr( gemv.validContext( testContext->getType() ) );
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "gemv< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )

    ValueType matrix_vals[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };  // row-major format

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

    const IndexType   lda         = n;

    HArray<ValueType> matrix( 6, matrix_vals, testContext );

    for ( IndexType incX = 1; incX <= 2; ++incX )
    {
        HArray<ValueType> y( 3, y_vals, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rMatrix( matrix, loc );
            ReadAccess<ValueType> rX( x, loc );
            WriteAccess<ValueType> wY( y, loc );
            gemv[loc->getType()]( trans, m, n, alpha, rMatrix.get(), lda, rX.get(), incX, beta, wY.get(), incY );
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

BOOST_AUTO_TEST_SUITE_END()
