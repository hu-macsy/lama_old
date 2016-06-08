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

    // CblasRowMajor and CblasNoTrans
    {
        ValueType matrix[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType x[] = { 2.0, -1.0, 4.0 };
        ValueType y[] = { 10.0, -20.0, 30.0 };
        const IndexType m = 2;
        const IndexType n = 3;
        const ValueType alpha = 17.0;
        const IndexType lda = 3;
        const IndexType incX = 1;
        const ValueType beta = 13.0;
        const IndexType incY = 2;
        const ValueType result[] = { -74.0, 33.0 };

        HArray<ValueType> Am( 6, matrix, testContext );
        HArray<ValueType> Ax( 3, x, testContext );
        HArray<ValueType> Ay( 3, y, testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAm( Am, loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            gemv[loc->getType()]( CblasRowMajor, CblasNoTrans, m, n, alpha, rAm.get(), lda, rAx.get(), incX, beta, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( result[0], rAy[0] );
            BOOST_CHECK_EQUAL( result[1], rAy[2] );
        }
    }
    // CblasColMajor and CblasNoTrans
    {
        ValueType matrix[] = { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        ValueType x[] = { 2.0, -1.0, 4.0 };
        ValueType y[] = { 10.0, -20.0, 30.0 };

        const IndexType m = 2;
        const IndexType n = 3;
        const ValueType alpha = 17.0;

        const IndexType lda = 2;
        const IndexType incX = 1;
        const ValueType beta = 13.0;
        const IndexType incY = 2;
        const ValueType result[] = { -74.0, 33.0 };

        HArray<ValueType> Am( 6, matrix, testContext );
        HArray<ValueType> Ax( 3, x, testContext );
        HArray<ValueType> Ay( 3, y, testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAm( Am, loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            gemv[loc->getType()]( CblasColMajor, CblasNoTrans, m, n, alpha, rAm.get(), lda, rAx.get(), incX, beta, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( result[0], rAy[0] );
            BOOST_CHECK_EQUAL( result[1], rAy[2] );
        }
    }
    // CblasRowMajor and CblasTrans
    {
        ValueType matrix[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        ValueType x[] =
        { 2.0, -1.0, 4.0 };
        ValueType y[] =
        { 10.0, -20.0, 30.0 };
        const IndexType m = 2;
        const IndexType n = 3;
        const ValueType alpha = 17.0;
        const IndexType lda = 3;
        const IndexType incX = 2;
        const ValueType beta = 13.0;
        const IndexType incY = 1;
        const ValueType result[] =
        { 436.0, 148.0, -120.0 };

        HArray<ValueType> Am( 6, matrix, testContext );
        HArray<ValueType> Ax( 3, x, testContext );
        HArray<ValueType> Ay( 3, y, testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAm( Am, loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            gemv[loc->getType()]( CblasRowMajor, CblasTrans, m, n, alpha, rAm.get(), lda, rAx.get(), incX, beta, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( result[0], rAy[0] );
            BOOST_CHECK_EQUAL( result[1], rAy[1] );
            BOOST_CHECK_EQUAL( result[2], rAy[2] );
        }
    }
    // CblasColMajor and CblasTrans
    {
        ValueType matrix[] = { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        ValueType x[] = { 2.0, -1.0, 4.0 };
        ValueType y[] = { 10.0, -20.0, 30.0 };

        const IndexType m = 2;
        const IndexType n = 3;
        const ValueType alpha = 17.0;
        const IndexType lda = 2;
        const IndexType incX = 2;
        const ValueType beta = 13.0;
        const IndexType incY = 1;
        const ValueType result[] = { 436.0, 148.0, -120.0 };

        HArray<ValueType> Am( 6, matrix, testContext );
        HArray<ValueType> Ax( 3, x, testContext );
        HArray<ValueType> Ay( 3, y, testContext );

        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAm( Am, loc );
            ReadAccess<ValueType> rAx( Ax, loc );
            WriteAccess<ValueType> wAy( Ay, loc );
            gemv[loc->getType()]( CblasColMajor, CblasTrans, m, n, alpha, rAm.get(), lda, rAx.get(), incX, beta, wAy.get(), incY );
        }
        {
            ReadAccess<ValueType> rAy( Ay );
            BOOST_CHECK_EQUAL( result[0], rAy[0] );
            BOOST_CHECK_EQUAL( result[1], rAy[1] );
            BOOST_CHECK_EQUAL( result[2], rAy[2] );
        }
    }
} // gemvTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
