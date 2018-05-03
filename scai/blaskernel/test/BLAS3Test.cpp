/**
 * @file BLAS3Test.cpp
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
 * @brief Contains tests for the blas3 methods.
 * @author Bea Hornef
 * @date 17.07.2013
 */

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/hmemo.hpp>
#include <scai/kregistry/KernelContextFunction.hpp>

#include <scai/blaskernel/test/TestMacros.hpp>

using namespace scai;
using namespace scai::hmemo;
using common::TypeTraits;
using common::MatrixOp;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( BLAS3Test )

/* --------------------------------------------------------------------- */

// use of Fixture ContextFix provides the testContext

BOOST_GLOBAL_FIXTURE( ContextFix );

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.BLAS3Test" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemmTest, ValueType, blas_test_types )
{
    ContextPtr testContext = ContextFix::testContext;
    //  input
    //                            (  2.0 3.0 )
    // 17.0 * ( 1.0  2.0 -3.0 ) * ( -1.0 1.0 ) - 13.0 * ( 15.0 13.0 ) =  (-9.0 -1.0)
    //        ( 4.0  5.0 -6.0 )   (  4.0 5.0 )          ( 27.0 17.0 )    (-6.0  0.0)
    kregistry::KernelTraitContextFunction<blaskernel::BLASKernelTrait::gemm<ValueType> > gemm;
    ContextPtr loc = Context::getContextPtr( gemm.validContext( testContext->getType() ) );
    // give warning if routine is not available at testContext
    BOOST_WARN_EQUAL( loc->getType(), testContext->getType() );
    SCAI_LOG_INFO( logger, "gemm< " << TypeTraits<ValueType>::id() << "> test for " << *testContext << " on " << *loc )
    const ValueType alpha = 17.0;
    const ValueType beta = 13.0;
    const ValueType resultRowMajor[] =
    { -9.0, -1.0, -6.0, 0.0 };
    const IndexType n = 2;
    const IndexType m = 2;
    const IndexType k = 3;
    // 2 x MatrixOp::NORMAL
    {
        const ValueType matrixA[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const ValueType matrixB[] = { 2.0, 3.0, -1.0, 1.0, 4.0, 5.0 };
        ValueType matrixC[]       = { 15.0, 13.0, 27.0, 17.0 };
        const IndexType lda = 3;
        const IndexType ldb = 2;
        const IndexType ldc = 2;
        const IndexType nA = sizeof( matrixA ) / sizeof( ValueType );
        const IndexType nB = sizeof( matrixB ) / sizeof( ValueType );
        const IndexType nC = sizeof( matrixC ) / sizeof( ValueType );
        HArray<ValueType> AmA( nA, matrixA, testContext );
        HArray<ValueType> AmB( nB, matrixB, testContext );
        HArray<ValueType> AmC( nC, matrixC, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc->getType()]( MatrixOp::NORMAL, MatrixOp::NORMAL, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
                                  wAmC.get(), ldc );
        }
        {
            ReadAccess<ValueType> rAmC( AmC );

            for ( int i = 0; i < 4; ++i )
            {
                BOOST_CHECK_EQUAL( resultRowMajor[i], rAmC[i] );
            }
        }
    }
    // MatrixOp::NORMAL for A and MatrixOp::TRANSPOSE for B
    {
        const ValueType matrixA[] = { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const ValueType matrixB[] = { 2.0, -1.0, 4.0, 3.0, 1.0, 5.0 };
        ValueType matrixC[]       = { 15.0, 13.0, 27.0, 17.0 };
        const IndexType lda = 3;
        const IndexType ldb = 3;
        const IndexType ldc = 2;
        const IndexType nA = sizeof( matrixA ) / sizeof( ValueType );
        const IndexType nB = sizeof( matrixB ) / sizeof( ValueType );
        const IndexType nC = sizeof( matrixC ) / sizeof( ValueType );
        HArray<ValueType> AmA( nA, matrixA, testContext );
        HArray<ValueType> AmB( nB, matrixB, testContext );
        HArray<ValueType> AmC( nC, matrixC, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc->getType()]( MatrixOp::NORMAL, MatrixOp::TRANSPOSE, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
                                  wAmC.get(), ldc );
        }
        {
            ReadAccess<ValueType> rAmC( AmC );

            for ( int i = 0; i < 4; ++i )
            {
                BOOST_CHECK_EQUAL( resultRowMajor[i], rAmC[i] );
            }
        }
    }

    // MatrixOp::TRANSPOSE for A and MatrixOp::NORMAL for B
    {
        const ValueType matrixA[] = { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        const ValueType matrixB[] = { 2.0, 3.0, -1.0, 1.0, 4.0, 5.0 };
        ValueType matrixC[]       = { 15.0, 13.0, 27.0, 17.0 };
        const IndexType lda = 2;
        const IndexType ldb = 2;
        const IndexType ldc = 2;
        HArray<ValueType> AmA( 6, matrixA, testContext );
        HArray<ValueType> AmB( 6, matrixB, testContext );
        HArray<ValueType> AmC( 4, matrixC, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc->getType()]( MatrixOp::TRANSPOSE, MatrixOp::NORMAL, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
                                  wAmC.get(), ldc );
        }
        {
            ReadAccess<ValueType> rAmC( AmC );

            for ( int i = 0; i < 4; ++i )
            {
                BOOST_CHECK_EQUAL( resultRowMajor[i], rAmC[i] );
            }
        }
    }
    // MatrixOp::TRANSPOSE for A and MatrixOp::TRANSPOSE for B
    {
        const ValueType matrixA[] =
        { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        const ValueType matrixB[] =
        { 2.0, -1.0, 4.0, 3.0, 1.0, 5.0 };
        ValueType matrixC[] =
        { 15.0, 13.0, 27.0, 17.0 };
        const IndexType lda = 2;
        const IndexType ldb = 3;
        const IndexType ldc = 2;
        const IndexType nA = sizeof( matrixA ) / sizeof( ValueType );
        const IndexType nB = sizeof( matrixB ) / sizeof( ValueType );
        const IndexType nC = sizeof( matrixC ) / sizeof( ValueType );
        HArray<ValueType> AmA( nA, matrixA, testContext );
        HArray<ValueType> AmB( nB, matrixB, testContext );
        HArray<ValueType> AmC( nC, matrixC, testContext );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc->getType()]( MatrixOp::TRANSPOSE, MatrixOp::TRANSPOSE, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
                                  wAmC.get(), ldc );
        }
        {
            ReadAccess<ValueType> rAmC( AmC );

            for ( int i = 0; i < 4; ++i )
            {
                BOOST_CHECK_EQUAL( resultRowMajor[i], rAmC[i] );
            }
        }
    }

} // gemmTest

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
