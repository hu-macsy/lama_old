/**
 * @file BLAS3Test.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains tests for the blas3 methods.
 * @author: Bea Hornef
 * @date 17.7.2013
 * @since 1.0.0
 **/

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/hmemo.hpp>
#include <scai/lama/LAMAKernel.hpp>
#include <scai/lama/BLASKernelTrait.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

namespace scai
{
namespace lama
{
namespace BLAS3Test
{

template<typename ValueType>
void gemmTest( ContextPtr loc )
{
    //  input
    //                            (  2.0 3.0 )
    // 17.0 * ( 1.0  2.0 -3.0 ) * ( -1.0 1.0 ) - 13.0 * ( 15.0 13.0 ) =  (-9.0 -1.0)
    //        ( 4.0  5.0 -6.0 )   (  4.0 5.0 )          ( 27.0 17.0 )    (-6.0  0.0)

    LAMAKernel<BLASKernelTrait::gemm<ValueType> > gemm;

    const ValueType alpha = 17.0;
    const ValueType beta = 13.0;
    const ValueType resultRowMajor[] =
    { -9.0, -1.0, -6.0, 0.0 };
    const ValueType resultColMajor[] =
    { -9.0, -6.0, -1.0, 0.0 };
    const IndexType n = 2;
    const IndexType m = 2;
    const IndexType k = 3;
    // CblasRowMajor and 2 x CblasNoTrans
    {
        const ValueType matrixA[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const ValueType matrixB[] =
        { 2.0, 3.0, -1.0, 1.0, 4.0, 5.0 };
        ValueType matrixC[] =
        { 15.0, 13.0, 27.0, 17.0 };
        const IndexType lda = 3;
        const IndexType ldb = 2;
        const IndexType ldc = 2;
        LAMAArray<ValueType> AmA( 6, matrixA );
        LAMAArray<ValueType> AmB( 6, matrixB );
        LAMAArray<ValueType> AmC( 4, matrixC );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc]( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
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
    // CblasColMajor and 2 x CblasNoTrans
    {
        const ValueType matrixA[] =
        { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        const ValueType matrixB[] =
        { 2.0, -1.0, 4.0, 3.0, 1.0, 5.0 };
        ValueType matrixC[] =
        { 15.0, 27.0, 13.0, 17.0 };
        const IndexType lda = 2;
        const IndexType ldb = 3;
        const IndexType ldc = 2;
        LAMAArray<ValueType> AmA( 6, matrixA );
        LAMAArray<ValueType> AmB( 6, matrixB );
        LAMAArray<ValueType> AmC( 4, matrixC );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc]( CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
                       wAmC.get(), ldc );
        }
        {
            ReadAccess<ValueType> rAmC( AmC );

            for ( int i = 0; i < 4; ++i )
            {
                BOOST_CHECK_EQUAL( resultColMajor[i], rAmC[i] );
            }
        }
    }
    // CblasRowMajor, CblasNoTrans for A and CblasTrans for B
    {
        const ValueType matrixA[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const ValueType matrixB[] =
        { 2.0, -1.0, 4.0, 3.0, 1.0, 5.0 };
        ValueType matrixC[] =
        { 15.0, 13.0, 27.0, 17.0 };
        const IndexType lda = 3;
        const IndexType ldb = 3;
        const IndexType ldc = 2;
        LAMAArray<ValueType> AmA( 6, matrixA );
        LAMAArray<ValueType> AmB( 6, matrixB );
        LAMAArray<ValueType> AmC( 4, matrixC );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc]( CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
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
    // CblasColMajor, CblasNoTrans for A and CblasTrans for B
    {
        const ValueType matrixA[] =
        { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        const ValueType matrixB[] =
        { 2.0, 3.0, -1.0, 1.0, 4.0, 5.0 };
        ValueType matrixC[] =
        { 15.0, 27.0, 13.0, 17.0 };
        const IndexType lda = 2;
        const IndexType ldb = 2;
        const IndexType ldc = 2;
        LAMAArray<ValueType> AmA( 6, matrixA );
        LAMAArray<ValueType> AmB( 6, matrixB );
        LAMAArray<ValueType> AmC( 4, matrixC );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc]( CblasColMajor, CblasNoTrans, CblasTrans, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
                  wAmC.get(), ldc );
        }
        {
            ReadAccess<ValueType> rAmC( AmC );

            for ( int i = 0; i < 4; ++i )
            {
                BOOST_CHECK_EQUAL( resultColMajor[i], rAmC[i] );
            }
        }
    }
    // CblasRowMajor, CblasTrans for A and CblasNoTrans for B
    {
        const ValueType matrixA[] =
        { 1.0, 4.0, 2.0, 5.0, -3.0, -6.0 };
        const ValueType matrixB[] =
        { 2.0, 3.0, -1.0, 1.0, 4.0, 5.0 };
        ValueType matrixC[] =
        { 15.0, 13.0, 27.0, 17.0 };
        const IndexType lda = 2;
        const IndexType ldb = 2;
        const IndexType ldc = 2;
        LAMAArray<ValueType> AmA( 6, matrixA );
        LAMAArray<ValueType> AmB( 6, matrixB );
        LAMAArray<ValueType> AmC( 4, matrixC );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc]( CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
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
    // CblasColMajor, CblasTrans for A and CblasNoTrans for B
    {
        const ValueType matrixA[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const ValueType matrixB[] =
        { 2.0, -1.0, 4.0, 3.0, 1.0, 5.0 };
        ValueType matrixC[] =
        { 15.0, 27.0, 13.0, 17.0 };
        const IndexType lda = 3;
        const IndexType ldb = 3;
        const IndexType ldc = 2;
        LAMAArray<ValueType> AmA( 6, matrixA );
        LAMAArray<ValueType> AmB( 6, matrixB );
        LAMAArray<ValueType> AmC( 4, matrixC );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc]( CblasColMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
                       wAmC.get(), ldc );
        }
        {
            ReadAccess<ValueType> rAmC( AmC );

            for ( int i = 0; i < 4; ++i )
            {
                BOOST_CHECK_EQUAL( resultColMajor[i], rAmC[i] );
            }
        }
    }
    // CblasRowMajor, CblasTrans for A and CblasTrans for B
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
        LAMAArray<ValueType> AmA( 6, matrixA );
        LAMAArray<ValueType> AmB( 6, matrixB );
        LAMAArray<ValueType> AmC( 4, matrixC );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc]( CblasRowMajor, CblasTrans, CblasTrans, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
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
    // CblasColMajor, CblasTrans for A and CblasTrans for B
    {
        const ValueType matrixA[] =
        { 1.0, 2.0, -3.0, 4.0, 5.0, -6.0 };
        const ValueType matrixB[] =
        { 2.0, 3.0, -1.0, 1.0, 4.0, 5.0 };
        ValueType matrixC[] =
        { 15.0, 27.0, 13.0, 17.0 };
        const IndexType lda = 3;
        const IndexType ldb = 2;
        const IndexType ldc = 2;
        LAMAArray<ValueType> AmA( 6, matrixA );
        LAMAArray<ValueType> AmB( 6, matrixB );
        LAMAArray<ValueType> AmC( 4, matrixC );
        {
            SCAI_CONTEXT_ACCESS( loc );
            ReadAccess<ValueType> rAmA( AmA, loc );
            ReadAccess<ValueType> rAmB( AmB, loc );
            WriteAccess<ValueType> wAmC( AmC, loc );
            gemm[loc]( CblasColMajor, CblasTrans, CblasTrans, m, n, k, alpha, rAmA.get(), lda, rAmB.get(), ldb, beta,
                       wAmC.get(), ldc );
        }
        {
            ReadAccess<ValueType> rAmC( AmC );

            for ( int i = 0; i < 4; ++i )
            {
                BOOST_CHECK_EQUAL( resultColMajor[i], rAmC[i] );
            }
        }
    }
} // gemmTest

} // namespace BLAS3Test
} /* end namespace lama */

} /* end namespace scai */
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE ( BLAS3Test )

SCAI_LOG_DEF_LOGGER( logger, "Test.BLAS3Test" )

LAMA_AUTO_TEST_CASE_CT( gemmTest, BLAS3Test, scai::lama )
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
