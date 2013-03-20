/**
 * @file CUDALAPACK.cpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief CUDALAPACK.cpp
 * @author lschubert
 * @date 06.07.2012
 * $Id$
 */

#include <lama/cuda/CUDABLAS1.hpp>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDALAPACK.hpp>

namespace lama
{

//TODO implement
//    template<>
//    T CUDALAPACK::lamch( const enum CBLAS_MACHINE_PARAM param );
//
//    template<>
//    IndexType CUDALAPACK::getrf(
//        const enum CBLAS_ORDER order,
//        const IndexType m,
//        const IndexType n,
//        float* const A,
//        const IndexType lda,
//        IndexType* const ipiv )
//
//    template<>
//    IndexType CUDALAPACK::getrf(
//        const enum CBLAS_ORDER order,
//        const IndexType m,
//        const IndexType n,
//        double* const A,
//        const IndexType lda,
//        IndexType* const ipiv )
//
//    template<>
//    IndexType CUDALAPACK::getrf(
//        const enum CBLAS_ORDER order,
//        const IndexType m,
//        const IndexType n,
//        float* const A,
//        const IndexType lda,
//        IndexType* const ipiv)
//
//    template<>
//    IndexType CUDALAPACK::getrf(
//        const enum CBLAS_ORDER order,
//        const IndexType m,
//        const IndexType n,
//        double* const A,
//        const IndexType lda,
//        IndexType* const ipiv)
//
//    template<>
//    IndexType CUDALAPACK::getri(
//        const enum CBLAS_ORDER order,
//        const IndexType n,
//        float* const A,
//        const IndexType lda,
//        IndexType* const ipiv)
//
//    template<>
//    IndexType CUDALAPACK::getri(
//        const enum CBLAS_ORDER order,
//        const IndexType n,
//        double* const A,
//        const IndexType lda,
//        IndexType* const ipiv)
//
//    template<>
//    IndexType CUDALAPACK::trtrs(
//        const enum CBLAS_ORDER order,
//        const enum CBLAS_UPLO uplo,
//        const enum CBLAS_TRANSPOSE trans,
//        const enum CBLAS_DIAG diag,
//        const IndexType n,
//        const IndexType nrhs,
//        const flaot* A,
//        const IndexType lda,
//        float* B,
//        const IndexType ldb)
//
//    template<>
//    IndexType CUDALAPACK::trtrs(
//        const enum CBLAS_ORDER order,
//        const enum CBLAS_UPLO uplo,
//        const enum CBLAS_TRANSPOSE trans,
//        const enum CBLAS_DIAG diag,
//        const IndexType n,
//        const IndexType nrhs,
//        const double* A,
//        const IndexType lda,
//        double* B,
//        const IndexType ldb)
//
//    template<>
//    IndexType CUDALAPACK::tptrs(
//        const enum CBLAS_ORDER order,
//        const enum CBLAS_UPLO uplo,
//        const enum CBLAS_TRANSPOSE trans,
//        const enum CBLAS_DIAG diag,
//        const IndexType n,
//        const IndexType nrhs,
//        const float* AP,
//        float* B,
//        const IndexType ldb)
//
//    template<>
//    IndexType CUDALAPACK::tptrs(
//        const enum CBLAS_ORDER order,
//        const enum CBLAS_UPLO uplo,
//        const enum CBLAS_TRANSPOSE trans,
//        const enum CBLAS_DIAG diag,
//        const IndexType n,
//        const IndexType nrhs,
//        const double* AP,
//        double* B,
//        const IndexType ldb)

//TODO look at info and feedback variables
template<>
void CUDALAPACK::laswp(
    const enum CBLAS_ORDER order,
    const IndexType n,
    float* A_d,
    const IndexType lda,
    const IndexType k1,
    const IndexType k2,
    const IndexType* ipiv_h,
    const IndexType incx,
    SyncToken* syncToken )
{
    int info = 0;
    int i = k1;

    if ( order == CblasRowMajor )
    {
        int feedback = 0;

        for ( i = k1; i < k2 /*&& feedback == LAMA_STATUS_SUCCESS*/; ++i )
        {
            if ( ipiv_h[i * incx] == i )
            {
                continue;
            }

            CUDABLAS1::swap( n, &A_d[ipiv_h[i * incx] * lda], incx, &A_d[i * lda], incx, syncToken );
            LAMA_CHECK_CUDA_ERROR
            ;
        }

        info = -1 * (IndexType) feedback;
    }
    else if ( order == CblasColMajor )
    {
        info = n + lda;
    }
    else
    {
        info = 1;
        BLASHelper::XERBLA_cpu( 0, info, "slaswp", "Illegal order setting." );
    }

    if ( info < 0 )
    {
        //TODO: throw exception
    }
//        return info;
}

template<>
void CUDALAPACK::laswp(
    const enum CBLAS_ORDER order,
    const IndexType n,
    double* A_d,
    const IndexType lda,
    const IndexType k1,
    const IndexType k2,
    const IndexType* ipiv_h,
    const IndexType incx,
    SyncToken* syncToken )
{
    int info = 0;
    int i = k1;

    if ( order == CblasRowMajor )
    {
        int feedback = 0;

        for ( i = k1; i < k2 /*&& feedback == LAMA_STATUS_SUCCESS*/; ++i )
        {
            if ( ipiv_h[i * incx] == i )
            {
                continue;
            }

            CUDABLAS1::swap( n, &A_d[ipiv_h[i * incx] * lda], incx, &A_d[i * lda], incx, syncToken );
            LAMA_CHECK_CUDA_ERROR
            ;
        }

        info = -1 * (IndexType) feedback;
    }
    else if ( order == CblasColMajor )
    {
        info = n + lda;
    }
    else
    {
        info = 1;
        BLASHelper::XERBLA_cpu( 0, info, "dlaswp", "Illegal order setting." );
    }

    if ( info < 0 )
    {
        //TODO: throw exception
    }
//        return info;
}

} /* namespace lama */

