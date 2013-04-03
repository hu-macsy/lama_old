/**
 * @file CUDABLAS2.cpp
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
 * @brief CUDABLAS2.cpp
 * @author lschubert
 * @date 05.07.2012
 * $Id$
 */

// hpp
#include <lama/cuda/CUDABLAS2.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>

// macros
#include <lama/macros/unused.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDABLAS2::logger, "CUDA.BLAS2" )

//TODO: throw exception instead of setlasterror

/** gemv */

template<>
void CUDABLAS2::gemv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE trans,
    const IndexType m,
    const IndexType n,
    const float alpha,
    const float* const A,
    const IndexType lda,
    const float* const x,
    const IndexType incx,
    const float beta,
    float* const y,
    const IndexType incy,
    SyncToken* syncToken )
{
    IndexType order_m = m;
    IndexType order_n = n;
    char trans_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }

        order_m = n;
        order_n = m;
    }
    else
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    LAMA_LOG_INFO( logger, "cublasSgemv: m = " << order_m << " x " << order_n )

    cublasSgemv( trans_char, order_m, order_n, alpha, A, lda, x, incx, beta, y, incy );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::gemv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE trans,
    const IndexType m,
    const IndexType n,
    const double alpha,
    const double* const A,
    const IndexType lda,
    const double* const x,
    const IndexType incx,
    const double beta,
    double* const y,
    const IndexType incy,
    SyncToken* syncToken )
{
    IndexType order_m = m;
    IndexType order_n = n;
    char trans_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }

        order_m = n;
        order_n = m;
    }
    else
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    LAMA_LOG_INFO( logger, "cublasDgemv: m = " << order_m << " x " << order_n )

    cublasDgemv( trans_char, order_m, order_n, alpha, A, lda, x, incx, beta, y, incy );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** symv */

template<>
void CUDABLAS2::symv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* x,
    const IndexType incX,
    const float beta,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasSsymv( uplo_char, n, alpha, A, lda, x, incX, beta, y, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::symv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* x,
    const IndexType incX,
    const double beta,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDsymv( uplo_char, n, alpha, A, lda, x, incX, beta, y, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** trmv */

template<>
void CUDABLAS2::trmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const float* A,
    const IndexType lda,
    float* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasStrmv( uplo_char, trans_char, diag_char, n, A, lda, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::trmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const double* A,
    const IndexType lda,
    double* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDtrmv( uplo_char, trans_char, diag_char, n, A, lda, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** trsv */

template<>
void CUDABLAS2::trsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const float* A,
    const IndexType lda,
    float* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasStrsv( uplo_char, trans_char, diag_char, n, A, lda, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::trsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const double* A,
    const IndexType lda,
    double* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDtrsv( uplo_char, trans_char, diag_char, n, A, lda, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** gbmv */

template<>
void CUDABLAS2::gbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE trans,
    const IndexType m,
    const IndexType n,
    IndexType kl,
    IndexType ku,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* x,
    const IndexType incX,
    const float beta,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    IndexType order_m = m;
    IndexType order_n = n;
    IndexType order_kl = kl;
    IndexType order_ku = ku;
    char trans_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }

        order_m = n;
        order_n = m;
        order_kl = ku;
        order_ku = kl;
    }
    else
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasSgbmv( trans_char, order_m, order_n, order_kl, order_ku, alpha, A, lda, x, incX, beta, y, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::gbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE trans,
    const IndexType m,
    const IndexType n,
    IndexType kl,
    IndexType ku,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* x,
    const IndexType incX,
    const double beta,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    IndexType order_m = m;
    IndexType order_n = n;
    IndexType order_kl = kl;
    IndexType order_ku = ku;
    char trans_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }

        order_m = n;
        order_n = m;
        order_kl = ku;
        order_ku = kl;
    }
    else
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDgbmv( trans_char, order_m, order_n, order_kl, order_ku, alpha, A, lda, x, incX, beta, y, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** sbmv */

template<>
void CUDABLAS2::sbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const IndexType k,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* x,
    const IndexType incX,
    const float beta,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasSsbmv( uplo_char, n, k, alpha, A, lda, x, incX, beta, y, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::sbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const IndexType k,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* x,
    const IndexType incX,
    const double beta,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDsbmv( uplo_char, n, k, alpha, A, lda, x, incX, beta, y, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** tbmv */

template<>
void CUDABLAS2::tbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const IndexType k,
    const float* A,
    const IndexType lda,
    float* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasStbmv( uplo_char, trans_char, diag_char, n, k, A, lda, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::tbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const IndexType k,
    const double* A,
    const IndexType lda,
    double* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDtbmv( uplo_char, trans_char, diag_char, n, k, A, lda, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** tbsv */

template<>
void CUDABLAS2::tbsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const IndexType k,
    const float* A,
    const IndexType lda,
    float* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasStbsv( uplo_char, trans_char, diag_char, n, k, A, lda, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::tbsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const IndexType k,
    const double* A,
    const IndexType lda,
    double* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDtbsv( uplo_char, trans_char, diag_char, n, k, A, lda, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** ger */

template<>
void CUDABLAS2::ger(
    const enum CBLAS_ORDER order,
    const IndexType m,
    const IndexType n,
    const float alpha,
    const float* x,
    const IndexType incX,
    const float* y,
    const IndexType incY,
    float* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    if ( order == CblasRowMajor )
    {
        cublasSger( n, m, alpha, y, incY, x, incX, A, lda );
        LAMA_CHECK_CUBLAS_ERROR
        ;
    }
    else
    {
        cublasSger( m, n, alpha, x, incX, y, incY, A, lda );
        LAMA_CHECK_CUBLAS_ERROR
        ;
    }

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::ger(
    const enum CBLAS_ORDER order,
    const IndexType m,
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    const double* y,
    const IndexType incY,
    double* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    if ( order == CblasRowMajor )
    {
        cublasDger( n, m, alpha, y, incY, x, incX, A, lda );
        LAMA_CHECK_CUBLAS_ERROR
        ;
    }
    else
    {
        cublasDger( m, n, alpha, x, incX, y, incY, A, lda );
        LAMA_CHECK_CUBLAS_ERROR
        ;
    }

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** syr */

template<>
void CUDABLAS2::syr(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const float alpha,
    const float* x,
    const IndexType incX,
    float* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasSsyr( uplo_char, n, alpha, x, incX, A, lda );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::syr(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    double* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDsyr( uplo_char, n, alpha, x, incX, A, lda );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** syr2 */

template<>
void CUDABLAS2::syr2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const float alpha,
    const float* x,
    const IndexType incX,
    const float* y,
    const IndexType incY,
    float* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasSsyr2( uplo_char, n, alpha, x, incX, y, incY, A, lda );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::syr2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    const double* y,
    const IndexType incY,
    double* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDsyr2( uplo_char, n, alpha, x, incX, y, incY, A, lda );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** spmv */

template<>
void CUDABLAS2::spmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const float alpha,
    const float* AP,
    const float* x,
    const IndexType incX,
    const float beta,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasSspmv( uplo_char, n, alpha, AP, x, incX, beta, y, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::spmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const double alpha,
    const double* AP,
    const double* x,
    const IndexType incX,
    const double beta,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor

    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDspmv( uplo_char, n, alpha, AP, x, incX, beta, y, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** spr */

template<>
void CUDABLAS2::spr(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const float alpha,
    const float* x,
    const IndexType incX,
    float* AP,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasSspr( uplo_char, n, alpha, x, incX, AP );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::spr(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    double* AP,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDspr( uplo_char, n, alpha, x, incX, AP );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** spr2 */

template<>
void CUDABLAS2::spr2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const float alpha,
    const float* x,
    const IndexType incX,
    const float* y,
    const IndexType incY,
    float* AP,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasSspr2( uplo_char, n, alpha, x, incX, y, incY, AP );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::spr2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    const double* y,
    const IndexType incY,
    double* AP,
    SyncToken* syncToken )
{
    char uplo_char = ' ';

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'U';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else if ( uplo == CblasLower )
        {
            uplo_char = 'L';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDspr2( uplo_char, n, alpha, x, incX, y, incY, AP );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** tpmv */

template<>
void CUDABLAS2::tpmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const float* AP,
    float* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasStpmv( uplo_char, trans_char, diag_char, n, AP, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::tpmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const double* AP,
    double* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDtpmv( uplo_char, trans_char, diag_char, n, AP, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/** tpsv */

template<>
void CUDABLAS2::tpsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const float* AP,
    float* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasStpsv( uplo_char, trans_char, diag_char, n, AP, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

template<>
void CUDABLAS2::tpsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const IndexType n,
    const double* AP,
    double* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    char uplo_char = ' ';
    char trans_char = ' ';
    char diag_char = ' ';

    if ( diag == CblasUnit )
    {
        diag_char = 'U';
    }
    else if ( diag == CblasNonUnit )
    {
        diag_char = 'N';
    }

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'L';
        }
        else
        {
            uplo_char = 'U';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'N';
        }
        else
        {
            trans_char = 'T';
        }
    }
    else
    {
        if ( uplo == CblasUpper )
        {
            uplo_char = 'U';
        }
        else
        {
            uplo_char = 'L';
        }

        if ( trans == CblasTrans )
        {
            trans_char = 'T';
        }
        else
        {
            trans_char = 'N';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDtpsv( uplo_char, trans_char, diag_char, n, AP, x, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDABLAS2::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, gemv, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, gemv, double )

    // other routines are not used by LAMA yet
}

} /* namespace lama */
