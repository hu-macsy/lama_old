/**
 * @file CUDABLAS3.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief CUDABLAS3.cpp
 * @author lschubert
 * @date 05.07.2012
 * $Id$
 */

// hpp
#include <lama/cuda/CUDABLAS3.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>

// macros
#include <lama/macros/unused.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDABLAS3::logger, "CUDA.BLAS3" )

/** gemm */

template<>
void CUDABLAS3::gemm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE transa,
    const enum CBLAS_TRANSPOSE transb,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const float alpha,
    const float* const A,
    const IndexType lda,
    const float* const B,
    const IndexType ldb,
    const float beta,
    float* const C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    char transA_char = 'N';
    char transB_char = 'N';

    //Swap matrix if RowMajor Order

    const int lda_call = ( order == CblasRowMajor ) ? ldb : lda;
    const int ldb_call = ( order == CblasRowMajor ) ? lda : ldb;
    const int m_call = ( order == CblasRowMajor ) ? n : m;
    const int n_call = ( order == CblasRowMajor ) ? m : n;
    const float* const A_call = ( order == CblasRowMajor ) ? B : A;
    const float* const B_call = ( order == CblasRowMajor ) ? A : B;

    LAMA_LOG_INFO( logger, "gemm<float>( m = " << m << ", n = " << n << ", k = " << k )

    if ( transa == CblasTrans )
    {
        transA_char = 'T';
    }

    if ( transb == CblasTrans )
    {
        transB_char = 'T';
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

    LAMA_LOG_INFO( logger, "cublasSgemm: m = " << m_call << " x " << n_call )

    cublasSgemm( transA_char, transB_char, m_call, n_call, k, alpha, A_call, lda_call, B_call, ldb_call, beta, C, ldc );

    // No error check here possible as kernel is started asynchronously in any case

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
    }
}

template<>
void CUDABLAS3::gemm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE transa,
    const enum CBLAS_TRANSPOSE transb,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const double alpha,
    const double* const A,
    const IndexType lda,
    const double* const B,
    const IndexType ldb,
    const double beta,
    double* const C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    char transA_char = 'N';
    char transB_char = 'N';

    //Swap matrix if RowMajor Order

    const int lda_call = ( order == CblasRowMajor ) ? ldb : lda;
    const int ldb_call = ( order == CblasRowMajor ) ? lda : ldb;
    const int m_call = ( order == CblasRowMajor ) ? n : m;
    const int n_call = ( order == CblasRowMajor ) ? m : n;
    const double* const A_call = ( order == CblasRowMajor ) ? B : A;
    const double* const B_call = ( order == CblasRowMajor ) ? A : B;

    if ( transa == CblasTrans )
    {
        transA_char = 'T';
    }

    if ( transb == CblasTrans )
    {
        transB_char = 'T';
    }

    LAMA_LOG_INFO( logger, "gemm<double>( m = " << m << ", n = " << n << ", k = " << k )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    LAMA_LOG_INFO( logger, "cublasDgemm: m = " << m_call << " x " << n_call )

    cublasDgemm( transA_char, transB_char, m_call, n_call, k, alpha, A_call, lda_call, B_call, ldb_call, beta, C, ldc );

    // No error check here possible as kernel is started asynchronously in any case

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
    }
}

/** symm */

/** trmm */

/** trsm */

template<>
void CUDABLAS3::trsm(
    const enum CBLAS_ORDER Order,
    const enum CBLAS_SIDE sidearg,
    const enum CBLAS_UPLO uploarg,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diagarg,
    const IndexType m,
    const IndexType n,
    const float alpha,
    const float* A,
    const IndexType lda,
    float* B,
    const IndexType ldb,
    SyncToken* syncToken )
{
    IndexType RowMajorStrg = 0;
    char side = ' ';
    char uplo = ' ';
    char transA = ' ';
    char diag = ' ';

    if ( trans == CblasTrans )
    {
        transA = 'T';
    }
    else if ( trans == CblasConjTrans )
    {
        transA = 'C';
    }
    else if ( trans == CblasNoTrans )
    {
        transA = 'N';
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strsm_cuda", "Illegal Trans setting, %d\n", transA );
        RowMajorStrg = 0;
        return;
    }

    if ( diagarg == CblasUnit )
    {
        diag = 'U';
    }
    else if ( diagarg == CblasNonUnit )
    {
        diag = 'N';
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_strsm_cuda", "Illegal Diag setting, %d\n", diagarg );
        RowMajorStrg = 0;
        return;
    }

    if ( Order == CblasColMajor )
    {
        if ( sidearg == CblasRight )
        {
            side = 'R';
        }
        else if ( sidearg == CblasLeft )
        {
            side = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strsm_cuda", "Illegal Side setting, %d\n", sidearg );
            RowMajorStrg = 0;
            return;
        }

        if ( uploarg == CblasUpper )
        {
            uplo = 'U';
        }
        else if ( uploarg == CblasLower )
        {
            uplo = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strsm_cuda", "Illegal Uplo setting, %d\n", uploarg );
            RowMajorStrg = 0;
            return;
        }
    }
    else if ( Order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( sidearg == CblasRight )
        {
            side = 'L';
        }
        else if ( sidearg == CblasLeft )
        {
            side = 'R';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strsm_cuda", "Illegal Side setting, %d\n", sidearg );
            RowMajorStrg = 0;
            return;
        }

        if ( uploarg == CblasUpper )
        {
            uplo = 'L';
        }
        else if ( uploarg == CblasLower )
        {
            uplo = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strsm_cuda", "Illegal Uplo setting, %d\n", uploarg );
            RowMajorStrg = 0;
            return;
        }
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_strsm_cuda", "Illegal order setting, %d\n", Order );
        RowMajorStrg = 0;
        return;
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

    cublasStrsm( side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb );

    // No error check here possible as kernel is started asynchronously in any case

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
    }
}

template<>
void CUDABLAS3::trsm(
    const enum CBLAS_ORDER Order,
    const enum CBLAS_SIDE sidearg,
    const enum CBLAS_UPLO uploarg,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diagarg,
    const IndexType m,
    const IndexType n,
    const double alpha,
    const double* A,
    const IndexType lda,
    double* B,
    const IndexType ldb,
    SyncToken* syncToken )
{
    char side = ' ';
    char uplo = ' ';
    char transA = ' ';
    char diag = ' ';

    if ( trans == CblasTrans )
    {
        transA = 'T';
    }
    else if ( trans == CblasConjTrans )
    {
        transA = 'C';
    }
    else if ( trans == CblasNoTrans )
    {
        transA = 'N';
    }

    if ( diagarg == CblasUnit )
    {
        diag = 'U';
    }
    else if ( diagarg == CblasNonUnit )
    {
        diag = 'N';
    }

    if ( Order == CblasColMajor )
    {
        if ( sidearg == CblasRight )
        {
            side = 'R';
        }
        else if ( sidearg == CblasLeft )
        {
            side = 'L';
        }

        if ( uploarg == CblasUpper )
        {
            uplo = 'U';
        }
        else if ( uploarg == CblasLower )
        {
            uplo = 'L';
        }
    }
    else if ( Order == CblasRowMajor )
    {
        if ( sidearg == CblasRight )
        {
            side = 'L';
        }
        else if ( sidearg == CblasLeft )
        {
            side = 'R';
        }

        if ( uploarg == CblasUpper )
        {
            uplo = 'L';
        }
        else if ( uploarg == CblasLower )
        {
            uplo = 'U';
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" );
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetKernelStream( stream ), "set cublas kernel stream = " << stream );

    cublasDtrsm( side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb );

    // No error check here possible as kernel is started asynchronously in any case

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDABLAS3::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, gemm, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, gemm, double )

    // trsm routines are not used yet by LAMA
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDABLAS3::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDABLAS3::initialized = registerInterface();

} /* namespace lama */
