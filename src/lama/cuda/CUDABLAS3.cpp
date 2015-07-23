/**
 * @file CUDABLAS3.cpp
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
 * @brief CUDABLAS3.cpp
 * @author lschubert
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/cuda/CUDABLAS3.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <cudamem/CUDAError.hpp>
#include <cudamem/CUDAStreamSyncToken.hpp>
#include <lama/openmp/BLASHelper.hpp>
#include <lama/cuda/lama_cublas.hpp>

// macros
#include <lama/macros/unused.hpp>

using namespace tasking;
using namespace memory;
using common::getScalarType;

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDABLAS3::logger, "CUDA.BLAS3" )

extern cublasHandle_t CUDAContext_cublasHandle;
/* ---------------------------------------------------------------------------------------*/
/*    gemm                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline
void cublasWrapperGemm(
    cublasOperation_t transA_char,
    cublasOperation_t transB_char,
    IndexType m,
    IndexType n,
    IndexType k,
    ValueType alpha,
    const ValueType* a,
    IndexType lda,
    const ValueType* b,
    IndexType ldb,
    ValueType beta,
    ValueType* c,
    IndexType ldc );

template<>
void cublasWrapperGemm(
    cublasOperation_t transA_char,
    cublasOperation_t transB_char,
    IndexType m,
    IndexType n,
    IndexType k,
    float alpha,
    const float* a,
    IndexType lda,
    const float* b,
    IndexType ldb,
    float beta,
    float* c,
    IndexType ldc )
{
    LAMA_CUBLAS_CALL(
        cublasSgemm( CUDAContext_cublasHandle, transA_char, transB_char, m, n, k, &alpha, a, lda, b, ldb,
                     &beta, c, ldc ),
        "cublasWrapperGemm<float>" );
}

template<>
void cublasWrapperGemm(
    cublasOperation_t transA_char,
    cublasOperation_t transB_char,
    IndexType m,
    IndexType n,
    IndexType k,
    double alpha,
    const double* a,
    IndexType lda,
    const double* b,
    IndexType ldb,
    double beta,
    double* c,
    IndexType ldc )
{
    LAMA_CUBLAS_CALL(
        cublasDgemm( CUDAContext_cublasHandle, transA_char, transB_char, m, n, k, &alpha, a, lda, b, ldb,
                     &beta, c, ldc ),
        "cublasWrapperGemm<dobule>" );
}

template<>
void cublasWrapperGemm(
    cublasOperation_t transA_char,
    cublasOperation_t transB_char,
    IndexType m,
    IndexType n,
    IndexType k,
    ComplexFloat alpha,
    const ComplexFloat* a,
    IndexType lda,
    const ComplexFloat* b,
    IndexType ldb,
    ComplexFloat beta,
    ComplexFloat* c,
    IndexType ldc )
{
    LAMA_CUBLAS_CALL(
        cublasCgemm( CUDAContext_cublasHandle, transA_char, transB_char, m, n, k, cublasCast( &alpha ),
                     cublasCast( a ), lda, cublasCast( b ), ldb, cublasCast( &beta ), cublasCast( c ),
                     ldc ),
        "cublasWrapperGemm<ComplexFloat>" );
}

template<>
void cublasWrapperGemm(
    cublasOperation_t transA_char,
    cublasOperation_t transB_char,
    IndexType m,
    IndexType n,
    IndexType k,
    ComplexDouble alpha,
    const ComplexDouble* a,
    IndexType lda,
    const ComplexDouble* b,
    IndexType ldb,
    ComplexDouble beta,
    ComplexDouble* c,
    IndexType ldc )
{
    LAMA_CUBLAS_CALL(
        cublasZgemm( CUDAContext_cublasHandle, transA_char, transB_char, m, n, k, cublasCast( &alpha ),
                     cublasCast( a ), lda, cublasCast( b ), ldb, cublasCast( &beta ), cublasCast( c ),
                     ldc ),
        "cublasWrapperGemm<double>" );
}

template<typename ValueType>
void CUDABLAS3::gemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transa,
    const CBLAS_TRANSPOSE transb,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const ValueType alpha,
    const ValueType* const A,
    const IndexType lda,
    const ValueType* const B,
    const IndexType ldb,
    const ValueType beta,
    ValueType* const C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    cublasOperation_t transA_char = CUBLAS_OP_N;
    cublasOperation_t transB_char = CUBLAS_OP_N;

    //Swap matrix if RowMajor Order

    const int lda_call = ( order == CblasRowMajor ) ? ldb : lda;
    const int ldb_call = ( order == CblasRowMajor ) ? lda : ldb;
    const int m_call = ( order == CblasRowMajor ) ? n : m;
    const int n_call = ( order == CblasRowMajor ) ? m : n;
    const ValueType* const A_call = ( order == CblasRowMajor ) ? B : A;
    const ValueType* const B_call = ( order == CblasRowMajor ) ? A : B;

    LAMA_LOG_INFO( logger, "gemm<" << getScalarType<ValueType>() << ">( m = " << m << ", n = " << n << ", k = " << k )

    if( transa == CblasTrans )
    {
        if( order == CblasRowMajor )
        {
//          transB_char = 'T';
            transB_char = CUBLAS_OP_T;
        }
        else
        {
//          transA_char = 'T';
            transA_char = CUBLAS_OP_T;
        }
    }

    if( transb == CblasTrans )
    {
        if( order == CblasRowMajor )
        {
//          transA_char = 'T';
            transA_char = CUBLAS_OP_T;
        }
        else
        {
//          transB_char = 'T';
            transB_char = CUBLAS_OP_T;
        }
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ),
                      "CUDABLAS3::gemm set cublas kernel stream = " << stream );

    LAMA_LOG_INFO( logger, "cublasSgemm: m = " << m_call << " x " << n_call )

    cublasWrapperGemm( transA_char, transB_char, m_call, n_call, k, alpha, A_call, lda_call, B_call, ldb_call, beta, C,
                       ldc );

    // No error check here possible as kernel is started asynchronously in any case

    if( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS3::gemm set stream to NULL" );
}

/** trsm */

template<typename ValueType>
static inline
void cublasWrapperTrsm(
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t transA,
    cublasDiagType_t diag,
    IndexType m,
    IndexType n,
    ValueType alpha,
    const ValueType* a,
    IndexType lda,
    ValueType* b,
    IndexType ldb );

template<>
void cublasWrapperTrsm(
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t transA,
    cublasDiagType_t diag,
    IndexType m,
    IndexType n,
    float alpha,
    const float* a,
    IndexType lda,
    float* b,
    IndexType ldb )
{
    cublasStrsm( CUDAContext_cublasHandle, side, uplo, transA, diag, m, n, &alpha, a, lda, b, ldb );
}

template<>
void cublasWrapperTrsm(
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t transA,
    cublasDiagType_t diag,
    IndexType m,
    IndexType n,
    double alpha,
    const double* a,
    IndexType lda,
    double* b,
    IndexType ldb )
{
    cublasDtrsm( CUDAContext_cublasHandle, side, uplo, transA, diag, m, n, &alpha, a, lda, b, ldb );
}

template<>
void cublasWrapperTrsm(
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t transA,
    cublasDiagType_t diag,
    IndexType m,
    IndexType n,
    ComplexFloat alpha,
    const ComplexFloat* a,
    IndexType lda,
    ComplexFloat* b,
    IndexType ldb )
{
    cublasCtrsm( CUDAContext_cublasHandle, side, uplo, transA, diag, m, n, cublasCast( &alpha ), cublasCast( a ), lda,
                 cublasCast( b ), ldb );
}

template<>
void cublasWrapperTrsm(
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t transA,
    cublasDiagType_t diag,
    IndexType m,
    IndexType n,
    ComplexDouble alpha,
    const ComplexDouble* a,
    IndexType lda,
    ComplexDouble* b,
    IndexType ldb )
{
    cublasZtrsm( CUDAContext_cublasHandle, side, uplo, transA, diag, m, n, cublasCast( &alpha ), cublasCast( a ), lda,
                 cublasCast( b ), ldb );
}

template<typename ValueType>
void CUDABLAS3::trsm(
    const CBLAS_ORDER Order,
    const CBLAS_SIDE sidearg,
    const CBLAS_UPLO uploarg,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diagarg,
    const IndexType m,
    const IndexType n,
    const ValueType alpha,
    const ValueType* A,
    const IndexType lda,
    ValueType* B,
    const IndexType ldb,
    SyncToken* syncToken )
{
    IndexType RowMajorStrg = 0;
    cublasSideMode_t side = ' ';
    cublasFillMode_t uplo = ' ';
    cublasOperation_t transA = ' ';
    cublasDiagType_t diag = ' ';

    if( trans == CblasTrans )
    {
//        transA = 'T';
        transA = CUBLAS_OP_T;
    }
    else if( trans == CblasConjTrans )
    {
//        transA = 'C';
        transA = CUBLAS_OP_C;
    }
    else if( trans == CblasNoTrans )
    {
//        transA = 'N';
        transA = CUBLAS_OP_N;
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strsm_cuda", "Illegal Trans setting, %d\n", transA );
        RowMajorStrg = 0;
        return;
    }

    if( diagarg == CblasUnit )
    {
//        diag = 'U';
        diag = CUBLAS_DIAG_UNIT;
    }
    else if( diagarg == CblasNonUnit )
    {
//        diag = 'N';
        diag = CUBLAS_DIAG_NON_UNIT;
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_strsm_cuda", "Illegal Diag setting, %d\n", diagarg );
        RowMajorStrg = 0;
        return;
    }

    if( Order == CblasColMajor )
    {
        if( sidearg == CblasRight )
        {
//            side = 'R';
            side = CUBLAS_SIDE_RIGHT;
        }
        else if( sidearg == CblasLeft )
        {
//            side = 'L';
            side = CUBLAS_SIDE_LEFT;
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strsm_cuda", "Illegal Side setting, %d\n", sidearg );
            RowMajorStrg = 0;
            return;
        }

        if( uploarg == CblasUpper )
        {
//            uplo = 'U';
            uplo = CUBLAS_FILL_MODE_UPPER;
        }
        else if( uploarg == CblasLower )
        {
//            uplo = 'L';
            uplo = CUBLAS_FILL_MODE_LOWER;
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strsm_cuda", "Illegal Uplo setting, %d\n", uploarg );
            RowMajorStrg = 0;
            return;
        }
    }
    else if( Order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if( sidearg == CblasRight )
        {
//            side = 'L';
            side = CUBLAS_SIDE_LEFT;
        }
        else if( sidearg == CblasLeft )
        {
//            side = 'R';
            side = CUBLAS_SIDE_RIGHT;
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strsm_cuda", "Illegal Side setting, %d\n", sidearg );
            RowMajorStrg = 0;
            return;
        }

        if( uploarg == CblasUpper )
        {
//            uplo = 'L';
            uplo = CUBLAS_FILL_MODE_LOWER;
        }
        else if( uploarg == CblasLower )
        {
//            uplo = 'U';
            uplo = CUBLAS_FILL_MODE_UPPER;
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

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ),
                      "CUDABLAS3::trsm set cublas kernel stream = " << stream );

    cublasWrapperTrsm( side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb );

    // No error check here possible as kernel is started asynchronously in any case

    if( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS3::trsm set stream NULL" );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDABLAS3::setInterface( BLASInterface& BLAS )
{
    LAMA_LOG_INFO( logger, "set BLAS3 routines for CUDA in Interface" )

    // Note: macro takes advantage of same name for routines and type definitions
    //       ( e.g. routine CUDABLAS1::sum<ValueType> is set for BLAS::BLAS1::sum variable

#define LAMA_BLAS3_REGISTER(z, I, _)                                            \
    LAMA_INTERFACE_REGISTER_T( BLAS, gemm, ARITHMETIC_TYPE##I )                 \

    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_BLAS3_REGISTER, _ )

#undef LAMA_BLAS3_REGISTER

    // trsm routines are not used yet by LAMA
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDABLAS3::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( context::CUDA );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDABLAS3::initialized = registerInterface();

} /* namespace lama */
