/**
 * @file CUDABLAS3.cpp
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
 * @brief CUDA implementations of BLAS3 routines.
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

// hpp
#include <scai/blaskernel/cuda/CUDABLAS3.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/cuda/CUBLASWrapper.hpp>

// internal scai library
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/tracing.hpp>

namespace scai
{

using namespace tasking;
using common::TypeTraits;

/* ---------------------------------------------------------------------------------------*/
namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( CUDABLAS3::logger, "CUDA.BLAS3" )

/*    gemm                                                                                */
/* ---------------------------------------------------------------------------------------*/

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
    const IndexType ldc )
{
    SCAI_REGION( "CUDA.BLAS3.gemm" )

    typedef CUBLASTrait::BLASTrans BLASTrans;
    BLASTrans transA_char = CUBLAS_OP_N;
    BLASTrans transB_char = CUBLAS_OP_N;
    //Swap matrix if RowMajor Order
    const IndexType lda_call = ( order == CblasRowMajor ) ? ldb : lda;
    const IndexType ldb_call = ( order == CblasRowMajor ) ? lda : ldb;
    const IndexType m_call = ( order == CblasRowMajor ) ? n : m;
    const IndexType n_call = ( order == CblasRowMajor ) ? m : n;
    const ValueType* const A_call = ( order == CblasRowMajor ) ? B : A;
    const ValueType* const B_call = ( order == CblasRowMajor ) ? A : B;
    SCAI_LOG_INFO( logger, "gemm<" << TypeTraits<ValueType>::id() << ">( m = " << m << ", n = " << n << ", k = " << k )

    if ( transa == CblasTrans )
    {
        if ( order == CblasRowMajor )
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

    if ( transb == CblasTrans )
    {
        if ( order == CblasRowMajor )
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

    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0; // default stream if no syncToken is given
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    cublasHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuBLASHandle();
    SCAI_CUBLAS_CALL( cublasSetStream( handle, stream ),
                      "CUDABLAS3::gemm set cublas kernel stream = " << stream );
    SCAI_LOG_INFO( logger, "cublasSgemm: m = " << m_call << " x " << n_call )
    CUBLASWrapper<ValueType>::gemm( handle, transA_char, transB_char,  m_call ,  n_call ,  k , alpha, A_call,  lda_call , B_call,  ldb_call , beta, C,
                                    ldc );

    // No error check here possible as kernel is started asynchronously in any case

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
    }

    SCAI_CUBLAS_CALL( cublasSetStream( handle, NULL ), "CUDABLAS3::gemm set stream to NULL" );
}

/** trsm */

//template<typename ValueType>
//static inline
//void cublasWrapperTrsm(
//    cublasSideMode_t side,
//    cublasFillMode_t uplo,
//    cublasOperation_t transA,
//    cublasDiagType_t diag,
//    IndexType m,
//    IndexType n,
//    ValueType alpha,
//    const ValueType* a,
//    IndexType lda,
//    ValueType* b,
//    IndexType ldb );
//
//template<>
//void cublasWrapperTrsm(
//    cublasSideMode_t side,
//    cublasFillMode_t uplo,
//    cublasOperation_t transA,
//    cublasDiagType_t diag,
//    IndexType m,
//    IndexType n,
//    float alpha,
//    const float* a,
//    IndexType lda,
//    float* b,
//    IndexType ldb )
//{
//    cublasStrsm( CUDAContext_cublasHandle, side, uplo, transA, diag, m, n, &alpha, a, lda, b, ldb );
//}
//
//template<>
//void cublasWrapperTrsm(
//    cublasSideMode_t side,
//    cublasFillMode_t uplo,
//    cublasOperation_t transA,
//    cublasDiagType_t diag,
//    IndexType m,
//    IndexType n,
//    double alpha,
//    const double* a,
//    IndexType lda,
//    double* b,
//    IndexType ldb )
//{
//    cublasDtrsm( CUDAContext_cublasHandle, side, uplo, transA, diag, m, n, &alpha, a, lda, b, ldb );
//}
//
//template<>
//void cublasWrapperTrsm(
//    cublasSideMode_t side,
//    cublasFillMode_t uplo,
//    cublasOperation_t transA,
//    cublasDiagType_t diag,
//    IndexType m,
//    IndexType n,
//    ComplexFloat alpha,
//    const ComplexFloat* a,
//    IndexType lda,
//    ComplexFloat* b,
//    IndexType ldb )
//{
//    cublasCtrsm( CUDAContext_cublasHandle, side, uplo, transA, diag, m, n, cublasCast( &alpha ), cublasCast( a ), lda,
//                 cublasCast( b ), ldb );
//}
//
//template<>
//void cublasWrapperTrsm(
//    cublasSideMode_t side,
//    cublasFillMode_t uplo,
//    cublasOperation_t transA,
//    cublasDiagType_t diag,
//    IndexType m,
//    IndexType n,
//    ComplexDouble alpha,
//    const ComplexDouble* a,
//    IndexType lda,
//    ComplexDouble* b,
//    IndexType ldb )
//{
//    cublasZtrsm( CUDAContext_cublasHandle, side, uplo, transA, diag, m, n, cublasCast( &alpha ), cublasCast( a ), lda,
//                 cublasCast( b ), ldb );
//}
//
//template<typename ValueType>
//void CUDABLAS3::trsm(
//    const CBLAS_ORDER Order,
//    const CBLAS_SIDE sidearg,
//    const CBLAS_UPLO uploarg,
//    const CBLAS_TRANSPOSE trans,
//    const CBLAS_DIAG diagarg,
//    const IndexType m,
//    const IndexType n,
//    const ValueType alpha,
//    const ValueType* A,
//    const IndexType lda,
//    ValueType* B,
//    const IndexType ldb )
//{
//    cublasSideMode_t side = ' ';
//    cublasFillMode_t uplo = ' ';
//    cublasOperation_t transA = ' ';
//    cublasDiagType_t diag = ' ';
//
//    if( trans == CblasTrans )
//    {
////        transA = 'T';
//        transA = CUBLAS_OP_T;
//    }
//    else if( trans == CblasConjTrans )
//    {
////        transA = 'C';
//        transA = CUBLAS_OP_C;
//    }
//    else if( trans == CblasNoTrans )
//    {
////        transA = 'N';
//        transA = CUBLAS_OP_N;
//    }
//    else
//    {
//        BLASHelper::XERBLA_cpu( 0, 4, "cblas_strsm_cuda", "Illegal Trans setting, %d\n", transA );
//        return;
//    }
//
//    if( diagarg == CblasUnit )
//    {
////        diag = 'U';
//        diag = CUBLAS_DIAG_UNIT;
//    }
//    else if( diagarg == CblasNonUnit )
//    {
////        diag = 'N';
//        diag = CUBLAS_DIAG_NON_UNIT;
//    }
//    else
//    {
//        BLASHelper::XERBLA_cpu( 0, 5, "cblas_strsm_cuda", "Illegal Diag setting, %d\n", diagarg );
//        return;
//    }
//
//    if( Order == CblasColMajor )
//    {
//        if( sidearg == CblasRight )
//        {
////            side = 'R';
//            side = CUBLAS_SIDE_RIGHT;
//        }
//        else if( sidearg == CblasLeft )
//        {
////            side = 'L';
//            side = CUBLAS_SIDE_LEFT;
//        }
//        else
//        {
//            BLASHelper::XERBLA_cpu( 0, 2, "cblas_strsm_cuda", "Illegal Side setting, %d\n", sidearg );
//            return;
//        }
//
//        if( uploarg == CblasUpper )
//        {
////            uplo = 'U';
//            uplo = CUBLAS_FILL_MODE_UPPER;
//        }
//        else if( uploarg == CblasLower )
//        {
////            uplo = 'L';
//            uplo = CUBLAS_FILL_MODE_LOWER;
//        }
//        else
//        {
//            BLASHelper::XERBLA_cpu( 0, 3, "cblas_strsm_cuda", "Illegal Uplo setting, %d\n", uploarg );
//            return;
//        }
//    }
//    else if( Order == CblasRowMajor )
//    {
//
//        if( sidearg == CblasRight )
//        {
////            side = 'L';
//            side = CUBLAS_SIDE_LEFT;
//        }
//        else if( sidearg == CblasLeft )
//        {
////            side = 'R';
//            side = CUBLAS_SIDE_RIGHT;
//        }
//        else
//        {
//            BLASHelper::XERBLA_cpu( 1, 2, "cblas_strsm_cuda", "Illegal Side setting, %d\n", sidearg );
//            return;
//        }
//
//        if( uploarg == CblasUpper )
//        {
////            uplo = 'L';
//            uplo = CUBLAS_FILL_MODE_LOWER;
//        }
//        else if( uploarg == CblasLower )
//        {
////            uplo = 'U';
//            uplo = CUBLAS_FILL_MODE_UPPER;
//        }
//        else
//        {
//            BLASHelper::XERBLA_cpu( 1, 3, "cblas_strsm_cuda", "Illegal Uplo setting, %d\n", uploarg );
//            return;
//        }
//    }
//    else
//    {
//        BLASHelper::XERBLA_cpu( 0, 1, "cblas_strsm_cuda", "Illegal order setting, %d\n", Order );
//        return;
//    }
//
//    SCAI_CHECK_CUDA_ACCESS
//
//    cudaStream_t stream = 0; // default stream if no syncToken is given
//
//    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();
//
//    if ( syncToken )
//    {
//        stream = syncToken->getCUDAStream();
//    }
//
//    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ),
//                      "CUDABLAS3::trsm set cublas kernel stream = " << stream );
//
//    cublasWrapperTrsm( side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb );
//
//    // No error check here possible as kernel is started asynchronously in any case
//
//    if( !syncToken )
//    {
//        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
//    }
//
//    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS3::trsm set stream NULL" );
//}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDABLAS3::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_INFO( logger, "register BLAS3 routines implemented by CuBLAS in KernelRegistry [" << flag << "]" )
    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( CUDABLAS3::gemm, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDABLAS3::CUDABLAS3()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

CUDABLAS3::~CUDABLAS3()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

CUDABLAS3 CUDABLAS3::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
