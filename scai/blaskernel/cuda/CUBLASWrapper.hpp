/**
 * @file blaskernel/cuda/CUBLASWrapper.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Wrapper class for calls of the cuBLAS library to resolve overloading.
 * @author eschricker
 * @date 24.08.2015
 */

#pragma once

// internal scai libraries
#include <scai/blaskernel/cuda/CUBLASTrait.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/TypeTraits.hpp>

// CUDA
#ifdef SCAI_COMPLEX_SUPPORTED
#include <cuComplex.h>
#endif

#include <cublas_v2.h>

namespace scai
{

namespace blaskernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT CUBLASWrapper;

#define CUBLASWRAPPER_DEF( ValueType, CUBLASValueType, prefix1, prefix2, prefix3, DOT )                                     \
    template<>                                                                                                              \
    class COMMON_DLL_IMPORTEXPORT CUBLASWrapper<ValueType>                                                                  \
    {                                                                                                                       \
    public:                                                                                                                 \
        \
        typedef CUBLASTrait::BLASIndexType BLASIndexType;                                                                   \
        typedef CUBLASTrait::BLASTrans BLASTrans;                                                                           \
        \
        static void scal( cublasHandle_t handle, const IndexType n, const ValueType alpha,                                  \
                          ValueType* x, const IndexType incX )                                                              \
        {                                                                                                                   \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            const CUBLASValueType* b_alpha = reinterpret_cast<const CUBLASValueType*>( &alpha );                            \
            CUBLASValueType* b_x = reinterpret_cast<CUBLASValueType*>( x );                                                 \
            \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( scal, prefix1 )( handle, b_n, b_alpha, b_x, b_incX),                        \
                              "CUBLASWrapper::scal<" #ValueType ">()" );                                                    \
        }                                                                                                                   \
        \
        static ValueType nrm2( cublasHandle_t handle, const IndexType n, const ValueType *x, const IndexType incX )         \
        {                                                                                                                   \
            common::TypeTraits<ValueType>::RealType nrm2 = 0;                                                                \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( nrm2, prefix2 )( handle,                                                    \
                              b_n, reinterpret_cast<const CUBLASValueType*>(x), b_incX,                                     \
                              reinterpret_cast<common::TypeTraits<ValueType>::RealType*>(&nrm2)),                            \
                              "CUBLASWrapper::nrm2<" #ValueType ">");                                                       \
            return nrm2;                                                                                                    \
        }                                                                                                                   \
        \
        static ValueType asum( cublasHandle_t handle, const IndexType n, const ValueType *x, IndexType incX)                \
        {                                                                                                                   \
            common::TypeTraits<ValueType>::RealType asum = 0;                                                                \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( asum, prefix2 )( handle, b_n,                                               \
                              reinterpret_cast<const CUBLASValueType*>(x), b_incX,                                          \
                              reinterpret_cast<common::TypeTraits<ValueType>::RealType*>( &asum )),                          \
                              "CUBLASWrapper::asum<" #ValueType ">");                                                       \
            return asum;                                                                                                    \
        }                                                                                                                   \
        \
        static BLASIndexType iamax( cublasHandle_t handle, const IndexType n, const ValueType *x, const IndexType incX)     \
        {                                                                                                                   \
            BLASIndexType iamax;                                                                                            \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            SCAI_CUBLAS_CALL(CUBLAS_BLAS_NAME( amax, prefix3 )( handle,                                                     \
                             b_n, reinterpret_cast<const CUBLASValueType*>(x), b_incX, &iamax),                             \
                             "CUBLASWrapper::iamax<" #ValueType ">");                                                       \
            return iamax;                                                                                                   \
        }                                                                                                                   \
        \
        static void swap( cublasHandle_t handle, const BLASIndexType n, ValueType *x, const IndexType incX, ValueType *y,   \
                          const IndexType incY)                                                                             \
        {                                                                                                                   \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            BLASIndexType b_incY = static_cast<BLASIndexType>( incY );                                                      \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( swap, prefix1 )( handle, b_n,                                               \
                              reinterpret_cast<CUBLASValueType*>(x), b_incX,                                                \
                              reinterpret_cast<CUBLASValueType*>(y), b_incY ),                                              \
                              "CUBLASWrapper::swap<" #ValueType ">");                                                       \
        }                                                                                                                   \
        \
        static void copy( cublasHandle_t handle, const IndexType n, const ValueType *x, const IndexType incX,               \
                          ValueType *y, const IndexType incY)                                                               \
        {                                                                                                                   \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            BLASIndexType b_incY = static_cast<BLASIndexType>( incY );                                                      \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( copy, prefix1 )( handle, b_n,                                               \
                              reinterpret_cast<const CUBLASValueType*>(x), b_incX,                                          \
                              reinterpret_cast<CUBLASValueType*>(y), b_incY ),                                              \
                              "CUBLASWrapper::copy<" #ValueType ">");                                                       \
        }                                                                                                                   \
        \
        static void axpy( cublasHandle_t handle, const IndexType n, const ValueType alpha,                                  \
                          const ValueType *x, const IndexType incX, ValueType *y, const IndexType incY)                     \
        {                                                                                                                   \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            BLASIndexType b_incY = static_cast<BLASIndexType>( incY );                                                      \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( axpy, prefix1 )( handle, b_n,                                               \
                              reinterpret_cast<const CUBLASValueType*>(&alpha),                                             \
                              reinterpret_cast<const CUBLASValueType*>(x),                                                  \
                              b_incX, reinterpret_cast<CUBLASValueType*>(y), b_incY),                                       \
                              "CUBLASWrapper::axpy<" #ValueType ">");                                                       \
        }                                                                                                                   \
        \
        static ValueType dot( cublasHandle_t handle, const IndexType n, const ValueType *x,                                 \
                              const IndexType incX, const ValueType *y, const IndexType incY)                               \
        {                                                                                                                   \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            BLASIndexType b_incY = static_cast<BLASIndexType>( incY );                                                      \
            ValueType dot;                                                                                                  \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( DOT, prefix1 )( handle, b_n,                                                \
                              reinterpret_cast<const CUBLASValueType*>(x), b_incX,                                          \
                              reinterpret_cast<const CUBLASValueType*>(y), b_incY,                                          \
                              reinterpret_cast<CUBLASValueType*>(&dot)),                                                    \
                              "CUBLASWrapper::dot<" #ValueType ">");                                                        \
            return dot;                                                                                                     \
        }                                                                                                                   \
        \
        static void gemv(  cublasHandle_t handle, const BLASTrans transA, const IndexType m, const IndexType n,             \
                           const ValueType alpha, const ValueType* A, const IndexType lda, const ValueType* x,              \
                           const IndexType incX, const ValueType beta, ValueType* y, const IndexType incY)                  \
        {                                                                                                                   \
            BLASIndexType b_m    = static_cast<BLASIndexType>( m );                                                         \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_lda  = static_cast<BLASIndexType>( lda );                                                       \
            BLASIndexType b_incX = static_cast<BLASIndexType>( incX );                                                      \
            BLASIndexType b_incY = static_cast<BLASIndexType>( incY );                                                      \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( gemv, prefix1 )( handle, transA,                                            \
                              b_m, b_n, reinterpret_cast<const CUBLASValueType*>(&alpha),                                   \
                              reinterpret_cast<const CUBLASValueType*>(A),                                                  \
                              b_lda, reinterpret_cast<const CUBLASValueType*>(x), b_incX,                                   \
                              reinterpret_cast<const CUBLASValueType*>(&beta), reinterpret_cast<CUBLASValueType*>(y),       \
                              b_incY), "CUBLASWrapper::gemv<" #ValueType ">");                                              \
        }                                                                                                                   \
        static void geam(  cublasHandle_t handle, const cublasOperation_t transA, const cublasOperation_t transB,           \
                           const IndexType m, const IndexType n,                                                            \
                           const ValueType alpha, const ValueType* A, const IndexType lda,                                  \
                           const ValueType beta, const ValueType* B, const IndexType ldb,                                   \
                           ValueType* C, const IndexType ldc )                                                              \
        {                                                                                                                   \
            BLASIndexType b_m    = static_cast<BLASIndexType>( m );                                                         \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_lda  = static_cast<BLASIndexType>( lda );                                                       \
            BLASIndexType b_ldb  = static_cast<BLASIndexType>( ldb );                                                       \
            BLASIndexType b_ldc  = static_cast<BLASIndexType>( ldc );                                                       \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( geam, prefix1 )( handle,                                                    \
                              transA, transB,                                                                               \
                              b_m, b_n,                                                                                     \
                              reinterpret_cast<const CUBLASValueType*>( &alpha ),                                           \
                              reinterpret_cast<const CUBLASValueType*>( A ), b_lda,                                         \
                              reinterpret_cast<const CUBLASValueType*>( &beta ),                                            \
                              reinterpret_cast<const CUBLASValueType*>( B ), b_ldb,                                         \
                              reinterpret_cast<CUBLASValueType*>( C ), b_ldc ),                                             \
                              "CUBLASWrapper::geam<" #ValueType ">");                                                       \
        }                                                                                                                   \
        \
        static void gemm(  cublasHandle_t handle, const BLASTrans transA, const BLASTrans transB,                           \
                           const IndexType m, const IndexType n, const IndexType k,                                         \
                           const ValueType alpha, const ValueType* A, const IndexType lda, const ValueType* B,              \
                           const IndexType ldb, const ValueType beta, ValueType* C, const IndexType ldc)                    \
        {                                                                                                                   \
            BLASIndexType b_m    = static_cast<BLASIndexType>( m );                                                         \
            BLASIndexType b_k    = static_cast<BLASIndexType>( k );                                                         \
            BLASIndexType b_n    = static_cast<BLASIndexType>( n );                                                         \
            BLASIndexType b_lda  = static_cast<BLASIndexType>( lda );                                                       \
            BLASIndexType b_ldb  = static_cast<BLASIndexType>( ldb );                                                       \
            BLASIndexType b_ldc  = static_cast<BLASIndexType>( ldc );                                                       \
            SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( gemm, prefix1 )( handle, transA, transB,                                    \
                              b_m, b_n, b_k,                                                                                \
                              reinterpret_cast<const CUBLASValueType*>( &alpha ),                                           \
                              reinterpret_cast<const CUBLASValueType*>( A ), b_lda,                                         \
                              reinterpret_cast<const CUBLASValueType*>( B ), b_ldb,                                         \
                              reinterpret_cast<const CUBLASValueType*>( &beta ),                                            \
                              reinterpret_cast<CUBLASValueType*>( C ), b_ldc ),                                             \
                              "CUBLASWrapper::gemm<" #ValueType ">");                                                       \
        }                                                                                                                   \
    };

CUBLASWRAPPER_DEF( float, float, S, S, Is, dot )
CUBLASWRAPPER_DEF( double, double, D, D, Id, dot )

#ifdef SCAI_COMPLEX_SUPPORTED
CUBLASWRAPPER_DEF( ComplexFloat, cuFloatComplex, C, Sc, Ic, dotc )
CUBLASWRAPPER_DEF( ComplexDouble, cuDoubleComplex, Z, Dz, Iz, dotc )
#endif

#undef CUBLASWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */

