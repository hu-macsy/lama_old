/**
 * @file CUDABLAS1.cpp
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
 * @brief Wrapper implementations for BLAS1 routines in CUDA using cuBLAS
 * @author Lauretta Schubert, Thomas Brandes, Eric Stricker
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/blaskernel/cuda/CUDABLAS1.hpp>

// local library
#include <scai/blaskernel/cuda/CUBLASTrait.hpp>
#include <scai/blaskernel/cuda/CUBLASWrapper.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/Container.hpp>

using namespace scai::tasking;
using namespace scai::hmemo;
using scai::common::TypeTraits;

namespace scai
{

extern cublasHandle_t CUDAContext_cublasHandle;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( CUDABLAS1::logger, "CUDA.BLAS1" )

/* ---------------------------------------------------------------------------------------*/
/*    sum                                                                                 */
/* ---------------------------------------------------------------------------------------*/

template<typename T>
__global__
void sum_kernel( const int n, T alpha, const T* x, T beta, const T* y, T* z )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        z[i] = alpha * x[i] + beta * y[i];
    }
}

/* ---------------------------------------------------------------------------------------*/
/*    sum                                                                                 */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void CUDABLAS1::sum(
    const IndexType n,
    ValueType alpha,
    const ValueType* x,
    ValueType beta,
    const ValueType* y,
    ValueType* z )
{
    SCAI_REGION( "CUDA.BLAS1.sum" )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger,
                    "sum<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", " << alpha << " * x + " << beta << " * y " )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    sum_kernel <<< dimGrid, dimBlock, 0, stream>>> ( n, alpha, x, beta, y, z );

    if( !syncToken )
    {
        cudaStreamSynchronize( stream );
        SCAI_CHECK_CUDA_ERROR
    }
}

/* ---------------------------------------------------------------------------------------*/
/*    scale                                                                               */
/* ---------------------------------------------------------------------------------------*/

// Note: the cublasWrapper routines could be static routines on its own. But using
//       a common template routine is helpful to guarantee correct syntax

template<typename ValueType>
void CUDABLAS1::scal( IndexType n, const ValueType alpha, ValueType* x_d, const IndexType incX )
{
    SCAI_REGION( "CUDA.BLAS1.scal" )

	typedef CUBLASTrait::BLASIndexType BLASIndexType;

    if( incX == 0 )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger, "scal<" << TypeTraits<ValueType>::id() << "> of x[" << n << "], alpha = " << alpha )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    SyncToken* syncToken = SyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::scal set stream" );

    CUBLASWrapper<ValueType>::scal(  static_cast<BLASIndexType>(n), alpha, x_d,  static_cast<BLASIndexType>(incX)  );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::scal set stream" );
}

/* ---------------------------------------------------------------------------------------*/
/*    nrm2                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType CUDABLAS1::nrm2( IndexType n, const ValueType* x_d, IndexType incX )
{
    SCAI_REGION( "CUDA.BLAS1.nrm2" )

	typedef CUBLASTrait::BLASIndexType BLASIndexType;

    if( incX <= 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    SCAI_LOG_DEBUG( logger, "nrm2<" << TypeTraits<ValueType>::id() << "> of x[" << n << "]" )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    // Note: we have to switch cublas Stream, this might be done globally later

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::nrm2 set stream" );

    ValueType res = CUBLASWrapper<ValueType>::nrm2(  static_cast<BLASIndexType>(n) , x_d,  static_cast<BLASIndexType>(incX)  );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::nrm2 set stream null" );

    return res;
}

/* ---------------------------------------------------------------------------------------*/
/*    asum                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType CUDABLAS1::asum( const IndexType n, const ValueType* x_d, const IndexType incX )
{
    SCAI_REGION( "CUDA.BLAS1.asum" )

	typedef CUBLASTrait::BLASIndexType BLASIndexType;

    if( incX <= 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    SCAI_LOG_DEBUG( logger, "asum<" << TypeTraits<ValueType>::id() << "> of x[" << n << "]" )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::asum set stream" );

    ValueType res = CUBLASWrapper<ValueType>::asum(  static_cast<BLASIndexType>(n) , x_d,  static_cast<BLASIndexType>(incX)  );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::asum set stream NULL" );
    return res;
}

/* ---------------------------------------------------------------------------------------*/
/*    iamax                                                                               */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType CUDABLAS1::iamax( const IndexType n, const ValueType* x_d, const IndexType incX )
{
    SCAI_REGION( "CUDA.BLAS1.iamax" )

	typedef CUBLASTrait::BLASIndexType BLASIndexType;

    SCAI_LOG_DEBUG( logger, "iamax<" << TypeTraits<ValueType>::id() << "> of x[" << n << "]" )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUABLAS1::iamax set stream" );

    IndexType iamax = CUBLASWrapper<ValueType>::iamax(  static_cast<BLASIndexType>(n) , x_d,  static_cast<BLASIndexType>(incX) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::iamax set stream NULL" );
    return iamax ? iamax - 1 : 0;
}

/* ---------------------------------------------------------------------------------------*/
/*    swap                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void CUDABLAS1::swap(
    const IndexType n,
    ValueType* x_d,
    const IndexType incX,
    ValueType* y_d,
    const IndexType incY )
{
    SCAI_REGION( "CUDA.BLAS1.swap" )

	typedef CUBLASTrait::BLASIndexType BLASIndexType;

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger, "swap<" << TypeTraits<ValueType>::id() << "> of x, y with size " << n )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS::swap set stream" );

    CUBLASWrapper<ValueType>::swap(  static_cast<BLASIndexType>(n) , x_d,  static_cast<BLASIndexType>(incX) , y_d,  static_cast<BLASIndexType>(incY)  );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUADABLAS1::swap set stream NULL" );
}

/* ---------------------------------------------------------------------------------------*/
/*    copy                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void CUDABLAS1::copy(
    IndexType n,
    const ValueType* x_d,
    IndexType incX,
    ValueType* y_d,
    IndexType incY )
{
    SCAI_REGION( "CUDA.BLAS1.copy" )

	typedef CUBLASTrait::BLASIndexType BLASIndexType;

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger, "copy<" << TypeTraits<ValueType>::id() << "> of x, y, n = " << n )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::copy set stream" );

    CUBLASWrapper<ValueType>::copy(  static_cast<BLASIndexType>(n) , x_d,  static_cast<BLASIndexType>(incX) , y_d,  static_cast<BLASIndexType>(incY) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::copy set stream NULL" );
}

/* ---------------------------------------------------------------------------------------*/
/*    axpy                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void CUDABLAS1::axpy(
    int n,
    ValueType alpha,
    const ValueType* x_d,
    int incX,
    ValueType* y_d,
    const int incY )
{
    SCAI_REGION( "CUDA.BLAS1.axpy" )

	typedef CUBLASTrait::BLASIndexType BLASIndexType;

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger, "axpy<" << TypeTraits<ValueType>::id() << "> of x, y, n = " << n << ", alpha = " << alpha )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::axpy set stream" );

    CUBLASWrapper<ValueType>::axpy(  static_cast<BLASIndexType>(n) , alpha, x_d,  static_cast<BLASIndexType>(incX) , y_d,  static_cast<BLASIndexType>(incY) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::axpy set stream NULL" );
}

/* ---------------------------------------------------------------------------------------*/
/*    dot                                                                                 */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType CUDABLAS1::dot(
    IndexType n,
    const ValueType* x_d,
    IndexType incX,
    const ValueType* y_d,
    IndexType incY )
{
    SCAI_REGION( "CUDA.BLAS1.dot" )

	typedef CUBLASTrait::BLASIndexType BLASIndexType;

    SCAI_LOG_DEBUG( logger,
                    "dot<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", incX = " << incX << ", incY = " << incY << ", x_d = " << x_d << ", y_d = " << y_d )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return static_cast<ValueType>(0.0);
    }

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::dot set stream" );

    ValueType res = CUBLASWrapper<ValueType>::dot(  static_cast<BLASIndexType>(n) , x_d,  static_cast<BLASIndexType>(incX) , y_d,
    		static_cast<BLASIndexType>(incY)  );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR
    }

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::dot set stream NULL" );
    return res;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDABLAS1::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::CUDA;

    SCAI_LOG_INFO( logger, "register BLAS1 routines implemented by CuBLAS in KernelRegistry [" << flag << "]" )

    KernelRegistry::set<BLASKernelTrait::sum<ValueType> >( CUDABLAS1::sum, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::scal<ValueType> >( CUDABLAS1::scal, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::nrm2<ValueType> >( CUDABLAS1::nrm2, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::asum<ValueType> >( CUDABLAS1::asum, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::iamax<ValueType> >( CUDABLAS1::iamax, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::swap<ValueType> >( CUDABLAS1::swap, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::copy<ValueType> >( CUDABLAS1::copy, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::axpy<ValueType> >( CUDABLAS1::axpy, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::dot<ValueType> >( CUDABLAS1::dot, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDABLAS1::CUDABLAS1()
{
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_CUDA_LIST>::call(
                    kregistry::KernelRegistry::KERNEL_ADD );
}

CUDABLAS1::~CUDABLAS1()
{
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_CUDA_LIST>::call(
                    kregistry::KernelRegistry::KERNEL_ERASE );
}

CUDABLAS1 CUDABLAS1::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
