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
#include <lama/cuda/CUDABLAS1.hpp>

// others
#include <cudamem/CUDAError.hpp>
#include <cudamem/CUDAStreamSyncToken.hpp>
#include <lama/cuda/lama_cublas.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// macros
#include <lama/macros/unused.hpp>

// tracing with LAMA_REGION
#include <tracing/tracing.hpp>

// blas
#include <boost/preprocessor.hpp>

using namespace tasking;
using namespace memory;
using common::getScalarType;

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDABLAS1::logger, "CUDA.BLAS1" )

extern cublasHandle_t CUDAContext_cublasHandle;

/* ---------------------------------------------------------------------------------------*/
/*    scale                                                                               */
/* ---------------------------------------------------------------------------------------*/

// Note: the cublasWrapper routines could be static routines on its own. But using
//       a common template routine is helpful to guarantee correct syntax
template<typename ValueType>
static inline void cublasWrapperScale( int n, ValueType alpha, ValueType* x_d, int incX );

template<>
void cublasWrapperScale( int n, float alpha, float* x_d, int incX )
{
    LAMA_CUBLAS_CALL( cublasSscal( CUDAContext_cublasHandle, n, &alpha, x_d, incX ), "cublasWrapperScale<float>" );
}

template<>
void cublasWrapperScale( int n, double alpha, double* x_d, int incX )
{
    LAMA_CUBLAS_CALL( cublasDscal( CUDAContext_cublasHandle, n, &alpha, x_d, incX ), "cublasWrapperScale<double>" );
}

template<>
void cublasWrapperScale( int n, ComplexFloat alpha, ComplexFloat* x_d, int incX )
{
    // use of cublasCast to convert ComplexFloat to cuComplex via reinterpret_cast
    LAMA_CUBLAS_CALL( cublasCscal( CUDAContext_cublasHandle, n, cublasCast( &alpha ), cublasCast( x_d ), incX ),
                      "cublasWrapperScale<ComplexFloat>" );
}

template<>
void cublasWrapperScale( int n, ComplexDouble alpha, ComplexDouble* x_d, int incX )
{
    // use of cublasCast to convert ComplexDouble to cuDoubleComplex via reinterpret_cast
    LAMA_CUBLAS_CALL( cublasZscal( CUDAContext_cublasHandle, n, cublasCast( &alpha ), cublasCast( x_d ), incX ),
                      "cublasWrapperScale<ComplexDouble>" );
}

template<typename ValueType>
void CUDABLAS1::scal( IndexType n, const ValueType alpha, ValueType* x_d, const IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.scal" )

    if( incX == 0 )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "scal<" << getScalarType<ValueType>() << "> of x[" << n << "], alpha = " << alpha )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::scal set stream" );

    cublasWrapperScale( static_cast<int>( n ), alpha, x_d, static_cast<int>( incX ) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::scal set stream" );
}

/* ---------------------------------------------------------------------------------------*/
/*    nrm2                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline ValueType cublasWrapperNrm2( int n, const ValueType* x_d, int incX );

template<>
float cublasWrapperNrm2( int n, const float* x_d, int incX )
{
    float nrm2;
    LAMA_CUBLAS_CALL( cublasSnrm2( CUDAContext_cublasHandle, n, x_d, incX, &nrm2 ), "cublasWrapperNrm2<float>" );
    return nrm2;
}

template<>
double cublasWrapperNrm2( int n, const double* x_d, int incX )
{
    double nrm2;
    LAMA_CUBLAS_CALL( cublasDnrm2( CUDAContext_cublasHandle, n, x_d, incX, &nrm2 ), "cublasWrapperNrm2<double>" );
    return nrm2;
}

template<>
ComplexFloat cublasWrapperNrm2( int n, const ComplexFloat* x_d, int incX )
{
    // CUBLAS returns only float result so we convert it back to Complex
    float nrm2;
    LAMA_CUBLAS_CALL( cublasScnrm2( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, &nrm2 ),
                      "cublasWrapperNrm2<ComplexFloat>" );
    return ComplexFloat( nrm2, 0.0f );
}

template<>
ComplexDouble cublasWrapperNrm2( int n, const ComplexDouble* x_d, int incX )
{
    // CUBLAS returns only double result so we convert it back to Complex
    double nrm2;
    LAMA_CUBLAS_CALL( cublasDznrm2( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, &nrm2 ),
                      "cublasWrapperNrm2<ComplexDouble>" );
    return ComplexDouble( nrm2, 0.0 );
}

template<typename ValueType>
ValueType CUDABLAS1::nrm2( IndexType n, const ValueType* x_d, IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.nrm2" )

    if( incX <= 0 )
    {
        return 0.0;
    }

    LAMA_LOG_DEBUG( logger, "nrm2<" << getScalarType<ValueType>() << "> of x[" << n << "]" )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::nrm2 set stream" );

    ValueType res = cublasWrapperNrm2( static_cast<int>( n ), x_d, static_cast<int>( incX ) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::nrm2 set stream null" );
    return res;
}

/* ---------------------------------------------------------------------------------------*/
/*    asum                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline ValueType cublasWrapperAsum( int n, const ValueType* x_d, int incX );

template<>
float cublasWrapperAsum( int n, const float* x_d, int incX )
{
    float asum;
    LAMA_CUBLAS_CALL( cublasSasum( CUDAContext_cublasHandle, n, x_d, incX, &asum ), "cublasWrapperAsum<float>" );
    return asum;
}

template<>
double cublasWrapperAsum( int n, const double* x_d, int incX )
{
    double asum;
    LAMA_CUBLAS_CALL( cublasDasum( CUDAContext_cublasHandle, n, x_d, incX, &asum ), "cublasWrapperAsum<double>" );
    return asum;
}

template<>
ComplexFloat cublasWrapperAsum( int n, const ComplexFloat* x_d, int incX )
{
    // CUBLAS returns only float result so we convert it back to Complex
    float asum;
    LAMA_CUBLAS_CALL( cublasScasum( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, &asum ),
                      "cublasWrapperAsum<ComplexFloat>" );
    return ComplexFloat( asum, 0.0f );
}

template<>
ComplexDouble cublasWrapperAsum( int n, const ComplexDouble* x_d, int incX )
{
    // CUBLAS returns only double result so we convert it back to Complex
    double asum;
    LAMA_CUBLAS_CALL( cublasDzasum( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, &asum ),
                      "cublasWrapperAsum<ComplexDouble>" );
    return ComplexDouble( asum, 0.0 );
}

template<typename ValueType>
ValueType CUDABLAS1::asum( const IndexType n, const ValueType* x_d, const IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.asum" )

    if( incX <= 0 )
    {
        return 0.0;
    }

    LAMA_LOG_DEBUG( logger, "asum<" << getScalarType<ValueType>() << "> of x[" << n << "]" )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::asum set stream" );

    ValueType res = cublasWrapperAsum( static_cast<int>( n ), x_d, static_cast<int>( incX ) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::asum set stream NULL" );
    return res;
}

/* ---------------------------------------------------------------------------------------*/
/*    iamax                                                                               */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static int cublasWrapperIamax( int n, const ValueType* x_d, int incX );

template<>
int cublasWrapperIamax( int n, const float* x_d, int incX )
{
    int iamax;
    LAMA_CUBLAS_CALL( cublasIsamax( CUDAContext_cublasHandle, n, x_d, incX, &iamax ), "cublasWrapperIamax<float>" );
    return iamax;
}

template<>
int cublasWrapperIamax( int n, const double* x_d, int incX )
{
    int iamax;
    LAMA_CUBLAS_CALL( cublasIdamax( CUDAContext_cublasHandle, n, x_d, incX, &iamax ), "cublasWrapperIamax<double>" );
    return iamax;
}

template<>
int cublasWrapperIamax( int n, const ComplexFloat* x_d, int incX )
{
    int iamax;
    LAMA_CUBLAS_CALL( cublasIcamax( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, &iamax ),
                      "cublasWrapperIamax<ComplexFloat>" );
    return iamax;
}

template<>
int cublasWrapperIamax( int n, const ComplexDouble* x_d, int incX )
{
    int iamax;
    LAMA_CUBLAS_CALL( cublasIzamax( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, &iamax ),
                      "cublasWrapperIamax<ComplexDouble>" );
    return iamax;
}

template<typename ValueType>
IndexType CUDABLAS1::iamax( const IndexType n, const ValueType* x_d, const IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.iamax" )

    LAMA_LOG_DEBUG( logger, "iamax<" << getScalarType<ValueType>() << "> of x[" << n << "]" )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUABLAS1::iamax set stream" );

    IndexType iamax = cublasWrapperIamax( n, x_d, incX );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::iamax set stream NULL" );
    return iamax ? iamax - 1 : 0;
}

/* ---------------------------------------------------------------------------------------*/
/*    swap                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline void cublasWrapperSwap( IndexType n, ValueType* x_d, IndexType incX, ValueType* y_d, IndexType incY );

template<>
void cublasWrapperSwap( int n, float* x_d, int incX, float* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasSswap( CUDAContext_cublasHandle, n, x_d, incX, y_d, incY ), "cublasWrapperSwap<float>" );
}

template<>
void cublasWrapperSwap( int n, double* x_d, int incX, double* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasDswap( CUDAContext_cublasHandle, n, x_d, incX, y_d, incY ), "cublasWrapperSwap<double>" );
}

template<>
void cublasWrapperSwap( int n, ComplexFloat* x_d, int incX, ComplexFloat* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasCswap( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, cublasCast( y_d ), incY ),
                      "cublasWrapperSwap<ComplexFloat>" );
}

template<>
void cublasWrapperSwap( int n, ComplexDouble* x_d, int incX, ComplexDouble* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasZswap( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, cublasCast( y_d ), incY ),
                      "cublasWrapperSwap<ComplexDouble>" );
}

template<typename ValueType>
void CUDABLAS1::swap(
    const IndexType n,
    ValueType* x_d,
    const IndexType incX,
    ValueType* y_d,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.swap" )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "swap<" << getScalarType<ValueType>() << "> of x, y with size " << n )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS::swap set stream" );

    cublasWrapperSwap( static_cast<int>( n ), x_d, static_cast<int>( incX ), y_d, static_cast<int>( incY ) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUADABLAS1::swap set stream NULL" );
}

/* ---------------------------------------------------------------------------------------*/
/*    copy                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline void cublasWrapperCopy( int n, const ValueType* x_d, int incX, ValueType* y_d, int incY );

template<>
void cublasWrapperCopy( int n, const float* x_d, int incX, float* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasScopy( CUDAContext_cublasHandle, n, x_d, incX, y_d, incY ), "cublasWrapperCopy<float>" );
}

template<>
void cublasWrapperCopy( int n, const double* x_d, int incX, double* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasDcopy( CUDAContext_cublasHandle, n, x_d, incX, y_d, incY ), "cublasWrapperCopy<double>" );
}

template<>
void cublasWrapperCopy( int n, const ComplexFloat* x_d, int incX, ComplexFloat* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasCcopy( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, cublasCast( y_d ), incY ),
                      "cublasWrapperCopy<ComplexFloat>" );
}

template<>
void cublasWrapperCopy( int n, const ComplexDouble* x_d, int incX, ComplexDouble* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasZcopy( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, cublasCast( y_d ), incY ),
                      "cublasWrapperCopy<ComplexDouble>" );
}

template<typename ValueType>
void CUDABLAS1::copy(
    IndexType n,
    const ValueType* x_d,
    IndexType incX,
    ValueType* y_d,
    IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.copy" )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "copy<" << getScalarType<ValueType>() << "> of x, y, n = " << n )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::copy set stream" );

    cublasWrapperCopy( static_cast<int>( n ), x_d, static_cast<int>( incX ), y_d, static_cast<int>( incY ) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::copy set stream NULL" );
}

/* ---------------------------------------------------------------------------------------*/
/*    axpy                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline void cublasWrapperAxpy(
    int n,
    ValueType alpha,
    const ValueType* x_d,
    int incX,
    ValueType* y_d,
    int incY );

template<>
void cublasWrapperAxpy( int n, float alpha, const float* x_d, int incX, float* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasSaxpy( CUDAContext_cublasHandle, n, &alpha, x_d, incX, y_d, incY ),
                      "cublasWrapperAxpy<float>" );
}

template<>
void cublasWrapperAxpy( int n, double alpha, const double* x_d, int incX, double* y_d, int incY )
{
    LAMA_CUBLAS_CALL( cublasDaxpy( CUDAContext_cublasHandle, n, &alpha, x_d, incX, y_d, incY ),
                      "cublasWrapperAxpy<double>" );
}

template<>
void cublasWrapperAxpy( int n, ComplexFloat alpha, const ComplexFloat* x_d, int incX, ComplexFloat* y_d, int incY )
{
    LAMA_CUBLAS_CALL(
        cublasCaxpy( CUDAContext_cublasHandle, n, cublasCast( &alpha ), cublasCast( x_d ), incX,
                     cublasCast( y_d ), incY ),
        "cublasWrapperAxpy<ComplexFloat>" );
}

template<>
void cublasWrapperAxpy( int n, ComplexDouble alpha, const ComplexDouble* x_d, int incX, ComplexDouble* y_d, int incY )
{
    LAMA_CUBLAS_CALL(
        cublasZaxpy( CUDAContext_cublasHandle, n, cublasCast( &alpha ), cublasCast( x_d ), incX,
                     cublasCast( y_d ), incY ),
        "cublasWrapperAxpy<ComplexDouble>" );
}

template<typename ValueType>
void CUDABLAS1::axpy(
    int n,
    ValueType alpha,
    const ValueType* x_d,
    int incX,
    ValueType* y_d,
    const int incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.axpy" )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "axpy<" << getScalarType<ValueType>() << "> of x, y, n = " << n << ", alpha = " << alpha )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::axpy set stream" );

    cublasWrapperAxpy( n, alpha, x_d, incX, y_d, incY );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::axpy set stream NULL" );
}

/* ---------------------------------------------------------------------------------------*/
/*    dot                                                                                 */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
static inline ValueType cublasWrapperDot( int n, const ValueType* x_d, int incX, const ValueType* y_d, int incY );

template<>
float cublasWrapperDot( int n, const float* x_d, int incX, const float* y_d, int incY )
{
    float dot;
    LAMA_CUBLAS_CALL( cublasSdot( CUDAContext_cublasHandle, n, x_d, incX, y_d, incY, &dot ), "cublasWrapperDot<float>" );
    return dot;
}

template<>
double cublasWrapperDot( int n, const double* x_d, int incX, const double* y_d, int incY )
{
    double dot;
    LAMA_CUBLAS_CALL( cublasDdot( CUDAContext_cublasHandle, n, x_d, incX, y_d, incY, &dot ),
                      "cublasWrapperDot<double>" );
    return dot;
}

template<>
ComplexFloat cublasWrapperDot( int n, const ComplexFloat* x_d, int incX, const ComplexFloat* y_d, int incY )
{
    ComplexFloat dot;
    LAMA_CUBLAS_CALL(
        cublasCdotu( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, cublasCast( y_d ), incY,
                     cublasCast( &dot ) ),
        "cublasWrapperDot<ComplexFloat>" );
    return dot;
}

template<>
ComplexDouble cublasWrapperDot( int n, const ComplexDouble* x_d, int incX, const ComplexDouble* y_d, int incY )
{
    ComplexDouble dot;
    LAMA_CUBLAS_CALL(
        cublasZdotu( CUDAContext_cublasHandle, n, cublasCast( x_d ), incX, cublasCast( y_d ), incY,
                     cublasCast( &dot ) ),
        "cublasWrapperDot<ComplexDouble>" );
    return dot;
}

template<typename ValueType>
ValueType CUDABLAS1::dot(
    IndexType n,
    const ValueType* x_d,
    IndexType incX,
    const ValueType* y_d,
    IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.dot" )

    LAMA_LOG_DEBUG( logger,
                    "dot<" << getScalarType<ValueType>() << ">, n = " << n << ", incX = " << incX << ", incY = " << incY << ", x_d = " << x_d << ", y_d = " << y_d )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return 0.0;
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "CUDABLAS1::dot set stream" );

    ValueType res = cublasWrapperDot( static_cast<int>( n ), x_d, static_cast<int>( incX ), y_d,
                                      static_cast<int>( incY ) );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, NULL ), "CUDABLAS1::dot set stream NULL" );
    return res;
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
    ValueType* z,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.sum" )

    if( n <= 0 )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger,
                    "sum<" << getScalarType<ValueType>() << ">, n = " << n << ", " << alpha << " * x + " << beta << " * y " )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    sum_launcher( n, alpha, x, beta, y, z, stream );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        cudaStreamSynchronize( stream );
        LAMA_CHECK_CUDA_ERROR
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDABLAS1::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions
    //       ( e.g. routine CUDABLAS1::sum<ValueType> is set for BLAS::BLAS1::sum variable

#define LAMA_BLAS1_REGISTER(z, I, _)                                            \
    LAMA_INTERFACE_REGISTER_T( BLAS, scal, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, nrm2, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, asum, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, iamax, ARITHMETIC_TYPE##I )                \
    LAMA_INTERFACE_REGISTER_T( BLAS, swap, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, copy, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, axpy, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, dot, ARITHMETIC_TYPE##I )                  \
    LAMA_INTERFACE_REGISTER_T( BLAS, sum, ARITHMETIC_TYPE##I )                  \

    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_BLAS1_REGISTER, _ )

#undef LAMA_BLAS1_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDABLAS1::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( context::CUDA );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDABLAS1::initialized = registerInterface();

} /* namespace lama */
