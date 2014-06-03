/**
 * @file CUDABLAS1.cpp
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
 * @brief CUDABLAS1.cpp
 * @author lschubert
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/cuda/CUDABLAS1.hpp>

// others
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// macros
#include <lama/macros/unused.hpp>

// tracing with LAMA_REGION
#include <lama/tracing.hpp>

// blas
#include <boost/preprocessor.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDABLAS1::logger, "CUDA.BLAS1" )

/** scale */

template<typename T>
static inline void wrapperScale( IndexType n, T alpha, T* x_d, IndexType incX );

template<>
void wrapperScale( IndexType n, float alpha, float* x_d, IndexType incX )
{
    cublasSscal( n, alpha, x_d, incX );
}

template<>
void wrapperScale( IndexType n, double alpha, double* x_d, IndexType incX )
{
    cublasDscal( n, alpha, x_d, incX );
}

template<>
void wrapperScale( IndexType n, ComplexFloat alpha, ComplexFloat* x_d, IndexType incX )
{
    cuFloatComplex* xc_d = reinterpret_cast<cuFloatComplex*>( x_d );

    const cuFloatComplex* alphac = reinterpret_cast<const cuFloatComplex*>( &alpha );

    cublasCscal( n, *alphac, xc_d, incX );
}

template<>
void wrapperScale( IndexType n, ComplexDouble alpha, ComplexDouble* x_d, IndexType incX )
{
    cuDoubleComplex* xc_d = reinterpret_cast<cuDoubleComplex*>( x_d );

    const cuDoubleComplex* alphac = reinterpret_cast<const cuDoubleComplex*>( &alpha );

    cublasZscal( n, *alphac, xc_d, incX );
}

template<typename T>
void CUDABLAS1::scal( IndexType n, const T alpha, T* x_d, const IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.scal" )

    if ( incX == 0 )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "scal<" << Scalar::getType<T>() << "> of x[" << n << "], alpha = " << alpha )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    cublasSetKernelStream( stream );
    LAMA_CHECK_CUBLAS_ERROR

    wrapperScale( n, alpha, x_d, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/** nrm2 */

template<typename T>
static inline T wrapperNrm2( IndexType n, const T* x_d, IndexType incX );

template<>
float wrapperNrm2( IndexType n, const float* x_d, IndexType incX )
{
    return cublasSnrm2( n, x_d, incX );
}

template<>
double wrapperNrm2( IndexType n, const double* x_d, IndexType incX )
{
    return cublasDnrm2( n, x_d, incX );
}

template<>
ComplexFloat wrapperNrm2( IndexType n, const ComplexFloat* x_d, IndexType incX )
{
    const cuFloatComplex* xc_d = reinterpret_cast<const cuFloatComplex*>( x_d );

    // CUBLAS returns only float result so we convert it back to Complex

    float res = ComplexFloat( cublasScnrm2( n, xc_d, incX ) );

    return ComplexFloat( res, 0.0f );
}

template<>
ComplexDouble wrapperNrm2( IndexType n, const ComplexDouble* x_d, IndexType incX )
{
    const cuDoubleComplex* xc_d = reinterpret_cast<const cuDoubleComplex*>( x_d );

    // CUBLAS returns only double result so we convert it back to Complex

    double res = ComplexFloat( cublasDznrm2( n, xc_d, incX ) );

    return ComplexDouble( res, 0.0 );
}

template<typename T>
T CUDABLAS1::nrm2( IndexType n, const T* x_d, IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.nrm2" )

    if ( incX <= 0 )
    {
        return 0.0;
    }

    LAMA_LOG_DEBUG( logger, "nrm2<" << Scalar::getType<T>() << "> of x[" << n << "]" )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    cublasSetKernelStream( stream );
    LAMA_CHECK_CUBLAS_ERROR

    T res = wrapperNrm2( n, x_d, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR

    return res;
}

/** asum */

template<typename T>
static inline T wrapperAsum( IndexType n, const T* x_d, IndexType incX );

template<>
float wrapperAsum( IndexType n, const float* x_d, IndexType incX )
{
    return cublasSasum( n, x_d, incX );
}

template<>
double wrapperAsum( IndexType n, const double* x_d, IndexType incX )
{
    return cublasDasum( n, x_d, incX );
}

template<>
ComplexFloat wrapperAsum( IndexType n, const ComplexFloat* x_d, IndexType incX )
{
    const cuFloatComplex* xc_d = reinterpret_cast<const cuFloatComplex*>( x_d );

    // CUBLAS returns only float result so we convert it back to Complex
    float res = ComplexFloat( cublasScasum( n, xc_d, incX ) );

    return ComplexFloat( res );
}

template<>
ComplexDouble wrapperAsum( IndexType n, const ComplexDouble* x_d, IndexType incX )
{
    const cuDoubleComplex* xc_d = reinterpret_cast<const cuDoubleComplex*>( x_d );

    // CUBLAS returns only double result so we convert it back to Complex

    double res = ComplexFloat( cublasDzasum( n, xc_d, incX ) );
    return ComplexDouble( res );
}

template<typename T>
T CUDABLAS1::asum( const IndexType n, const T* x_d, const IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.asum" )

    if ( incX <= 0 )
    {
        return 0.0;
    }

    LAMA_LOG_DEBUG( logger, "asum<" << Scalar::getType<T>() << "> of x[" << n << "]" )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    cublasSetKernelStream( stream );
    LAMA_CHECK_CUBLAS_ERROR

    T res = wrapperAsum( n, x_d, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR

    return res;
}

/** iamax */

template<typename T>
static IndexType wrapperIamax( IndexType n, const T* x_d, IndexType incX );

template<>
IndexType wrapperIamax( IndexType n, const float* x_d, IndexType incX )
{
    IndexType iamax = cublasIsamax( n, x_d, incX );
    return iamax;
}

template<>
IndexType wrapperIamax( IndexType n, const double* x_d, IndexType incX )
{
    IndexType iamax = cublasIdamax( n, x_d, incX );
    return iamax;
}

template<>
IndexType wrapperIamax( IndexType n, const ComplexFloat* x_d, IndexType incX )
{
    const cuFloatComplex* xc_d = reinterpret_cast<const cuFloatComplex*>( x_d );
    IndexType iamax = cublasIcamax( n, xc_d, incX );
    return iamax;
}

template<>
IndexType wrapperIamax( IndexType n, const ComplexDouble* x_d, IndexType incX )
{
    const cuDoubleComplex* xc_d = reinterpret_cast<const cuDoubleComplex*>( x_d );
    IndexType iamax = cublasIzamax( n, xc_d, incX );
    return iamax;
}

template<typename T>
IndexType CUDABLAS1::iamax( const IndexType n, const T* x_d, const IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.iamax" )

    LAMA_LOG_DEBUG( logger, "iamax<" << Scalar::getType<T>() << "> of x[" << n << "]" )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    cublasSetKernelStream( stream );
    LAMA_CHECK_CUBLAS_ERROR

    IndexType iamax = wrapperIamax( n, x_d, incX );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR

    return iamax ? iamax - 1 : 0;
}

/** swap */

template<typename T>
static inline void wrapperSwap( IndexType n, T* x_d, IndexType incX, T* y_d, IndexType incY  );

template<>
void wrapperSwap( IndexType n, float* x_d, IndexType incX, float* y_d, IndexType incY )
{
    cublasSswap( n, x_d, incX, y_d, incY );
}

template<>
void wrapperSwap( IndexType n, double* x_d, IndexType incX, double* y_d, IndexType incY )
{
    cublasDswap( n, x_d, incX, y_d, incY );
}

template<>
void wrapperSwap( IndexType n, ComplexFloat* x_d, IndexType incX, ComplexFloat* y_d, IndexType incY )
{
    cuFloatComplex* xc_d = reinterpret_cast<cuFloatComplex*>( x_d );
    cuFloatComplex* yc_d = reinterpret_cast<cuFloatComplex*>( y_d );

    cublasCswap( n, xc_d, incX, yc_d, incY );
}

template<>
void wrapperSwap( IndexType n, ComplexDouble* x_d, IndexType incX, ComplexDouble* y_d, IndexType incY )
{
    cuDoubleComplex* xc_d = reinterpret_cast<cuDoubleComplex*>( x_d );
    cuDoubleComplex* yc_d = reinterpret_cast<cuDoubleComplex*>( y_d );

    cublasZswap( n, xc_d, incX, yc_d, incY );
}

template<typename T>
void CUDABLAS1::swap(
    const IndexType n,
    T* x_d,
    const IndexType incX,
    T* y_d,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.swap" )

    if ( (incX <= 0) || (incY <= 0) )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "swap<" << Scalar::getType<T>() << "> of x, y with size " << n )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    cublasSetKernelStream( stream );
    LAMA_CHECK_CUBLAS_ERROR

    wrapperSwap( n, x_d, incX, y_d, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/** copy */

template<typename T>
static inline void wrapperCopy( IndexType n, const T* x_d, IndexType incX, T* y_d, IndexType incY  );

template<>
void wrapperCopy( IndexType n, const float* x_d, IndexType incX, float* y_d, IndexType incY )
{
    cublasScopy( n, x_d, incX, y_d, incY );
}

template<>
void wrapperCopy( IndexType n, const double* x_d, IndexType incX, double* y_d, IndexType incY )
{
    cublasDcopy( n, x_d, incX, y_d, incY );
}

template<>
void wrapperCopy( IndexType n, const ComplexFloat* x_d, IndexType incX, ComplexFloat* y_d, IndexType incY )
{
    const cuFloatComplex* xc_d = reinterpret_cast<const cuFloatComplex*>( x_d );
    cuFloatComplex* yc_d = reinterpret_cast<cuFloatComplex*>( y_d );

    cublasCcopy( n, xc_d, incX, yc_d, incY );
}

template<>
void wrapperCopy( IndexType n, const ComplexDouble* x_d, IndexType incX, ComplexDouble* y_d, IndexType incY )
{
    const cuDoubleComplex* xc_d = reinterpret_cast<const cuDoubleComplex*>( x_d );
    cuDoubleComplex* yc_d = reinterpret_cast<cuDoubleComplex*>( y_d );

    cublasZcopy( n, xc_d, incX, yc_d, incY );
}

template<typename T>
void CUDABLAS1::copy( IndexType n, const T* x_d, IndexType incX, T* y_d, IndexType incY, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.copy" )

    if ( (incX <= 0) || (incY <= 0) )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "copy<" << Scalar::getType<T>() << "> of x, y, n = " << n )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    cublasSetKernelStream( stream );
    LAMA_CHECK_CUBLAS_ERROR

    wrapperCopy( n, x_d, incX, y_d, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/** axpy */

template<typename T>
static inline void wrapperAxpy( IndexType n, T alpha, const T* x_d, IndexType incX, T* y_d, IndexType incY  );

template<>
void wrapperAxpy( IndexType n, float alpha, const float* x_d, IndexType incX, float* y_d, IndexType incY )
{
    cublasSaxpy( n, alpha, x_d, incX, y_d, incY );
}

template<>
void wrapperAxpy( IndexType n, double alpha, const double* x_d, IndexType incX, double* y_d, IndexType incY )
{
    cublasDaxpy( n, alpha, x_d, incX, y_d, incY );
}

template<>
void wrapperAxpy( IndexType n, ComplexFloat alpha, const ComplexFloat* x_d, IndexType incX, 
                                                   ComplexFloat* y_d, IndexType incY )
{
    const cuFloatComplex* alphac = reinterpret_cast<const cuFloatComplex*>( &alpha );

    const cuFloatComplex* xc_d = reinterpret_cast<const cuFloatComplex*>( x_d );
    cuFloatComplex* yc_d = reinterpret_cast<cuFloatComplex*>( y_d );

    cublasCaxpy( n, *alphac, xc_d, incX, yc_d, incY );
}

template<>
void wrapperAxpy( IndexType n, ComplexDouble alpha, const ComplexDouble* x_d, IndexType incX, 
                                                   ComplexDouble* y_d, IndexType incY )
{
    const cuDoubleComplex* alphac = reinterpret_cast<const cuDoubleComplex*>( &alpha );

    const cuDoubleComplex* xc_d = reinterpret_cast<const cuDoubleComplex*>( x_d );
    cuDoubleComplex* yc_d = reinterpret_cast<cuDoubleComplex*>( y_d );

    cublasZaxpy( n, *alphac, xc_d, incX, yc_d, incY );
}

template<typename T>
void CUDABLAS1::axpy( IndexType n, T alpha,
                      const T* x_d, IndexType incX,
                      T* y_d, const IndexType incY,
                      SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.axpy" )

    if ( (incX <= 0) || (incY <= 0) )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "axpy<" << Scalar::getType<T>() << "> of x, y, n = " << n
                            << ", alpha = " << alpha )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    cublasSetKernelStream( stream );
    LAMA_CHECK_CUBLAS_ERROR

    wrapperAxpy( n, alpha, x_d, incX, y_d, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/** dot */

template<typename T>
static inline T wrapperDot( IndexType n, const T* x_d, IndexType incX, const T* y_d, IndexType incY  );

template<>
float wrapperDot( IndexType n, const float* x_d, IndexType incX, const float* y_d, IndexType incY )
{
    return cublasSdot( n, x_d, incX, y_d, incY );
}

template<>
double wrapperDot( IndexType n, const double* x_d, IndexType incX, const double* y_d, IndexType incY )
{
    return cublasDdot( n, x_d, incX, y_d, incY );
}

template<>
ComplexFloat wrapperDot( IndexType n, const ComplexFloat* x_d, IndexType incX, 
                                      const ComplexFloat* y_d, IndexType incY )
{
    const cuFloatComplex* xc_d = reinterpret_cast<const cuFloatComplex*>( x_d );
    const cuFloatComplex* yc_d = reinterpret_cast<const cuFloatComplex*>( y_d );

    cuFloatComplex dotResult = cublasCdotu ( n, xc_d, incX, yc_d, incY );

    return ComplexFloat( dotResult.x, dotResult.y );
}

template<>
ComplexDouble wrapperDot( IndexType n, const ComplexDouble* x_d, IndexType incX, 
                                       const ComplexDouble* y_d, IndexType incY )
{
    const cuDoubleComplex* xc_d = reinterpret_cast<const cuDoubleComplex*>( x_d );
    const cuDoubleComplex* yc_d = reinterpret_cast<const cuDoubleComplex*>( y_d );

    cuDoubleComplex dotResult = cublasZdotu( n, xc_d, incX, yc_d, incY );

    return ComplexDouble( dotResult.x, dotResult.y );
}

template<typename T>
T CUDABLAS1::dot(
    IndexType n,
    const T* x_d,
    IndexType incX,
    const T* y_d,
    IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.dot" )

    LAMA_LOG_DEBUG( logger, "dot<" << Scalar::getType<T>() << ">, n = " << n 
                             << ", incX = " << incX << ", incY = " << incY 
                             << ", x_d = " << x_d << ", y_d = " << y_d )

    if ( (incX <= 0) || (incY <= 0) )
    {
        return 0.0;
    }

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    cublasSetKernelStream( stream );
    LAMA_CHECK_CUBLAS_ERROR

    T res = wrapperDot( n, x_d, incX, y_d, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR

    return res;
}

/** sum */
template<typename T>
void CUDABLAS1::sum( const IndexType n, T alpha, const T* x, T beta, const T* y, T* z, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.sum" )

    if ( n <= 0 )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "sum<" << Scalar::getType<T>() << ">, n = " << n 
                            << ", " << alpha << " * x + " << beta << " * y " )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    sum_launcher( n, alpha, x, beta, y, z, stream );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
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
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

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
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDABLAS1::initialized = registerInterface();

} /* namespace lama */
