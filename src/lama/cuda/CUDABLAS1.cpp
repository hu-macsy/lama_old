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
 * @brief Wrapper implementations for BLAS1 routines in CUDA using cuBLAS
 * @author Lauretta Schubert, Thomas Brandes, Eric Stricker
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/cuda/CUDABLAS1.hpp>

// others
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>
#include <lama/cuda/lama_cublas.hpp>
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

/* ---------------------------------------------------------------------------------------*/
/*    scale                                                                               */
/* ---------------------------------------------------------------------------------------*/

// Note: the cublasWrapper routines could be static routines on its own. But using
//       a common template routine is helpful to guarantee correct syntax

template<typename T>
static inline void cublasWrapperScale( int n, T alpha, T* x_d, int incX );

template<>
void cublasWrapperScale( int n, float alpha, float* x_d, int incX )
{
    cublasSscal( n, alpha, x_d, incX );
}

template<>
void cublasWrapperScale( int n, double alpha, double* x_d, int incX )
{
    cublasDscal( n, alpha, x_d, incX );
}

template<>
void cublasWrapperScale( int n, ComplexFloat alpha, ComplexFloat* x_d, int incX )
{
    // use of cublasCast to convert ComplexFloat to cuFloatComplex via reinterpret_cast
    cublasCscal( n, cublasCast( alpha ), cublasCast( x_d ), incX );
}

template<>
void cublasWrapperScale( int n, ComplexDouble alpha, ComplexDouble* x_d, int incX )
{
    // use of cublasCast to convert ComplexDouble to cuDoubleComplex via reinterpret_cast
    cublasZscal( n, cublasCast( alpha ), cublasCast( x_d ), incX );
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

    cublasWrapperScale( static_cast<int>( n ), alpha, x_d, static_cast<int>( incX ) );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/* ---------------------------------------------------------------------------------------*/
/*    nrm2                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename T>
static inline T cublasWrapperNrm2( int n, const T* x_d, int incX );

template<>
float cublasWrapperNrm2( int n, const float* x_d, int incX )
{
    return cublasSnrm2( n, x_d, incX );
}

template<>
double cublasWrapperNrm2( int n, const double* x_d, int incX )
{
    return cublasDnrm2( n, x_d, incX );
}

template<>
ComplexFloat cublasWrapperNrm2( int n, const ComplexFloat* x_d, int incX )
{
    // CUBLAS returns only float result so we convert it back to Complex

    float res = ComplexFloat( cublasScnrm2( n, cublasCast( x_d ), incX ) );

    return ComplexFloat( res, 0.0f );
}

template<>
ComplexDouble cublasWrapperNrm2( int n, const ComplexDouble* x_d, int incX )
{
    // CUBLAS returns only double result so we convert it back to Complex

    double res = ComplexFloat( cublasDznrm2( n, cublasCast( x_d ), incX ) );

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

    T res = cublasWrapperNrm2( static_cast<int>( n ), x_d, static_cast<int>( incX ) );

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

/* ---------------------------------------------------------------------------------------*/
/*    asum                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename T>
static inline T cublasWrapperAsum( int n, const T* x_d, int incX );

template<>
float cublasWrapperAsum( int n, const float* x_d, int incX )
{
    return cublasSasum( n, x_d, incX );
}

template<>
double cublasWrapperAsum( int n, const double* x_d, int incX )
{
    return cublasDasum( n, x_d, incX );
}

template<>
ComplexFloat cublasWrapperAsum( int n, const ComplexFloat* x_d, int incX )
{
    // CUBLAS returns only float result so we convert it back to Complex

    float res = ComplexFloat( cublasScasum( n, cublasCast( x_d ), incX ) );

    return ComplexFloat( res );
}

template<>
ComplexDouble cublasWrapperAsum( int n, const ComplexDouble* x_d, int incX )
{
    // CUBLAS returns only double result so we convert it back to Complex

    double res = ComplexFloat( cublasDzasum( n, cublasCast( x_d ), incX ) );

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

    T res = cublasWrapperAsum( static_cast<int>( n ), x_d, static_cast<int>( incX ) );

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

/* ---------------------------------------------------------------------------------------*/
/*    iamax                                                                               */
/* ---------------------------------------------------------------------------------------*/

template<typename T>
static int cublasWrapperIamax( int n, const T* x_d, int incX );

template<>
int cublasWrapperIamax( int n, const float* x_d, int incX )
{
    int iamax = cublasIsamax( n, x_d, incX );
    return iamax;
}

template<>
int cublasWrapperIamax( int n, const double* x_d, int incX )
{
    int iamax = cublasIdamax( n, x_d, incX );
    return iamax;
}

template<>
int cublasWrapperIamax( int n, const ComplexFloat* x_d, int incX )
{
    int iamax = cublasIcamax( n, cublasCast( x_d ), incX );
    return iamax;
}

template<>
int cublasWrapperIamax( int n, const ComplexDouble* x_d, int incX )
{
    int iamax = cublasIzamax( n, cublasCast( x_d ), incX );
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

    IndexType iamax = cublasWrapperIamax( n, x_d, incX );

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

/* ---------------------------------------------------------------------------------------*/
/*    swap                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename T>
static inline void cublasWrapperSwap( IndexType n, T* x_d, IndexType incX, T* y_d, IndexType incY  );

template<>
void cublasWrapperSwap( int n, float* x_d, int incX, float* y_d, int incY )
{
    cublasSswap( n, x_d, incX, y_d, incY );
}

template<>
void cublasWrapperSwap( int n, double* x_d, int incX, double* y_d, int incY )
{
    cublasDswap( n, x_d, incX, y_d, incY );
}

template<>
void cublasWrapperSwap( int n, ComplexFloat* x_d, int incX, ComplexFloat* y_d, int incY )
{
    cublasCswap( n, cublasCast( x_d ), incX, cublasCast( y_d ), incY );
}

template<>
void cublasWrapperSwap( int n, ComplexDouble* x_d, int incX, ComplexDouble* y_d, int incY )
{
    cublasZswap( n, cublasCast( x_d ), incX, cublasCast( y_d ), incY );
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

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
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

    cublasWrapperSwap( static_cast<int>( n ), x_d, static_cast<int>( incX ), y_d, static_cast<int>( incY ) );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/* ---------------------------------------------------------------------------------------*/
/*    copy                                                                                */
/* ---------------------------------------------------------------------------------------*/


template<typename T>
static inline void cublasWrapperCopy( int n, const T* x_d, int incX, T* y_d, int incY  );

template<>
void cublasWrapperCopy( int n, const float* x_d, int incX, float* y_d, int incY )
{
    cublasScopy( n, x_d, incX, y_d, incY );
}

template<>
void cublasWrapperCopy( int n, const double* x_d, int incX, double* y_d, int incY )
{
    cublasDcopy( n, x_d, incX, y_d, incY );
}

template<>
void cublasWrapperCopy( int n, const ComplexFloat* x_d, int incX, ComplexFloat* y_d, int incY )
{
    cublasCcopy( n, cublasCast( x_d ), incX, cublasCast( y_d ), incY );
}

template<>
void cublasWrapperCopy( int n, const ComplexDouble* x_d, int incX, ComplexDouble* y_d, int incY )
{
    cublasZcopy( n, cublasCast( x_d ), incX, cublasCast( y_d ), incY );
}

template<typename T>
void CUDABLAS1::copy( IndexType n, const T* x_d, IndexType incX, T* y_d, IndexType incY, SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.copy" )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
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

    cublasWrapperCopy( static_cast<int>( n ), x_d, static_cast<int>( incX ), y_d, static_cast<int>( incY ) );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/* ---------------------------------------------------------------------------------------*/
/*    axpy                                                                                */
/* ---------------------------------------------------------------------------------------*/


template<typename T>
static inline void cublasWrapperAxpy( int n, T alpha, const T* x_d, int incX, T* y_d, int incY  );

template<>
void cublasWrapperAxpy( int n, float alpha, const float* x_d, int incX, float* y_d, int incY )
{
    cublasSaxpy( n, alpha, x_d, incX, y_d, incY );
}

template<>
void cublasWrapperAxpy( int n, double alpha, const double* x_d, int incX, double* y_d, int incY )
{
    cublasDaxpy( n, alpha, x_d, incX, y_d, incY );
}

template<>
void cublasWrapperAxpy( int n, ComplexFloat alpha, 
                  const ComplexFloat* x_d, int incX,
                  ComplexFloat* y_d, int incY )
{
    cublasCaxpy( n, cublasCast( alpha ), cublasCast( x_d ), incX, cublasCast( y_d ), incY );
}

template<>
void cublasWrapperAxpy( int n, ComplexDouble alpha, 
                  const ComplexDouble* x_d, int incX,
                  ComplexDouble* y_d, int incY )
{
    cublasZaxpy( n, cublasCast( alpha ), cublasCast( x_d ), incX, cublasCast( y_d ), incY );
}

template<typename T>
void CUDABLAS1::axpy( int n, T alpha,
                      const T* x_d, int incX,
                      T* y_d, const int incY,
                      SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.axpy" )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
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

    cublasWrapperAxpy( n, alpha, x_d, incX, y_d, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/* ---------------------------------------------------------------------------------------*/
/*    dot                                                                                 */
/* ---------------------------------------------------------------------------------------*/

template<typename T>
static inline T cublasWrapperDot( int n, const T* x_d, int incX, const T* y_d, int incY  );

template<>
float cublasWrapperDot( int n, const float* x_d, int incX, const float* y_d, int incY )
{
    return cublasSdot( n, x_d, incX, y_d, incY );
}

template<>
double cublasWrapperDot( int n, const double* x_d, int incX, const double* y_d, int incY )
{
    return cublasDdot( n, x_d, incX, y_d, incY );
}

template<>
ComplexFloat cublasWrapperDot( int n, 
                         const ComplexFloat* x_d, int incX,
                         const ComplexFloat* y_d, int incY )
{
    cuFloatComplex dotResult = cublasCdotu ( n, cublasCast( x_d ), incX, cublasCast( y_d ), incY );

    return ComplexFloat( dotResult.x, dotResult.y );
}

template<>
ComplexDouble cublasWrapperDot( int n, const ComplexDouble* x_d, int incX,
                          const ComplexDouble* y_d, int incY )
{
    cuDoubleComplex dotResult = cublasZdotu( n, cublasCast( x_d ), incX, cublasCast( y_d ), incY );

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

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
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

    T res = cublasWrapperDot( static_cast<int>( n ), x_d, static_cast<int>( incX ), y_d, static_cast<int>( incY ) );

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

/* ---------------------------------------------------------------------------------------*/
/*    sum                                                                                 */
/* ---------------------------------------------------------------------------------------*/

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
