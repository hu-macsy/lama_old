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

template<>
void CUDABLAS1::scal( IndexType n, const float alpha, float* x_d, const IndexType incx, SyncToken* syncToken )
{
    if ( incx == 0 )
    {
        return;
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

    cublasSscal( n, alpha, x_d, incx );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

template<>
void CUDABLAS1::scal( IndexType n, const double alpha, double* x_d, const IndexType incx, SyncToken* syncToken )
{
    if ( incx <= 0 )
    {
        return;
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

    cublasDscal( n, alpha, x_d, incx );

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

template<>
float CUDABLAS1::nrm2( IndexType n, const float* x_d, IndexType incx, SyncToken* syncToken )
{
    if ( incx <= 0 )
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

    float res = cublasSnrm2( n, x_d, incx );

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

template<>
double CUDABLAS1::nrm2( IndexType n, const double* x_d, IndexType incx, SyncToken* syncToken )
{
    if ( incx <= 0 )
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

    double res = cublasDnrm2( n, x_d, incx );

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
template<>
float CUDABLAS1::asum( const IndexType n, const float* x_d, const IndexType incX, SyncToken* syncToken )
{
    if ( incX <= 0 )
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

    float res = cublasSasum( n, x_d, incX );

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

template<>
double CUDABLAS1::asum( const IndexType n, const double* x_d, const IndexType incX, SyncToken* syncToken )
{
    if ( incX <= 0 )
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

    double res = cublasDasum( n, x_d, incX );

    // No error check here possible as kernel is started asynchronously

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR

    return res;
}

/** iamax */

template<>
IndexType CUDABLAS1::iamax( const IndexType n, const float* x_d, const IndexType incX, SyncToken* syncToken )
{
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

    IndexType iamax = cublasIsamax( n, x_d, incX );

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

template<>
IndexType CUDABLAS1::iamax( const IndexType n, const double* x_d, const IndexType incX, SyncToken* syncToken )
{
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

    IndexType iamax = cublasIdamax( n, x_d, incX );

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

template<>
void CUDABLAS1::swap(
    const IndexType n,
    float* x_d,
    const IndexType incX,
    float* y_d,
    const IndexType incY,
    SyncToken* syncToken )
{
    if ( (incX <= 0) || (incY <= 0) )
    {
        return;
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

    cublasSswap( n, x_d, incX, y_d, incY );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

template<>
void CUDABLAS1::swap(
    const IndexType n,
    double* x_d,
    const IndexType incX,
    double* y_d,
    const IndexType incY,
    SyncToken* syncToken )
{
    if ( (incX <= 0) || (incY <= 0) )
    {
        return;
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

    cublasDswap( n, x_d, incX, y_d, incY );

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

template<>
void CUDABLAS1::copy( IndexType n, const float* x_d, IndexType incx, float* y_d, IndexType incy, SyncToken* syncToken )
{
    if ( (incx <= 0) || (incy <= 0) )
    {
        return;
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

    cublasScopy( n, x_d, incx, y_d, incy );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

template<>
void CUDABLAS1::copy(
    IndexType n,
    const double* x_d,
    IndexType incx,
    double* y_d,
    IndexType incy,
    SyncToken* syncToken )
{
    if ( (incx <= 0) || (incy <= 0) )
    {
        return;
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

    cublasDcopy( n, x_d, incx, y_d, incy );

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

template<>
void CUDABLAS1::axpy(
    IndexType n,
    float alpha,
    const float* x_d,
    IndexType incx,
    float* y_d,
    const IndexType incy,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.saxpy" )

    if ( (incx <= 0) || (incy <= 0) )
    {
        return;
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

    cublasSaxpy( n, alpha, x_d, incx, y_d, incy );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

template<>
void CUDABLAS1::axpy(
    IndexType n,
    double alpha,
    const double* x_d,
    IndexType incx,
    double* y_d,
    const IndexType incy,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.daxpy" )

    if ( (incx <= 0) || (incy <= 0) )
    {
        return;
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

    cublasDaxpy( n, alpha, x_d, incx, y_d, incy );

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

template<>
float CUDABLAS1::dot(
    IndexType n,
    const float* x_d,
    IndexType incx,
    const float* y_d,
    IndexType incy,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.sdot" )

    LAMA_LOG_DEBUG( logger, "CUDABLAS1:sdot, n = " << n << ", incx = " << incx << ", incy = " << incy 
                             << ", x_d = " << x_d << ", y_d = " << y_d )

    if ( (incx <= 0) || (incy <= 0) )
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

    float res = cublasSdot( n, x_d, incx, y_d, incy );

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

template<>
double CUDABLAS1::dot(
    IndexType n,
    const double* x_d,
    IndexType incx,
    const double* y_d,
    IndexType incy,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.BLAS1.ddot" )

    LAMA_LOG_DEBUG( logger, "CUDABLAS1:ddot, n = " << n << ", incx = " << incx << ", incy = " << incy 
                            << ", x_d = " << x_d << ", y_d = " << y_d )

    if ( (incx <= 0) || (incy <= 0) )
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

    double res = cublasDdot( n, x_d, incx, y_d, incy );

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
    if ( n <= 0 )
    {
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

    sum_launcher( n, alpha, x, beta, y, z, stream );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( stream );
        LAMA_CHECK_CUDA_ERROR
    }
}

/** rot */

template<>
void CUDABLAS1::rot(
    const IndexType n,
    float* x_d,
    const IndexType incX,
    float* y_d,
    const IndexType incY,
    const float c,
    const float s,
    SyncToken* syncToken )
{
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

    cublasSrot( n, x_d, incX, y_d, incY, c, s );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

template<>
void CUDABLAS1::rot(
    const IndexType n,
    double* x_d,
    const IndexType incX,
    double* y_d,
    const IndexType incY,
    const double c,
    const double s,
    SyncToken* syncToken )
{
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

    cublasDrot( n, x_d, incX, y_d, incY, c, s );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

/** rotm */

template<>
void CUDABLAS1::rotm(
    const IndexType n,
    float* x_d,
    const IndexType incX,
    float* y_d,
    const IndexType incY,
    const float* p_d,
    SyncToken* syncToken )
{
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

    cublasSrotm( n, x_d, incX, y_d, incY, p_d );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
}

template<>
void CUDABLAS1::rotm(
    const IndexType n,
    double* x_d,
    const IndexType incX,
    double* y_d,
    const IndexType incY,
    const double* p_d,
    SyncToken* syncToken )
{
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

    cublasDrotm( n, x_d, incX, y_d, incY, p_d );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }

    cublasSetKernelStream( NULL );
    LAMA_CHECK_CUBLAS_ERROR
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
