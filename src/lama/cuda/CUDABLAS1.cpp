/**
 * @file CUDABLAS1.cpp
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
 * @brief CUDABLAS1.cpp
 * @author lschubert
 * @date 05.07.2012
 * $Id$
 */

// hpp
#include <lama/cuda/CUDABLAS1.hpp>

// others
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>

// macros
#include <lama/macros/unused.hpp>

namespace lama
{

/** scale */

template<>
void CUDABLAS1::scal( IndexType n, const float alpha, float* x_d, const IndexType incx, SyncToken* syncToken )
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

    double res = cublasSasum( n, x_d, incX );

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

// instantiation
template void CUDABLAS1::sum<float>(
    const IndexType n,
    float alpha,
    const float* x,
    float beta,
    const float* y,
    float* z,
    SyncToken* syncToken );
template void CUDABLAS1::sum<double>(
    const IndexType n,
    double alpha,
    const double* x,
    double beta,
    const double* y,
    double* z,
    SyncToken* syncToken );

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

/** rotg */
// TODO: implement
/** rotmg */
// TODO: implement
template<typename T>
void CUDABLAS1::ass( const IndexType n, const T value, T* x, SyncToken* syncToken )
{
    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    ass_launcher( n, value, x, stream );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        cudaStreamSynchronize( stream );
        LAMA_CHECK_CUDA_ERROR
    }
}

// instantiation
template void CUDABLAS1::ass<float>( const IndexType n, const float value, float* x, SyncToken* syncToken );
template void CUDABLAS1::ass<double>( const IndexType n, const double value, double* x, SyncToken* syncToken );

template<typename T>
T CUDABLAS1::viamax( const IndexType n, const T* x, const IndexType incx, SyncToken* syncToken )
{
    int maxIdx = iamax( n, x, incx, syncToken );

    float max = -1.0;
    cudaMemcpy( &max, x + maxIdx, sizeof(T), cudaMemcpyDeviceToHost );

    LAMA_CHECK_CUDA_ERROR

    return max;
}

// instantiation
template float CUDABLAS1::viamax<float>(
    const IndexType n,
    const float* x,
    const IndexType incx,
    SyncToken* syncToken );
template double CUDABLAS1::viamax<double>(
    const IndexType n,
    const double* x,
    const IndexType incx,
    SyncToken* syncToken );

} /* namespace lama */
