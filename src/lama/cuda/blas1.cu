/**
 * @file blas1.cu
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
 * @brief blas1.cu
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */

#include <lama/cuda/CUDABLAS1.hpp>
#include <lama/cuda/launchHelper.hpp>
#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>

#include <cuda_runtime.h>

namespace lama
{

/** kernel */

/** assign */
template<typename T1,typename T2>
__global__
void ass_kernel( const int n, T1* dst_d, const T2 value )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( i < n )
    {
        dst_d[i] = value;
    }
}

/** sum */
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

/**
 * @brief Adds the stride to the value of the given array at the position of
 *        threadId, if the threadId is smaller than the given limit.
 * @param[in,out] sum   The array holding the values to be added.
 * @param[in] threadId  The id of the thread.
 * @param[in] limit     The first threadId not to add.
 * @param[in] stride    The stride to add.
 */
template<typename T>
__inline__ __device__
void limitedIntervalAdd( T* sum, const unsigned int threadId, const unsigned int limit, const unsigned int stride )
{
    if ( threadId < limit )
    {
        sum[threadId] += sum[threadId + stride];
    }
}

/** launcher */

template<typename T>
void CUDABLAS1::ass_launcher( const int n, const T value, T* x, cudaStream_t stream )
{
    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    ass_kernel<<< dimGrid, dimBlock, 0, stream>>> ( n, x, value );
}

// instantiation
template void CUDABLAS1::ass_launcher<float>( const IndexType n, const float value, float* x, cudaStream_t stream );
template void CUDABLAS1::ass_launcher<double>( const IndexType n, const double value, double* x, cudaStream_t stream );

template<typename T>
void CUDABLAS1::sum_launcher( const int n, T alpha, const T* x, T beta, const T* y, T* z, cudaStream_t stream )
{
    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    sum_kernel<<< dimGrid, dimBlock, 0, stream>>> ( n, alpha, x, beta, y, z );
}

// instantiation
template void CUDABLAS1::sum_launcher<float>(
    const int n,
    float alpha,
    const float* x,
    float beta,
    const float* y,
    float* z,
    cudaStream_t stream );
template void CUDABLAS1::sum_launcher<double>(
    const int n,
    double alpha,
    const double* x,
    double beta,
    const double* y,
    double* z,
    cudaStream_t stream );

} /* namespace lama */
