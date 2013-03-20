/**
 * @file blas1.cu
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
 * @brief blas1.cu
 * @author lschubert
 * @date 06.07.2012
 * $Id$
 */

#include <lama/cuda/CUDABLAS1.hpp>
#include <lama/cuda/launchHelper.hpp>
#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>

#include <cuda_runtime.h>

const int block_size = 512;

namespace lama
{

/** kernel */

/** assign */
template<typename T1,typename T2>
__global__
void ass_kernel( const int n, T1* dst_d, const T2 value )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if( i < n )
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
    if( i < n )
    {
        z[i] = alpha * x[i] + beta * y[i];
    }
}

/** vamax */
template<typename T>
__global__
void vamax_kernel(
    const int n,
    const int numThreads,
    const int numElemsperThread,
    const T* x_d,
    const int incx,
    T* scratch_d )
{
    //TODO assert(blockDim.x%2 == 0)
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int tId = threadIdx.x;
    int stride = blockDim.x / 2;
    __shared__
    T localMax[block_size];
    if( i < numThreads )
    {
        int idx = i;
        T tmpMax = 0;
        if( idx % incx == 0 )
        {
            tmpMax = abs( x_d[idx] );
        }
        for( int j = 1; j < numElemsperThread; ++j )
        {
            idx += numThreads;
            if( idx >= n )
            {
                break;
            }
            if( idx % incx != 0 )
            {
                continue;
            }
            T value = abs( x_d[idx] );
            if( value > tmpMax )
            {
                tmpMax = value;
            }
        }

        localMax[tId] = tmpMax;
    }
    else
    {
        localMax[tId] = 0;
    }
    __syncthreads();
    while( stride > 0 )
    {
        if( tId < stride && localMax[tId + stride] > localMax[tId] )
        {
            localMax[tId] = localMax[tId + stride];
        }
        stride /= 2;
        __syncthreads();
    }
    if( tId == 0 )
    {
        atomicMax( scratch_d, localMax[tId] );
    }
}

/** max */

template<typename T>
__global__
void max_kernel( const int n, T* scratch_d )
{
    //TODO assert(blockDim.x%2 == 0)
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int tId = threadIdx.x;
    int bId = blockIdx.x;
    int stride = blockDim.x / 2;
    __shared__
    T localMax[block_size];
    if( i < n )
    {
        localMax[tId] = scratch_d[i];
    }
    else
    {
        localMax[tId] = 0;
    }
    __syncthreads();
    while( stride > 0 )
    {
        if( tId < stride && localMax[tId + stride] > localMax[tId] )
        {
            localMax[tId] = localMax[tId + stride];
        }
        stride /= 2;
        __syncthreads();
    }
    if( tId == 0 )
    {
        scratch_d[bId] = localMax[tId];
    }
}

/** gather */
template<typename T1,typename T2,typename T3>
__global__
void gather_kernel( T1* dst_d, const T2* src_d, const T3* indexes_d, const int n )

{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if( i < n )
    {
        dst_d[i] = src_d[indexes_d[i]];
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
    if( threadId < limit )
    {
        sum[threadId] += sum[threadId + stride];
    }
}

/** launcher */

template<typename T>
void CUDABLAS1::ass_launcher( const int n, const T value, T* x, cudaStream_t stream )
{
    LAMA_CHECK_CUDA_ACCESS
    ;

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
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
    ;

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
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
