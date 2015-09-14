/**
 * @file utils.cu.h
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
 * @brief utils.cu.h
 * @author jiri
 * @date 15.06.2010
 * @since 1.0.0
 */

#include <math.h>

//defines __global__ __device__ and __shared__ for the Eclipse CDT_Parser
//so we don't get a syntax error in each kernel.
#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__
#define __host__
#define __umul24(a,b)
#endif

/**
 * @brief the maximum grid size in one dimension.
 */
const unsigned int lama_maxGridSize_cuda = 65535;

/**
 * @brief makeGrid creates a grid large enough to start the passed number of
 *        threads with the given blockSize.
 *
 * @param[in] numThreads the number of threads that are needed.
 * @param[in] blockSize  the total number of threads in a block.
 * @return               a 1D or 2D grid that will start enough blocks to have
 *                       at least numThreads running.
 */
inline dim3 makeGrid( const unsigned int numThreads, const unsigned int blockSize )
{
    const unsigned int numBlocks = ( numThreads + blockSize - 1 ) / blockSize;

    if( numBlocks <= lama_maxGridSize_cuda )
    {
        //fits in a 1D grid
        return dim3( numBlocks );
    }
    else
    {
        //2D grid is required
        const unsigned int side = (unsigned int) ceil( sqrt( (double) numBlocks ) );
        return dim3( side, side );
    }
}

/**
 * @brief threadId calculates the global one dimensional x thread id from the
 *        passed parameters.
 *
 * threadId calculates the global one dimensional x thread id from the passed
 * parameters. threadId is a device function that can only be called from a CUDA
 * kernel.
 *
 * @param[in] gridDim   the size of the grid.
 * @param[in] blockIdx  the id of the block
 * @param[in] blockDim  the size of a block.
 * @param[in] threadIdx the local id of the thread.
 * @return              a global one dimensional x thread id.
 */
__inline__ __device__
unsigned int threadId( const dim3 gridDim, const dim3 blockIdx, const dim3 blockDim, const dim3 threadIdx )
{
    return __umul24( blockDim.x, blockIdx.
                    x + __umul24( blockIdx.y, gridDim.x ) ) +threadIdx.x;
}

/**
 * @brief localThreadId3D calculates the local three-dimensional thread id from
 *        the passed parameters.
 *
 * localThreadId3D calculates the local three-dimensional thread id from the
 * passed parameters. localThreadId3D is a device function that can only be
 * called from a CUDA kernel.
 *
 * @param[in] blockDim  the size of a block.
 * @param[in] threadIdx the local id of the thread.
 * @return              a local three-dimensional thread id.
 */
__inline__ __device__
unsigned int localThreadId3D( const dim3 blockDim, const dim3 threadIdx )
{
    return __umul24( blockDim.y, __umul24( blockDim.x, threadIdx.z ) )
    + __umul24( blockDim.x, threadIdx.y )
    +threadIdx.x;
}

/** * @brief blockId calculates the one dimensional block id of a thread in a two
 *        dimensional grid from the passed parameters.
 *
 * blockId calculates the two dimensional block id of a thread in a two
 * dimensional grid from the passed parameters. It is working for one
 * dimensional blocks also, but not save for work with three dimensional blocks.
 * blockId is a device function that can only be called from a CUDA kernel.
 *
 * @param[in] gridDim   the size of the grid.
 * @param[in] blockIdx  the id of the block.
 * @return              the two dimensional id of the block.
 */
__inline__ __device__
unsigned int blockId( const dim3 gridDim, const dim3 blockIdx )
{
    return __umul24( gridDim.x, blockIdx.y ) +blockIdx.x;
}

/**
 * @brief returns the halve of val without any rest.
 *
 * returns the halve of val without any rest. For example if val was 4, 2 is
 * returned, if val was 5, 2 is returned, too. div2 is a device function that
 * can only be called from a CUDA kernel.
 *
 * @param[in] val   the value to halve.
 *
 * @return the halve of val without any rest.
 */
__inline__ __device__
unsigned int div2( const int val )
{
    return val >> 1;
}

/**
 * @brief returns the rest of val/2.
 *
 * returns the rest of val/2. For example if val was 4, 0 is returned, if val
 * is 5, 1 is returned. mod2 is a device function that can only be called from
 * a CUDA kernel.
 *
 * @param[in] val   the value to return the rest of val/2 from.
 *
 * @return the rest of val/2.
 */
__inline__ __device__
unsigned int mod2( const int val )
{
    return val & 1;
}

/**
 * @brief returns the rounded up halve of val.
 *
 * returns the rounded up halve of val. For example if val was 4, 2 is
 * returned, if val was 5, 3 is returned. halve is a device function that can
 * only be called from a CUDA kernel.
 *
 * @param[in] val   the value to return the rounded up halve from.
 *
 * @return the rounded up halve of val.
 */
__inline__ __device__
unsigned int halve( const int val )
{
    return div2( val ) + mod2( val );
}

template<typename ValueType>
__inline__       __device__ ValueType lamaDeviceFMA(
    const ValueType alpha,
    const ValueType v1,
    const ValueType beta,
    const ValueType v2 )
{
    return alpha * v1 + beta * v2;
}

template<typename ValueType,int beta>
__inline__       __device__ ValueType lamaDeviceFMA( const ValueType alpha, const ValueType v1, const ValueType v2 )
{
    return alpha * v1 + beta * v2;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float,1>( const float alpha, const float v1, const float v2 )
{
    return alpha * v1 + v2;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float, -1>( const float alpha, const float v1, const float v2 )
{
    return alpha * v1 - v2;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float,0>( const float alpha, const float v1, const float v2 )
{
    return alpha * v1;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double,1>( const double alpha, const double v1, const double v2 )
{
    return alpha * v1 + v2;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double, -1>( const double alpha, const double v1, const double v2 )
{
    return alpha * v1 - v2;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double,0>( const double alpha, const double v1, const double v2 )
{
    return alpha * v1;
}

template<typename ValueType,int alpha,int beta>
__inline__       __device__ ValueType lamaDeviceFMA( const ValueType v1, const ValueType v2 )
{
    return alpha * v1 + beta * v2;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float,1,0>( const float v1, const float )
{
    return v1;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float,1,1>( const float v1, const float v2 )
{
    return v1 + v2;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float,1, -1>( const float v1, const float v2 )
{
    return v1 - v2;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double,1,0>( const double v1, const double )
{
    return v1;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double,1,1>( const double v1, const double v2 )
{
    return v1 + v2;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double,1, -1>( const double v1, const double v2 )
{
    return v1 - v2;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float, -1,0>( const float v1, const float )
{
    return -v1;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float, -1,1>( const float v1, const float v2 )
{
    return v2 - v1;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float, -1, -1>( const float v1, const float v2 )
{
    return -v1 - v2;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double, -1,0>( const double v1, const double )
{
    return -v1;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double, -1,1>( const double v1, const double v2 )
{
    return v2 - v1;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double, -1, -1>( const double v1, const double v2 )
{
    return -v1 - v2;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float,0,0>( const float, const float )
{
    return 0.0f;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float,0,1>( const float, const float v2 )
{
    return v2;
}

template<>
__inline__ __device__
float lamaDeviceFMA<float,0, -1>( const float, const float v2 )
{
    return -v2;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double,0,0>( const double, const double )
{
    return 0.0;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double,0,1>( const double, const double v2 )
{
    return v2;
}

template<>
__inline__ __device__
double lamaDeviceFMA<double,0, -1>( const double, const double v2 )
{
    return -v2;
}
