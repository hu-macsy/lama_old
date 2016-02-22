/**
 * @file launchHelper.hpp
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
 * @brief launchHelper.hpp
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */

#pragma once

#include <cuda_runtime.h>

namespace scai
{

//defines some CUDA specific things for the Eclipse CDT_Parser
//so we don't get syntax errors.
#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__
#define __host__
#define __syncthreads()
#define __threadfence_block()
dim3 blockIdx;
dim3 blockDim;
dim3 threadIdx;
dim3 gridDim;
int warpSize;
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
    // launch a kernel with numBlocks == 0 results in runtime error

    const unsigned int numBlocks = numThreads > 0 ? ( numThreads + blockSize - 1 ) / blockSize : 1;

    // will be safe as there is always check for legal value 0 <= ithread < numThreads

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

} /* end namespace scai */
