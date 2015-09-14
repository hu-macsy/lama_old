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
 * @brief Implementations of BLAS1 routines in CUDA itself, here sum
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Assert.hpp>

#include <scai/common/cuda/CUDAError.hpp>

#include <scai/lama/cuda/CUDABLAS1.hpp>
#include <scai/lama/cuda/launchHelper.hpp>
#include <scai/lama/cuda/utils.cu.h>

#include <cuda_runtime.h>

#include <boost/preprocessor.hpp>

namespace scai
{

namespace lama
{

    /** kernel */

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
    void CUDABLAS1::sum_launcher( const int n, T alpha, const T* x, T beta, const T* y, T* z, cudaStream_t stream )
    {
        SCAI_CHECK_CUDA_ACCESS

        const int blockSize = 256;
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( n, dimBlock.x );

        sum_kernel<<< dimGrid, dimBlock, 0, stream>>> ( n, alpha, x, beta, y, z );
    }

    /* ========================================================================= */
    /*       Template Instantiations                                             */
    /* ========================================================================= */

#define SCAI_CUDA_BLAS1_INSTANTIATE(z, I, _)                 \
                                                             \
template void CUDABLAS1::sum_launcher<ARITHMETIC_CUDA_TYPE_##I>(   \
    const int,                                                     \
    ARITHMETIC_CUDA_TYPE_##I,                                      \
    const ARITHMETIC_CUDA_TYPE_##I*,                               \
    ARITHMETIC_CUDA_TYPE_##I,                                      \
    const ARITHMETIC_CUDA_TYPE_##I*,                               \
    ARITHMETIC_CUDA_TYPE_##I*,                                     \
    cudaStream_t );   

    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT, SCAI_CUDA_BLAS1_INSTANTIATE, _ )

#undef SCAI_CUDA_BLAS1_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
