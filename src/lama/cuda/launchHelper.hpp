/**
 * @file launchHelper.hpp
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
 * @brief launchHelper.hpp
 * @author lschubert
 * @date 06.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_LAUNCHHELPER_HPP_
#define LAMA_LAUNCHHELPER_HPP_

#include <cuda_runtime.h>

namespace lama
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
#define min(a,b) a
#define max(a,b) a
#endif

#define MIN(a,b) ((a>b)?b:a)
#define MAX(a,b) ((a<b)?b:a)

} /* namespace lama */

#endif // LAMA_LAUNCHHELPER_HPP_
