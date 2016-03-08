/**
 * @file common/test/cuda/CUDAKernel.cpp
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
 * @brief Basic tests for LAMA arrays with context/memory at CUDA devices
 * @author: Thomas Brandes
 * @date 08.07.2015
 **/

#include <scai/common/cuda/CUDAError.hpp>

#include <thrust/device_vector.h>
#include <thrust/fill.h>

/* --------------------------------------------------------------------- */

float sum( const float array[], const int n )
{
    thrust::device_ptr<float> data( const_cast<float*>( array ) );

    float zero = static_cast<float>( 0 );

    float result = thrust::reduce( data, data + n, zero, thrust::plus<float>() );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );

    return result;
}

void init( const float array[], const int n, const float value )
{
    thrust::device_ptr<float> data( const_cast<float*>( array ) );
 
    thrust::fill( data, data + n, value );
}

