/**
 * @file tasking/test/cuda/CUDAKernel.cu
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Basic tests for LAMA arrays with context/memory at CUDA devices
 * @author Thomas Brandes
 * @date 08.07.2015
 */

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/cuda/CUDASettings.hpp>

#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <thrust/device_vector.h>
#include <thrust/fill.h>

using namespace scai;
using namespace tasking;

/* --------------------------------------------------------------------- */

float sum( const float array[], const int n )
{
    thrust::device_ptr<float> data( const_cast<float*>( array ) );
    float zero = static_cast<float>( 0 );
    float result = thrust::reduce( data, data + n, zero, thrust::plus<float>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    return result;
}

__global__
void initKernel( float* out, const int n, const float value )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = value;
    }
}

void init( float array[], const int n, const float value )
{
    SCAI_CHECK_CUDA_ACCESS
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();
    cudaStream_t stream = 0;

    if ( syncToken )
    {
        // asynchronous execution takes other stream and will not synchronize later
        stream = syncToken->getCUDAStream();
    }

    const int blockSize = common::CUDASettings::getBlockSize( n );

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( n, dimBlock.x );

    initKernel <<< dimGrid, dimBlock, 0, stream>>>( array, n, value );
}

