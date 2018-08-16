/**
 * @file examples/cuda/CUDAExample.cu
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
 * @brief ToDo: Missing description in ./examples/cuda/CUDAExample.cu
 * @author Thomas Brandes
 * @date 09.03.2016
 */

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>
#include <iostream>


#include <thrust/device_vector.h>
#include <thrust/fill.h>

/* --------------------------------------------------------------------- */

float sum( const float array[], const int n )
{
    thrust::device_ptr<float> data( const_cast<float*>( array ) );
    float zero = static_cast<float>( 0 );
    float result = thrust::reduce( data, data + n, zero, thrust::plus<float>() );
    return result;
}

void init( const float array[], const int n, const float value )
{
    thrust::device_ptr<float> data( const_cast<float*>( array ) );
    thrust::fill( data, data + n, value );
}

/* --------------------------------------------------------------------- */

using namespace scai;
using namespace common;

int main( int argc, const char** argv )
{
    // at least --SCAI_DEVICE=id may be specified
    Settings::parseArgs( argc, argv );
    int nr = 0;   // take this as default
    Settings::getEnvironment( nr, "SCAI_DEVICE" );
    CUDACtx device( nr );
    std::cout << "CUDA device available." << std::endl;
    CUDAAccess access( device );
    std::cout << "CUDA device accessed." << std::endl;
    // allocate memory
    const int N = 119;
    CUdeviceptr pointer = 0;
    size_t size = sizeof( float ) * N;
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, size ), "cuMemAlloc( size = " << size << " ) failed." )
    std::cout << "CUDA memory allocated." << std::endl;
    float* fpointer = reinterpret_cast<float*>( pointer );
    init( fpointer, N, 3.0 );
    std::cout << "CUDA memory initialzed." << std::endl;
    float s = sum( fpointer, N );
    std::cout << "Computed sum = " << s << ", should be " << N * 3.0 << std::endl;
    // free memory
    SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << pointer << " ) failed" )
}
