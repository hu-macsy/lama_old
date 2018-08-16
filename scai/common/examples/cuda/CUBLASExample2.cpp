/**
 * @file examples/cuda/CUBLASExample2.cpp
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
 * @brief ToDo: Missing description in ./examples/cuda/CUBLASExample2.cpp
 * @author Thomas Brandes
 * @date 11.03.2016
 */

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>
#include <iostream>

/* --------------------------------------------------------------------- */

using namespace scai;
using namespace common;

/* --------------------------------------------------------------------- */

float* myAllocate( const float h_data[], int N )
{
    // allocate memory on the accessed device and copy host data to it
    CUdeviceptr pointer = 0;
    size_t size = sizeof( float ) * N;
    // allocate memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, sizeof( float ) * N ), "cuMemAlloc( size = " << size << " ) failed." )
    // transfer host data
    SCAI_CUDA_DRV_CALL( cuMemcpyHtoD( pointer, h_data, size ), "tranfer host->device" )
    return reinterpret_cast<float*>( pointer );
}

/* --------------------------------------------------------------------- */

void myFree( const float* d_data )
{
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( d_data );
    SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << d_data << " ) failed" )
}

/* --------------------------------------------------------------------- */

float myDot( float* d_a, float* d_b, int n )
{
    float dot = 0.0;   // result argument
    const CUDACtx& device = CUDAAccess::getCurrentCUDACtx();
    SCAI_CUBLAS_CALL( cublasSdot( device.getcuBLASHandle(), n, d_a, 1, d_b, 1, &dot ),
                      "cublasSDot for float" );
    return dot;
}

/* --------------------------------------------------------------------- */

int main( int argc, const char** argv )
{
    // at least --SCAI_DEVICE=id may be specified
    Settings::parseArgs( argc, argv );
    int nr = 0;   // take this as default
    Settings::getEnvironment( nr, "SCAI_DEVICE" );
    CUDACtx device( nr );
    float a[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    float b[] = { 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 };
    int n  = sizeof( a ) / sizeof( float );
    int n1 = sizeof( b ) / sizeof( float );
    SCAI_ASSERT_EQUAL( n, n1, "mismatch of arrays for dot product" )
    CUDAAccess access( device );
    float* d_a = myAllocate( a, n );
    float* d_b = myAllocate( b, n1 );
    float dot = myDot( d_a, d_b, n );
    std::cout << "dot product a = [ " ;

    for ( int i = 0; i < n ; ++i )
    {
        std::cout << a[i] << " ";
    }

    std::cout << "] x b = [ " ;

    for ( int i = 0; i < n ; ++i )
    {
        std::cout << b[i] << " ";
    }

    std::cout << "]  = " << dot;
    float r = 0;

    for ( int i = 0; i < n ; ++i )
    {
        r += a[i] * b[i];
    }

    std::cout << ", should be " << r << std::endl;
    // free memory
    myFree( d_a );
    myFree( d_b );
}
