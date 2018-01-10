/**
 * @file tasking/examples/cuda/DemoStream.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief ToDo: Missing description in ./tasking/examples/cuda/DemoStream.cpp
 * @author Thomas Brandes
 * @date 14.03.2016
 */

#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>

#include <iostream>

using namespace scai;
using namespace common;
using namespace tasking;

/* --------------------------------------------------------------------- */

float* myAllocate( int N )
{
    // allocate memory for the current context
    SCAI_CHECK_CUDA_ACCESS
    CUdeviceptr d_pointer = 0;
    size_t size = sizeof( float ) * N;
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &d_pointer, sizeof( float ) * N ), "cuMemAlloc( size = " << size << " ) failed." )
    return reinterpret_cast<float*>( d_pointer );
}

/* --------------------------------------------------------------------- */

void myFree( const float* d_data )
{
    SCAI_CHECK_CUDA_ACCESS
    // free memory for the current context
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( d_data );
    SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << d_data << " ) failed" )
}

/* --------------------------------------------------------------------- */

void myInit( float* d_data, float val, int N )
{
    SCAI_CHECK_CUDA_ACCESS
    // asynchronous execution
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();
    SCAI_ASSERT( syncToken, "No CUDA stream sync token set" )
    cudaStream_t stream = syncToken->getCUDAStream();
    unsigned int* ival = reinterpret_cast<unsigned int*>( &val );
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( d_data );
    SCAI_CUDA_DRV_CALL( cuMemsetD32Async( pointer, *ival, N, stream ), "cuMemsetAsync" )
}

/* --------------------------------------------------------------------- */

void mySum( float* sum, const float d_array[], const int n )
{
    // asynchronous execution
    const CUDACtx& ctx = CUDAAccess::getCurrentCUDACtx();
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();
    SCAI_ASSERT( syncToken, "No CUDA stream sync token set" )
    cudaStream_t stream = syncToken->getCUDAStream();
    cublasHandle_t handle = ctx.getcuBLASHandle();
    std::cout << "Run cublas Sasum asynchronously via stream " << stream << std::endl;
    SCAI_CUBLAS_CALL( cublasSetStream( handle, stream ), "set stream" );
    SCAI_CUBLAS_CALL( cublasSasum( ctx.getcuBLASHandle(), n, d_array, 1, sum ),
                      "cublasSasum for float" );
}

/* --------------------------------------------------------------------- */


int main( int argc, const char** argv )
{
    // using SCAI_DEVICE, SCAI_THREADPOOL_SIZE, ....
    Settings::parseArgs( argc, argv );
    int nr = 0;   // take this as default
    Settings::getEnvironment( nr, "SCAI_DEVICE" );
    const int NSIZE = 16 * 1024 * 1024;   // problem size of one task
    const float VAL = 2.0;
    CUDACtx ctx( nr );
    CUDAAccess access( ctx );
    float* d_data = myAllocate( NSIZE );
    std::cout << "Allocated data on device, d_data = " << d_data << std::endl;
    {
        CUDAStreamSyncToken token( ctx, StreamType::ComputeStream );
        {
            // make sync token available for running compute kernels asynchronously
            SyncToken::ScopedAsynchronous scope( token );
            myInit( d_data, VAL, NSIZE );
        }
        std::cout << "Initialize data on device runs async" << std::endl;
    }
    std::cout << "Initialized data on device done" << std::endl;
    float s = 0.0;
    {
        CUDAStreamSyncToken token( ctx, StreamType::ComputeStream );
        {
            SyncToken::ScopedAsynchronous scope( token );
            mySum( &s, d_data, NSIZE );
        }
    }
    std::cout << "Final result, s= " << s << std::endl;
    myFree( d_data );
    std::cout << "Freed memory" << std::endl;
}
