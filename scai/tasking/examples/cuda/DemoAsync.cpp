/**
 * @file tasking/examples/cuda/DemoAsync.cpp
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
 * @brief ToDo: Missing description in ./tasking/examples/cuda/DemoAsync.cpp
 * @author Thomas Brandes
 * @date 14.03.2016
 */

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/thread.hpp>

#include <iostream>
#include <functional>

using namespace scai;
using namespace common;
using namespace tasking;

/* --------------------------------------------------------------------- */

float* myAllocate( int N )
{
    // allocate memory on the accessed device and copy host data to it
    CUdeviceptr d_pointer = 0;
    size_t size = sizeof( float ) * N;
    // allocate memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &d_pointer, sizeof( float ) * N ), "cuMemAlloc( size = " << size << " ) failed." )
    return reinterpret_cast<float*>( d_pointer );
}

/* --------------------------------------------------------------------- */

void myFree( const float* d_data )
{
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( d_data );
    SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << d_data << " ) failed" )
}

/* --------------------------------------------------------------------- */

void myInit( float* d_data, float val, int N )
{
    unsigned int* ival = reinterpret_cast<unsigned int*>( &val );
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( d_data );
    SCAI_CUDA_DRV_CALL( cuMemsetD32( pointer, *ival, N ), "cuMemset" )
}

/* --------------------------------------------------------------------- */

float mySum( const float d_array[], const int n )
{
    const CUDACtx& ctx = CUDAAccess::getCurrentCUDACtx();
    float sum;
    SCAI_CUBLAS_CALL( cublasSasum( ctx.getcuBLASHandle(), n, d_array, 1, &sum ),
                      "cublasSasum for float" );
    return sum;
}

/* --------------------------------------------------------------------- */

void task( int nr, int device, int N )
{
    std::cout << *thread::getCurrentThreadName() << ": running task " << nr << " on device " << device << std::endl;
    const float VAL = 2.0;
    CUDACtx ctx( device );
    CUDAAccess access( ctx );
    float* d_data = myAllocate( N );
    myInit( d_data, VAL, N );
    float s = mySum( d_data, N );
    myFree( d_data );
    std::cout << "Ready task " << nr << " on device " << device;
    std::cout << ", result = " << s << ", should be " << N* VAL << std::endl;
}

/* --------------------------------------------------------------------- */

int main( int argc, const char** argv )
{
    // using SCAI_DEVICE, SCAI_THREADPOOL_SIZE, ....
    Settings::parseArgs( argc, argv );
    int nr = 0;   // take this as default
    Settings::getEnvironment( nr, "SCAI_DEVICE" );
    // run NT tasks where each thread does full task
    const int NT    = 10;           // number of tasks
    const int NSIZE = 256 * 1024;   // problem size of one task
    std::cout << "Running " << NT << " tasks in thread pool" << std::endl;
    SyncToken* tokenArray[NT];

    for ( int i = 0; i < NT; ++i )
    {
        tokenArray[i] = new TaskSyncToken( std::bind( &task, i, nr, NSIZE ) );
    }

    for ( int i = 0; i < NT; ++i )
    {
        delete tokenArray[i];
    }
}
