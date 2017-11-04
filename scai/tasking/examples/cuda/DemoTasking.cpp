/**
 * @file tasking/examples/cuda/DemoTasking.cpp
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
 * @brief ToDo: Missing description in ./tasking/examples/cuda/DemoTasking.cpp
 * @author Thomas Brandes
 * @date 14.03.2016
 */

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>

#include <iostream>
#include <functional>

using namespace scai;
using namespace common;
using namespace tasking;

using std::bind;
using std::cref;

/* --------------------------------------------------------------------- */

float* myAllocate( const CUDACtx& ctx, int N )
{
    CUDAAccess access( ctx );
    // allocate memory on the accessed device and copy host data to it
    CUdeviceptr d_pointer = 0;
    size_t size = sizeof( float ) * N;
    // allocate memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &d_pointer, sizeof( float ) * N ), "cuMemAlloc( size = " << size << " ) failed." )
    return reinterpret_cast<float*>( d_pointer );
}

/* --------------------------------------------------------------------- */

void myFree( const CUDACtx& ctx, const float* d_data )
{
    CUDAAccess access( ctx );
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( d_data );
    SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << d_data << " ) failed" )
}

/* --------------------------------------------------------------------- */

void myInit( const CUDACtx& ctx, float* d_data, float val, int N )
{
    CUDAAccess access( ctx );
    unsigned int* ival = reinterpret_cast<unsigned int*>( &val );
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( d_data );
    SCAI_CUDA_DRV_CALL( cuMemsetD32( pointer, *ival, N ), "cuMemset" )
}

/* --------------------------------------------------------------------- */

float mySum( const CUDACtx& ctx, const float d_array[], const int n )
{
    CUDAAccess access( ctx );
    float sum;
    SCAI_CUBLAS_CALL( cublasSasum( ctx.getcuBLASHandle(), n, d_array, 1, &sum ),
                      "cublasSasum for float" );
    return sum;
}

/* --------------------------------------------------------------------- */

int main( int argc, const char** argv )
{
    // using SCAI_DEVICE, SCAI_THREADPOOL_SIZE, ....
    Settings::parseArgs( argc, argv );
    int nr = 0;   // take this as default
    Settings::getEnvironment( nr, "SCAI_DEVICE" );
    CUDACtx ctx( nr );
    const int NSIZE = 256 * 1024;   // problem size of one task
    const float VAL = 2.0;
    float* d_data = myAllocate( ctx, NSIZE );
    // Let other thread do the initialization
    {
        TaskSyncToken( bind( &myInit, std::cref( ctx ), d_data, VAL, NSIZE ) );
    }
    float s = mySum( ctx, d_data, NSIZE );
    // Let other thread do the free
    {
        TaskSyncToken( bind( &myFree, std::cref( ctx ), d_data ) );
    }
    std::cout << "Ready task " << nr << " on device " << nr;
    std::cout << ", result = " << s << ", should be " << NSIZE* VAL << std::endl;
}

/* --------------------------------------------------------------------- */
