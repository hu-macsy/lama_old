/**
 * @file CUDA/main.cu
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
 * @brief Simple CUDA example printing the device properties.
 * @author Jan Ecker
 * @date 20.03.2013
 */

#include <cuda_runtime_api.h>
#include <iostream>

using namespace std;

int main( int argc, char** argv)
{

    int deviceCount = 0;

    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    if ( error_id != cudaSuccess )
    {
        cerr << "cudaGetDeviceCount failed" << endl;
        return -1;
    }

    int driverVersion = 0;
    int runtimeVersion = 0;

    cudaDriverGetVersion( &driverVersion );
    cudaRuntimeGetVersion( &runtimeVersion );

    cout << "CUDA Driver Version = " << driverVersion << endl;
    cout << "CUDA Runtime Version = " << runtimeVersion << endl;

    for (int dev = 0; dev < deviceCount; ++dev) 
    {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties( &deviceProp, dev );
        cout << "Device " << dev << ": " << deviceProp.name << endl;
        cout << "  total amount of global memory = " << deviceProp.totalGlobalMem << " Byte"
             << " = " << deviceProp.totalGlobalMem / ( 1024.0 * 1024.0 ) << " MByte" << endl;
        cout << "  compute capability = " << deviceProp.major << "." << deviceProp.minor << endl;
    }
}
