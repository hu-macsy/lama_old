#include <cuda_runtime_api.h>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

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
    }
}
