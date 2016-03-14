
#include <scai/common/cuda/CUDADevice.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>
#include <iostream>

/* --------------------------------------------------------------------- */

using namespace scai;
using namespace common;

/* --------------------------------------------------------------------- */

float* myAllocate( float hostdata[], int N )
{
    // allocate memory on the accessed device and copy host data to it

    CUdeviceptr pointer = 0;

    size_t size = sizeof( float ) * N;

    // allocate memory
    SCAI_CUDA_DRV_CALL( cuMemAlloc( &pointer, sizeof( float ) * N ), "cuMemAlloc( size = " << size << " ) failed." )

    // transfer host data
    SCAI_CUDA_DRV_CALL( cuMemcpyHtoD( pointer, hostdata, size ), "tranfer host->device" )

    return reinterpret_cast<float*>( pointer );
}

/* --------------------------------------------------------------------- */

void myFree( const float* data )
{
    CUdeviceptr pointer = reinterpret_cast<CUdeviceptr>( data );

    SCAI_CUDA_DRV_CALL( cuMemFree( pointer ), "cuMemFree( " << data << " ) failed" )
}

/* --------------------------------------------------------------------- */

float myDot( float* d_a, float* d_b, int n )
{
    float dot = 0.0;   // result argument

    const CUDADevice& device = CUDAAccess::getCurrentCUDADevice();

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

    CUDADevice device( nr );

    float a[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    float b[] = { 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 };

    int n  = sizeof( a ) / sizeof( float );
    int n1 = sizeof( b ) / sizeof( float );

    SCAI_ASSERT_EQUAL( n, n1, "mismatch of arrays for dot product" )

    CUDAAccess access( device );  

    float* d_a = myAllocate( a, n );
    float* d_b = myAllocate( b, n1 );

    float dot = myDot( d_a, d_b, n );

    std::cout << "Result = " << dot << std::endl;

    // free memory

    myFree( d_a );
    myFree( d_b );
}
