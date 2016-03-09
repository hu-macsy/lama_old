
#include <scai/common/cuda/CUDADevice.hpp>
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

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );

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

    CUDADevice device;

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
