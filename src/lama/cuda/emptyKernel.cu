#include <lama/Walltime.hpp>

__global__
void empty_kernel()
{
}

extern "C" double getKernelLaunchTime( int devNo )
{
    cudaSetDevice( devNo );

    empty_kernel<<<1, 1>>>( );

    double time = lama::Walltime::get();
    cudaDeviceSynchronize();
    empty_kernel<<<1, 1>>>();
    cudaDeviceSynchronize();
    return lama::Walltime::get() - time;
}
