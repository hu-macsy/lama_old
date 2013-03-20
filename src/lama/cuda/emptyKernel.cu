#include <omp.h>

__global__
void empty_kernel()
{
}

extern "C" double getKernelLaunchTime( int devNo )
{
    cudaSetDevice( devNo );

    empty_kernel<<<1, 1>>>( );

    double time = omp_get_wtime();
    cudaDeviceSynchronize();
    empty_kernel<<<1, 1>>>();
    cudaDeviceSynchronize();
    return omp_get_wtime() - time;
}
