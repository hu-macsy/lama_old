
#include <memory.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <logging.hpp>

#include <common/cuda/CUDAError.hpp>
#include <common/Walltime.hpp>

#include <iostream>

using namespace memory;

LAMA_LOG_DEF_LOGGER( logger, "CudaExample" )

template<typename ValueType>
ValueType sum( const ValueType array[], const IndexType n )
{
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

    ValueType zero = static_cast<ValueType>( 0 );

    ValueType result = thrust::reduce( data, data + n, zero, thrust::plus<ValueType>() );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );

    LAMA_LOG_INFO( logger, "sum of " << n << " values = " << result )

    return result;
}

template<typename ValueType>
__global__
void add_kernel( ValueType* array, IndexType n )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    ValueType one = 1;

    if ( i < n )
    {
        array[i] += one;
    }
}

template<typename ValueType>
void add( ValueType array[], const IndexType n )
{
    const int blockSize = 256;
    const int nblocks   = ( n + blockSize - 1 ) / blockSize;

    dim3 block( blockSize, 1, 1 );
    dim3 grid( nblocks, 1, 1 );

    add_kernel<<<grid, block>>>( array, n );
 
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cuda failure" );
}

template<typename ValueType>
void addHost( ValueType array[], const IndexType n )
{
    for ( IndexType i = 0; i < n; ++i )
    {
       array[i] += 1;
    }
}

template<typename ValueType>
void doBench( LAMAArray<ValueType>& array, const IndexType N )
{
    ContextPtr hostContext = Context::getContextPtr( context::Host );
    ContextPtr cudaContext = Context::getContextPtr( context::CUDA );

    int nhost = 1;
    int ncuda = 1;
    int niter = 50;

    double time = common::Walltime::get();

    // init on host

    {
        WriteOnlyAccess<double> write( array, hostContext, N );
        double* v = write.get();
        for ( IndexType i = 0; i < N; ++i )
        {
            v[i] = 0.0;
        }
    }

    for ( int iter = 0; iter < niter; ++iter )
    {
        // do some work on cuda

        for ( int k = 0; k < ncuda; ++k )
        {
            WriteAccess<double> write( array, cudaContext );
            LAMA_CONTEXT_ACCESS( cudaContext )
            add( write.get(), N );
        }

        // do some work on host

        for ( int k = 0; k < nhost; ++k )
        {
            WriteAccess<double> write( array, hostContext );
            addHost( write.get(), N );
        }
    }

    // compute result

    double res = 0.0;

    {
        ReadAccess<double> read( array, cudaContext );
        LAMA_CONTEXT_ACCESS( cudaContext )
        res = sum( read.get(), N );
    }

    double resExpected = N;
    resExpected *= double ( niter * ( ncuda + nhost ) );

    COMMON_ASSERT_EQUAL( res, resExpected, "wrong result, N = " << N 
        << ", niter = " << niter << ", ncuda = " << ncuda << ", nhost = " << nhost )

    time = common::Walltime::get() - time;

    std::cout << "Time = " << time << " seconds" << std::endl;
}

int main()
{
    const IndexType N = 8 * 1024 * 1024;  // 8 MB data

    ContextPtr hostContextPtr = Context::getContextPtr( context::Host );
    ContextPtr cudaContextPtr = Context::getContextPtr( context::CUDA );

    // First touch on host memory, never uses CUDA host memory

    std::cout << "Benchmark for array, first touch on host memory" << std::endl;

    LAMAArray<double> data1( hostContextPtr->getMemoryPtr() );
    doBench( data1, N );

    std::cout << "Benchmark for array, first touch on cuda memory" << std::endl;
 
    LAMAArray<double> data2( cudaContextPtr->getMemoryPtr() );
    doBench( data2, N );

    std::cout << "Benchmark for array, first touch on cuda host memory" << std::endl;
 
    LAMAArray<double> data3( cudaContextPtr->getHostMemoryPtr() );
    doBench( data3, N );
}

