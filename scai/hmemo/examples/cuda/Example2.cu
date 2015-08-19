
#include <scai/hmemo.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <scai/logging.hpp>

#include <scai/common/cuda/CUDAError.hpp>

#include <iostream>

using namespace hmemo;

SCAI_LOG_DEF_LOGGER( logger, "CudaExample" )

template<typename ValueType>
ValueType sum( const ValueType array[], const IndexType n )
{
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

    ValueType zero = static_cast<ValueType>( 0 );

    ValueType result = thrust::reduce( data, data + n, zero, thrust::plus<ValueType>() );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );

    SCAI_LOG_INFO( logger, "sum of " << n << " values = " << result )

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
}

int main()
{
    std::cout << "try to get " << context::CUDA << " context from factory" << std::endl;
    ContextPtr cudaContext1 = Context::getContextPtr( context::CUDA, 0 );
    std::cout << "cudaContext1 = " << *cudaContext1 << std::endl;

    std::cout << "try to get " << context::CUDA << " context from factory" << std::endl;
    ContextPtr cudaContext2 = Context::getContextPtr( context::CUDA, 1 );
    std::cout << "cudaContext2 = " << *cudaContext2 << std::endl;

    std::cout << "try to get " << context::Host << " context from factory" << std::endl;
    ContextPtr hostContext = Context::getContextPtr( context::Host, 1 );
    std::cout << "hostContext = " << *hostContext << std::endl;

    const IndexType N = 100;

    LAMAArray<double> data;
    
    std::cout << "data = " << data << std::endl;

    {
        SCAI_LOG_INFO( logger, "write only on host" )
        WriteOnlyAccess<double> write( data, hostContext, N );
        double* v = write.get();
        for ( IndexType i = 0; i < N; ++i )
        {
            v[i] = 1.0;
        }
    }

    {
        SCAI_LOG_INFO( logger, "read on gpu 0" )
        ReadAccess<double> read( data, cudaContext1 );
        SCAI_CONTEXT_ACCESS( cudaContext1 )
        double s = sum( read.get(), data.size() );
        std::cout << "sum = " << s << ", should be " << N  << std::endl;
    }

    {
        SCAI_LOG_INFO( logger, "write on gpu1" )
        WriteAccess<double> write( data, cudaContext1 );
        SCAI_CONTEXT_ACCESS( cudaContext1 )
        add( write.get(), data.size() );
    }

    {
        SCAI_LOG_INFO( logger, "write on gpu2" )
        WriteAccess<double> write( data, cudaContext2 );
        SCAI_CONTEXT_ACCESS( cudaContext2 )
        add( write.get(), data.size() );
    }

    {
        SCAI_LOG_INFO( logger, "read on host" )
        ReadAccess<double> read( data );
        const double* values = read.get();
        for ( IndexType i = 0; i < N; ++i )
        {
            SCAI_ASSERT_EQUAL( values[i], 3.0, "wrong value after add" )
        }
    }
}

