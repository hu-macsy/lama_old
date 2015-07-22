
#include <memory/memory.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <logging/logging.hpp>

#include <cudamem/CUDAError.hpp>

#include <iostream>
#include <unistd.h>

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
void add( ValueType* array, const IndexType n )
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
    ContextPtr cudaContext = Context::getContextPtr( context::CUDA );
    COMMON_ASSERT( cudaContext, "NULL context" )
    std::cout << "cudaContext = " << *cudaContext << std::endl;

    MemoryPtr cudaMemory = cudaContext->getMemoryPtr();
    COMMON_ASSERT( cudaMemory, "NULL memory" )
    std::cout << "cudaMemory = " << *cudaMemory << std::endl;

    std::cout << "try to get " << context::Host << " context from factory" << std::endl;
    ContextPtr hostContext = Context::getContextPtr( context::Host );
    COMMON_ASSERT( hostContext, "NULL context" )
    std::cout << "hostContext = " << *hostContext << std::endl;

    MemoryPtr hostMemory = hostContext->getMemoryPtr();
    COMMON_ASSERT( hostMemory, "NULL memory" )
    std::cout << "hostMemory = " << *hostMemory << std::endl;

    MemoryPtr cudaHostMemory = cudaContext->getHostMemoryPtr();
    COMMON_ASSERT( cudaHostMemory, "NULL memory" )
    std::cout << "cudaHostMemory = " << *cudaHostMemory << std::endl;

    const IndexType N = 100;

    LAMAArray<double> data( cudaContext );
    
    std::cout << "data = " << data << std::endl;

    {
        LAMA_LOG_INFO( logger, "write only on cuda host" )
        // HostWriteOnlyAccess<double> write( data,  N );

        WriteOnlyAccess<double> writeData( data, N );

        double* dataHost = writeData;

        for ( IndexType i = 0; i < N; ++i )
        {
            writeData[i] = 1.0;
        }
    }

    std::cout << "After host write: data = " << data << std::endl;

    {
        LAMA_LOG_INFO( logger, "read on cuda" )
        ReadAccess<double> read( data, cudaContext );
        LAMA_CONTEXT_ACCESS( cudaContext )
        double s = sum( read.get(), data.size() );
        std::cout << "sum = " << s << ", should be " << N  << std::endl;
    }

    std::cout << "After cuda read: data = " << data << std::endl;

    {
        LAMA_LOG_INFO( logger, "write on cuda" )
        WriteAccess<double> write( data, cudaContext );
        LAMA_CONTEXT_ACCESS( cudaContext )
        add( static_cast<double*>( write ), data.size() );
    }

    std::cout << "After cuda write: data = " << data << std::endl;

    {
        LAMA_LOG_INFO( logger, "read on host" )
        ReadAccess<double> read( data );
        sleep( 1 );
        for ( IndexType i = 0; i < N; ++i )
        {
            COMMON_ASSERT_EQUAL( read[i], 2 * 1.0, "wrong value after add, i = " << i )
        }
    }

    std::cout << "After host read: data = " << data << std::endl;
}

