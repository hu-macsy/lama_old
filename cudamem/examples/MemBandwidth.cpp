
#include <memory/memory.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <logging/logging.hpp>

#include <cudamem/CUDAError.hpp>
#include <common/Walltime.hpp>

#include <iostream>

using namespace memory;

LAMA_LOG_DEF_LOGGER( logger, "MemBandwidth" )

template <typename ValueType> 
void bench( LAMAArray<ValueType>& array )
{
    ContextPtr cudaContext = Context::getContext( context::CUDA );
    ContextPtr hostContext = Context::getContext( context::Host );

    const IndexType N = 1024;
    const IndexType NITER = 16;

    {
        HostWriteOnlyAccess<ValueType> read( array, N );
        for ( IndexType i = 0; i < N; ++i )
        {
            read[i] = 1.0;
        }
    }

    double time = common::Walltime::get();

    for ( int iter = 0; iter < NITER; ++iter )
    {
        // Transfer: Host->CUDA by WriteAccess on CUDA, invalidates Host
        {
            WriteAccess<ValueType> write( array, cudaContext );
        }
        // Transfer: CUDA->Host by WriteAccess on HOST, invalidates CUDA
        {
            HostWriteAccess<ValueType> write( array );
        }
    }

    time = common::Walltime::get() - time ;

    double bytes = N;
    bytes *= sizeof( ValueType );
    bytes *= NITER;
    double mbytes = bytes / ( 1024.0 * 1024.0 );

    std::cout << "Transfer " << mbytes << " MBytes in " << time << " seconds." << std::endl;
}


int main()
{
    ContextPtr cudaContext = Context::getContext( context::CUDA );
    ContextPtr hostContext = Context::getContext( context::Host );

    LAMAArray<float> A1;
    LAMAArray<float> A2( hostContext );

    bench( A1 );
    bench( A2 );
}

