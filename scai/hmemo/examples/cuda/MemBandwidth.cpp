
#include <scai/hmemo.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <scai/logging.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>

using namespace scai;
using namespace scai::hmemo;

SCAI_LOG_DEF_LOGGER( logger, "MemBandwidth" )

template <typename ValueType> 
void bench( LAMAArray<ValueType>& array )
{
    ContextPtr cudaContext = Context::getContextPtr( context::CUDA );
    ContextPtr hostContext = Context::getContextPtr( context::Host );

    const IndexType N = 8 * 1024 * 1024;
    const IndexType NITER = 128;

    {
        WriteOnlyAccess<ValueType> write( array, N );
        ValueType* data = write.get();
        for ( IndexType i = 0; i < N; ++i )
        {
            data[i] = 1.0;
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
            WriteAccess<ValueType> write( array, hostContext );
        }
    }

    time = common::Walltime::get() - time ;

    double bytes = N;
    bytes *= sizeof( ValueType );
    bytes *= NITER;
    double mbytes = bytes / ( 1024.0 * 1024.0 );

    double gBytePerSecond = ( mbytes / 1024.0 ) / time;

    std::cout << "Transfer " << mbytes << " MBytes in " << time << " seconds." << std::endl;
    std::cout << "This is " << gBytePerSecond << " GByte/s (round trip)" << std::endl;
}


int main()
{
    ContextPtr cudaContext = Context::getContextPtr( context::CUDA );
    ContextPtr hostContext = Context::getContextPtr( context::Host );

    LAMAArray<float> A1( hostContext );  // same as LAMAArray<float> A1;
    LAMAArray<float> A2( cudaContext );

    bench( A1 );
    bench( A2 );
}
