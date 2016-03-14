
#include <scai/common/Walltime.hpp>

#include <scai/common/cuda/CUDADevice.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <scai/common/Settings.hpp>

#include <scai/tasking/cuda/CUDAStreamPool.hpp>
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <iostream>

using namespace scai;
using namespace tasking;
using namespace common;

int main( int argc, const char** argv )
{
    // at least --SCAI_DEVICE=id may be specified

    Settings::parseArgs( argc, argv );

    const int N_USES = 100000;   // number of stream uses

    int deviceNr = 0;

    Settings::getEnvironment( deviceNr, "SCAI_DEVICE" );

    CUDADevice device( deviceNr );  

    CUDAAccess access( device );

    double t0 = Walltime::get();

    for ( int i = 0; i < N_USES; ++i )
    {
        CUstream stream;
        int flags = 0;
        SCAI_CUDA_DRV_CALL( cuStreamCreate( &stream, flags ), "cuStreamCreate failed" )
        SCAI_CUDA_DRV_CALL( cuStreamDestroy( stream ), "cuStreamDestroy failed" )
    }

    double t1 =  Walltime::get() - t0;

    t0 =  Walltime::get();

    CUDAStreamPool& pool = CUDAStreamPool::getPool( device );

    for ( int i = 0; i < N_USES; ++i )
    {
        CUstream str = pool.reserveStream( CUDAStreamSyncToken::ComputeStream );
        pool.releaseStream( str );
    }

    CUDAStreamPool::freePool( device );

    double t2 =  Walltime::get() - t0;

    t0 =  Walltime::get();

    for ( int i = 0; i < N_USES; ++i )
    {
        CUDAStreamSyncToken( device, CUDAStreamSyncToken::ComputeStream );
    }

    double t3 =  Walltime::get() - t0;

    std::cout << "Measure time for " << N_USES << " stream accesses." << std::endl;
    std::cout << "Time create/destroy stream of CUDA: " << ( t1 / N_USES ) * 1000.0 * 1000.0 << " µs" << std::endl;
    std::cout << "Time get/release stream of pool:    " <<  (t2 / N_USES ) * 1000.0 * 1000.0 << " µs" << std::endl;
    std::cout << "Time construct/destruct SyncToken:  " <<  (t3 / N_USES ) * 1000.0 * 1000.0 << " µs" << std::endl;
}
