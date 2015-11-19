/**
 * @file Allocate.cpp
 *
 * @brief Benchmark allocate on Host vs CUDAHost memory
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#include <scai/hmemo.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/shared_ptr.hpp>

#include <iostream>

using namespace scai;
using namespace scai::hmemo;
using namespace scai::common;

void doit( int NITER, IndexType NSIZE, ContextPtr context )
{
    std::vector<shared_ptr<LAMAArray<double> > > stack;

    ContextPtr host = Context::getHostPtr();

    for ( int iter = 0; iter < NITER; ++iter )
    {
        shared_ptr<LAMAArray<double> > X;
        X.reset( new LAMAArray<double>( context ) );

        // first touch on context, but allocate it on host

        WriteOnlyAccess<double> wX( *X, host, NSIZE );
        stack.push_back( X );
    }
}

/*---------------------------------------------------------------------------*
 * Main program                                                              *
 *---------------------------------------------------------------------------*/

int main( int, char** )
{
    using namespace std;

    ContextPtr cudaContext = Context::getContextPtr( common::context::CUDA );
    ContextPtr hostContext = Context::getContextPtr( common::context::Host );

    static int ITER_VEC[]    = { 10000, 10000, 10000, 3000, 2000,  1000,   700,    500,    300,    200 };
    static IndexType N_VEC[] = {     1,    10,   100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000 };

    int NCASES = sizeof( ITER_VEC ) / sizeof( int );

    for ( int k = 0; k < NCASES; ++k )
    {
        int NITER    = ITER_VEC[k];
        IndexType N  = N_VEC[k];

        double t = common::Walltime::get();

        doit( NITER, N, hostContext );

        double hostTime = common::Walltime::get() - t;

        t = common::Walltime::get();

        doit( NITER, N, cudaContext );

        double cudaTime = common::Walltime::get() - t;

        hostTime *= 1000.0 * 1000.0  / NITER;
        cudaTime *= 1000.0 * 1000.0  / NITER;

        cout << "NITER = " << NITER << ", N = " << N << endl;
        cout << "Time for one allocate on HostMemory: " << hostTime << " us" << endl;
        cout << "Time for one allocate on CUDAHostMemory: " << cudaTime << " us" << endl;
        cout << "Ratio = " << ( cudaTime / hostTime ) << endl;
    }
}
