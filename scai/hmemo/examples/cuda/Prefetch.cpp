/**
 * @file Prefetch.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Demo of benefits for prefetch
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#include <scai/hmemo.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>

using namespace scai;
using namespace scai::hmemo;

SCAI_LOG_DEF_LOGGER( logger, "Prefetch" )

/*---------------------------------------------------------------------------*
 * workload: subroutine that gives CPU some work to do                       *
 *---------------------------------------------------------------------------*/

void workload( double& A, int NITER )
{
    static const int N = 1024;
    double v = 1.0;

    for ( int iter = 0; iter < NITER; iter++ )
    {
        for ( int i = 0; i < N; ++i )
        {
            double x = v + 1.0 / static_cast<double>( i );
            A = A / x ;
        }

        v -= 1.0 / static_cast<double>( 2 * NITER );
    }
}

/*---------------------------------------------------------------------------*
 * setval: Initialize HArray on Host with a given value                   *
 *---------------------------------------------------------------------------*/

void setval( HArray<double>& array, double val, IndexType N )
{
    ContextPtr hostContext = Context::getContextPtr( common::ContextType::Host );
    WriteOnlyAccess<double> writeArray( array, hostContext, N );

    for ( IndexType i = 0; i < N; ++i )
    {
        writeArray[i] = val;
    }
}

/*---------------------------------------------------------------------------*
 * Main program                                                              *
 *---------------------------------------------------------------------------*/

int main( int argc, char** )
{
    using namespace std;
    double dummy = 1.0;          // for workload
    int    NWORK = 1024;
    int    N     = 4 * 1024 * 1024;  // array size
    // use pinned Host memory in case of at least one argument
    bool isPinned = argc > 1;
    ContextPtr cudaContext = Context::getContextPtr( common::ContextType::CUDA );
    ContextPtr hostContext = Context::getContextPtr( common::ContextType::Host );
    ContextPtr homeContext = isPinned ? cudaContext : hostContext;
    cudaContext->enableZeroCopy( true );
    // HArray<double> A( homeContext );  // uses pinned Host memory if first touch is on CUDA
    HArray<double> A( homeContext );  // uses pinned Host memory if first touch is on CUDA
    // A little bit warm up to initialize every thing
    setval( A, 0.0, N );   // initialize on host
    {
        WriteAccess<double> rA( A, cudaContext );
    }
    setval( A, 1.0, N );   // initialize on host
    cout << "LAMA array A: " << A << endl;
    // transferToGPU: time to transfer content of A from CPU to GPU
    double t = common::Walltime::get();
    {
        // valid data on Host: transfer from Host -> GPU
        WriteAccess<double> wA( A, cudaContext );
    }
    double transferToGPU = common::Walltime::get() - t;
    cout << "Transfer time to GPU: " << transferToGPU << endl;
    // transferFromGPU: time to transfer content of A from GPU to CPU
    t = common::Walltime::get();
    {
        // valid data on GPU: transfer from GPU -> Host
        WriteAccess<double> wA( A, hostContext );
    }
    double transferFromGPU = common::Walltime::get() - t;
    cout << "Transfer time from GPU: " << transferFromGPU << endl;
    // worktime: time to run a certain workload on CPU
    t = common::Walltime::get();
    workload( dummy, NWORK );
    double worktime = common::Walltime::get() - t;
    cout << "Workload time on Host " << worktime << endl;
    setval( A, 0.5, N );
    // prefetchTo: overlap transfer Host->GPU with workload
    t = common::Walltime::get();
    A.prefetch( cudaContext );
    workload( dummy, NWORK );
    {
        WriteAccess<double> rA( A, cudaContext );
    }
    double prefetchTo = common::Walltime::get() - t;
    cout << "Time for transfer to GPU and workload on Host " << prefetchTo << endl;
    // prefetchFrom: overlap transfer GPU->Host with workload
    t = common::Walltime::get();
    A.prefetch( hostContext );
    workload( dummy, NWORK );
    {
        WriteAccess<double> wA( A, hostContext );
    }
    double prefetchFrom = common::Walltime::get() - t;
    cout << "Time for transfer from GPU and workload on Host " << prefetchFrom << endl;
    cout << "Output of dummy = " << dummy << " avoids dead code elimination." << endl;
    double overheadTo = prefetchTo - worktime;
    double overheadFrom = prefetchFrom - worktime;
    // scale to ms
    overheadTo *= 1000.0;
    overheadFrom *= 1000.0;
    transferToGPU *= 1000.0;
    transferFromGPU *= 1000.0;
    cout << endl;
    cout << "Summary of benefits using prefetch, pinned = " << isPinned  << endl;
    cout << "==============================================" << endl;
    cout << endl;
    cout << "Effective time to GPU: " << overheadTo << " instead of " << transferToGPU << " ms" << endl;
    cout << "Effective time from GPU: " << overheadFrom << " instead of " << transferFromGPU << " ms" << endl;
}
