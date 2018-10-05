/**
 * @file BenchContext.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Benchmarking of memory transfers HOST <-> context
 * @author Thomas Brandes
 * @date 14.09.2015
 */

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <scai/logging.hpp>

using namespace std;
using namespace scai;
using namespace hmemo;
using common::ContextType;
using tasking::TaskSyncToken;

const IndexType NITER = 10;

void runDirect( hmemo::HArray<double>& array, ContextPtr ctx1, ContextPtr ctx2 )
{
    for ( IndexType k = 0; k < NITER; k++ )
    {
        {
            WriteAccess<double>( array, ctx2 );
        }
        {
            WriteAccess<double>( array, ctx1 );
        }
    }
}

void runHost( hmemo::HArray<double>& array, ContextPtr ctx1, ContextPtr ctx2 )
{
    ContextPtr host = Context::getHostPtr();

    for ( IndexType k = 0; k < NITER; k++ )
    {
        {
            WriteAccess<double>( array, host );
        }
        {
            WriteAccess<double>( array, ctx2 );
        }
        {
            WriteAccess<double>( array, host );
        }
        {
            WriteAccess<double>( array, ctx1 );
        }
    }
}

int main()
{
    const IndexType N = 1024 * 1024 * 64;

    ContextPtr host = Context::getHostPtr();
    ContextPtr gpu0 = Context::getContextPtr( ContextType::CUDA, 0 );
    ContextPtr gpu1 = Context::getContextPtr( ContextType::CUDA, 1 );
    ContextPtr gpu2 = Context::getContextPtr( ContextType::CUDA, 2 );
    ContextPtr gpu3 = Context::getContextPtr( ContextType::CUDA, 3 );

    std::cout << "host = " << *host << std::endl;
    std::cout << "gpu0 = " << *gpu0 << std::endl;
    std::cout << "gpu1 = " << *gpu1 << std::endl;
    std::cout << "gpu2 = " << *gpu2 << std::endl;
    std::cout << "gpu3 = " << *gpu3 << std::endl;

    HArray<double> array0( N, 1.0, gpu0 );
    HArray<double> array1( N, 1.0, gpu2 );

    std::cout << "Data allocated, start bench." << std::endl;

    double time0 = common::Walltime::get();

    for ( IndexType k = 0; k < NITER; k++ )
    {
        {
            WriteAccess<double>( array0, host );
        }
        {
            WriteAccess<double>( array0, gpu0 );
        }
    }

    time0 = common::Walltime::get() - time0;
    time0 = time0 / NITER;

    double transfer = double( N ) * sizeof( double ) / ( 1024.0 * 1024.0 * 1024.0 * time0 );

    std::cout << "Time for transfer CPU0 - CPU1: " <<  time0 << std::endl;
    std::cout << "Rate: " << transfer << " GByte/s" << std::endl;

    time0 = common::Walltime::get();

    for ( IndexType k = 0; k < NITER; k++ )
    {
        {
            WriteAccess<double>( array0, host );
        }
        {
            WriteAccess<double>( array0, gpu1 );
        }
        {
            WriteAccess<double>( array0, host );
        }
        {
            WriteAccess<double>( array0, gpu0 );
        }
    }

    time0 = common::Walltime::get() - time0;
    time0 = time0 / NITER;

    transfer = double( N ) * sizeof( double ) / ( 1024.0 * 1024.0 * 1024.0 * time0 );

    std::cout << "Time for transfer GPU0 - GPU1: " <<  time0 << std::endl;
    std::cout << "Rate: " << transfer << " GByte/s" << std::endl;

    time0 = common::Walltime::get();

    for ( IndexType k = 0; k < NITER; k++ )
    {
        {
            WriteAccess<double>( array0, gpu1 );
        }
        {
            WriteAccess<double>( array0, gpu0 );
        }
    }

    time0 = common::Walltime::get() - time0;
    time0 = time0 / NITER;

    transfer = double( N ) * sizeof( double ) / ( 1024.0 * 1024.0 * 1024.0 * time0 );

    std::cout << "Time for transfer GPU0 - GPU1 direct: " <<  time0 << std::endl;
    std::cout << "Rate: " << transfer << " GByte/s" << std::endl;
 
    time0 = common::Walltime::get();

    {
        runDirect( array1, gpu2, gpu3 );
    }

    time0 = common::Walltime::get() - time0;
    time0 = time0 / NITER;

    transfer = double( N ) * sizeof( double ) / ( 1024.0 * 1024.0 * 1024.0 * time0 );

    std::cout << "Time for transfer GPU2 - GPU3 direct: " <<  time0 << std::endl;
    std::cout << "Rate: " << transfer << " GByte/s" << std::endl;

    time0 = common::Walltime::get();

    {
        TaskSyncToken* t = new TaskSyncToken( std::bind( runHost, array0, gpu0, gpu1 ) );
        runHost( array1, gpu2, gpu3 );
        delete t;
    }

    time0 = common::Walltime::get() - time0;
    time0 = time0 / NITER;

    transfer = double( N ) * sizeof( double ) / ( 1024.0 * 1024.0 * 1024.0 * time0 );

    std::cout << "Time for transfer GPU0 - GPU1, GPU2 - GPU3 host: " <<  time0 << std::endl;
    std::cout << "Rate: " << transfer << " GByte/s" << std::endl;

    time0 = common::Walltime::get();

    {
        TaskSyncToken* t = new TaskSyncToken( std::bind( runDirect, array0, gpu0, gpu1 ) );
        runDirect( array1, gpu2, gpu3 );
        delete t;
    }

    time0 = common::Walltime::get() - time0;
    time0 = time0 / NITER;

    transfer = double( N ) * sizeof( double ) / ( 1024.0 * 1024.0 * 1024.0 * time0 );

    std::cout << "Time for transfer GPU0 - GPU1, GPU2 - GPU3 direct: " <<  time0 << std::endl;
    std::cout << "Rate: " << transfer << " GByte/s" << std::endl;
}

/* Example output:

    host = HostContext( #Threads = 32 )
    gpu0 = CUDAContext(0: Tesla P100-PCIE-16GB)
    gpu1 = CUDAContext(1: Tesla P100-PCIE-16GB)
    gpu2 = CUDAContext(2: Tesla P100-PCIE-16GB)
    gpu3 = CUDAContext(3: Tesla P100-PCIE-16GB)
    Data allocated, start bench.
    Time for transfer CPU0 - CPU1: 0.0844445
    Rate: 5.92105 GByte/s
    Time for transfer GPU0 - GPU1: 0.168475
    Rate: 2.96779 GByte/s
    Time for transfer GPU0 - GPU1 direct: 0.0903717
    Rate: 5.5327 GByte/s
    Time for transfer GPU2 - GPU3 direct: 0.0910597
    Rate: 5.49091 GByte/s
    Time for transfer GPU0 - GPU1, GPU2 - GPU3 host: 0.366881
    Rate: 1.36284 GByte/s
    Time for transfer GPU0 - GPU1, GPU2 - GPU3 direct: 0.102825
    Rate: 4.86265 GByte/s

*/
