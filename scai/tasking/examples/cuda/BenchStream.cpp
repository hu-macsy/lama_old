/**
 * @file tasking/examples/cuda/BenchStream.cpp
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
 * @brief ToDo: Missing description in ./tasking/examples/cuda/BenchStream.cpp
 * @author Thomas Brandes
 * @date 09.03.2016
 */

#include <scai/common/Walltime.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
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
    CUDACtx device( deviceNr );
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
        CUstream str = pool.reserveStream( StreamType::ComputeStream );
        pool.releaseStream( str );
    }

    CUDAStreamPool::freePool( device );
    double t2 =  Walltime::get() - t0;
    t0 =  Walltime::get();

    for ( int i = 0; i < N_USES; ++i )
    {
        CUDAStreamSyncToken( device, StreamType::ComputeStream );
    }

    double t3 =  Walltime::get() - t0;
    std::cout << "Measure time for " << N_USES << " stream accesses." << std::endl;
    std::cout << "Time create/destroy stream of CUDA: " << ( t1 / N_USES ) * 1000.0 * 1000.0 << " µs" << std::endl;
    std::cout << "Time get/release stream of pool:    " <<  ( t2 / N_USES ) * 1000.0 * 1000.0 << " µs" << std::endl;
    std::cout << "Time construct/destruct SyncToken:  " <<  ( t3 / N_USES ) * 1000.0 * 1000.0 << " µs" << std::endl;
}
