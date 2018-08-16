/**
 * @file common/examples/cuda/BenchCUDA.cpp
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
 * @brief Benchmark of some CUDA stuff
 * @author Thomas Brandes
 * @date 11.02.2016
 */

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <scai/common/cuda/CUDACtx.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/CUDAError.hpp>

#include <iostream>

using namespace scai::common;

int main( int argc, const char** argv )
{
    // at least --SCAI_DEVICE=id may be specified
    Settings::parseArgs( argc, argv );
    int nr = 0;   // take this as default
    Settings::getEnvironment( nr, "SCAI_DEVICE" );
    {
        // do not measure overhead for 1st allocation
        CUDACtx dev( nr );
    }
    double t0 = Walltime::get();
    const int N_CONTEXT = 10;

    for ( int i = 0; i < N_CONTEXT; ++i )
    {
        CUDACtx dev( nr );
    }

    double t1 = Walltime::get() - t0;
    t1 = t1 / N_CONTEXT;
    std::cout << "Time for CUDACtx(..) = " << ( t1 * 1000.0 * 1000.0 ) << " µs" << std::endl;
    const int N_ACCESS = 10000;
    t0 = Walltime::get();
    {
        CUDACtx dev( nr );

        for ( int i = 0; i < N_ACCESS; ++i )
        {
            SCAI_CUDA_DRV_CALL( cuCtxPushCurrent( dev.getCUcontext() ), "could not push context" )
            CUcontext tmp; // temporary for last context, not necessary to save it
            SCAI_CUDA_DRV_CALL( cuCtxPopCurrent( &tmp ), "could not pop context" )
        }
    }
    t1 = Walltime::get() - t0;
    t1 = t1 / N_ACCESS;
    std::cout << "Time for push/pop context = " << ( t1 * 1000.0 * 1000.0 ) << " µs" << std::endl;
    t0 = Walltime::get();
    {
        CUDACtx dev( nr );

        for ( int i = 0; i < N_ACCESS; ++i )
        {
            CUDAAccess tmp( dev );
        }
    }
    t1 = Walltime::get() - t0;
    t1 = t1 / N_ACCESS;
    std::cout << "Time for CUDAAccess(..) = " << ( t1 * 1000.0 * 1000.0 ) << " µs" << std::endl;
}
