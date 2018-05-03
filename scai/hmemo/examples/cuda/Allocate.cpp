/**
 * @file Allocate.cpp
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
 * @brief Benchmark allocate on Host vs CUDAHost memory
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#include <scai/hmemo.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <memory>

using namespace scai;
using namespace scai::hmemo;
using namespace scai::common;

using std::shared_ptr;

void doit( int NITER, IndexType NSIZE, ContextPtr context )
{
    std::vector<shared_ptr<HArray<double> > > stack;
    ContextPtr host = Context::getHostPtr();

    for ( int iter = 0; iter < NITER; ++iter )
    {
        shared_ptr<HArray<double> > X;
        X.reset( new HArray<double>( context ) );
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
    ContextPtr cudaContext = Context::getContextPtr( common::ContextType::CUDA );
    ContextPtr hostContext = Context::getContextPtr( common::ContextType::Host );
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
