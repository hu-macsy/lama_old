/**
 * @file MemBandwidth.cpp
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
 * @brief Measuring the CUDA Memory Bandwidth
 * @author Thomas Brandes
 * @date 16.07.2015
 */

#include <scai/hmemo.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <scai/logging.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>

using namespace scai;
using namespace hmemo;
using common::ContextType;

SCAI_LOG_DEF_LOGGER( logger, "MemBandwidth" )

template <typename ValueType>
void bench( HArray<ValueType>& array )
{
    ContextPtr cudaContext = Context::getContextPtr( common::ContextType::CUDA );
    ContextPtr hostContext = Context::getContextPtr( common::ContextType::Host );
    const IndexType N = 8 * 1024 * 1024;
    const int NITER = 128;
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
    ContextPtr cudaContext = Context::getContextPtr( ContextType::CUDA );
    ContextPtr hostContext = Context::getContextPtr( ContextType::Host );
    HArray<float> A1( hostContext );  // same as HArray<float> A1;
    HArray<float> A2( cudaContext );
    bench( A1 );
    bench( A2 );
}

