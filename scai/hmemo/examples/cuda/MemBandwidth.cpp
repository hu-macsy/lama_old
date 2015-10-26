/**
 * @file MemBandwidth.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
using namespace scai::hmemo;

SCAI_LOG_DEF_LOGGER( logger, "MemBandwidth" )

template <typename ValueType> 
void bench( LAMAArray<ValueType>& array )
{
    ContextPtr cudaContext = Context::getContextPtr( common::context::CUDA );
    ContextPtr hostContext = Context::getContextPtr( common::context::Host );

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
    ContextPtr cudaContext = Context::getContextPtr( common::context::CUDA );
    ContextPtr hostContext = Context::getContextPtr( common::context::Host );

    LAMAArray<float> A1( hostContext );  // same as LAMAArray<float> A1;
    LAMAArray<float> A2( cudaContext );

    bench( A1 );
    bench( A2 );
}

