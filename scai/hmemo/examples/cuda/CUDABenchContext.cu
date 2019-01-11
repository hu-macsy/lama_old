/**
 * @file hmemo/examples/cuda/CUDABenchContext.cu
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
 * @brief ToDo: Missing description in ./hmemo/examples/cuda/CUDABenchContext.cu
 * @author Thomas Brandes
 * @date 21.07.2015
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

SCAI_LOG_DEF_LOGGER( logger, "CudaExample" )

template<typename ValueType>
ValueType sum( const ValueType array[], const IndexType n )
{
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
    ValueType zero = static_cast<ValueType>( 0 );
    ValueType result = thrust::reduce( data, data + n, zero, thrust::plus<ValueType>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "sum of " << n << " values = " << result )
    return result;
}

template<typename ValueType>
__global__
void add_kernel( ValueType* array, IndexType n )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;
    ValueType one = 1;

    if ( i < n )
    {
        array[i] += one;
    }
}

template<typename ValueType>
void add( ValueType array[], const IndexType n )
{
    const int blockSize = 256;
    const int nblocks   = ( n + blockSize - 1 ) / blockSize;
    dim3 block( blockSize, 1, 1 );
    dim3 grid( nblocks, 1, 1 );
    add_kernel <<< grid, block>>>( array, n );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cuda failure" );
}

template<typename ValueType>
void addHost( ValueType array[], const IndexType n )
{
    for ( IndexType i = 0; i < n; ++i )
    {
        array[i] += 1;
    }
}

template<typename ValueType>
void doBench( HArray<ValueType>& array, const IndexType N )
{
    ContextPtr hostContext = Context::getContextPtr( common::ContextType::Host );
    ContextPtr cudaContext = Context::getContextPtr( common::ContextType::CUDA );
    int nhost = 1;
    int ncuda = 1;
    int niter = 50;
    double time = common::Walltime::get();
    // init on host
    {
        WriteOnlyAccess<double> write( array, hostContext, N );
        double* v = write.get();

        for ( IndexType i = 0; i < N; ++i )
        {
            v[i] = 0.0;
        }
    }

    for ( int iter = 0; iter < niter; ++iter )
    {
        // do some work on cuda
        for ( int k = 0; k < ncuda; ++k )
        {
            WriteAccess<double> write( array, cudaContext );
            SCAI_CONTEXT_ACCESS( cudaContext )
            add( write.get(), N );
        }

        // do some work on host

        for ( int k = 0; k < nhost; ++k )
        {
            WriteAccess<double> write( array, hostContext );
            addHost( write.get(), N );
        }
    }

    // compute result
    double res = 0.0;
    {
        ReadAccess<double> read( array, cudaContext );
        SCAI_CONTEXT_ACCESS( cudaContext )
        res = sum( read.get(), N );
    }
    double resExpected = N;
    resExpected *= double ( niter * ( ncuda + nhost ) );
    SCAI_ASSERT_EQUAL( res, resExpected, "wrong result, N = " << N
                       << ", niter = " << niter << ", ncuda = " << ncuda << ", nhost = " << nhost )
    time = common::Walltime::get() - time;
    std::cout << "Time = " << time << " seconds" << std::endl;
}

int main()
{
    const IndexType N = 8 * 1024 * 1024;  // 8 MB data
    ContextPtr hostContextPtr = Context::getContextPtr( common::ContextType::Host );
    ContextPtr cudaContextPtr = Context::getContextPtr( common::ContextType::CUDA );
    // First touch on host memory, never uses CUDA host memory
    std::cout << "Benchmark for array, first touch on host memory" << std::endl;
    HArray<double> data1( hostContextPtr->getMemoryPtr() );
    doBench( data1, N );
    std::cout << "Benchmark for array, first touch on cuda memory" << std::endl;
    HArray<double> data2( cudaContextPtr->getMemoryPtr() );
    doBench( data2, N );
    std::cout << "Benchmark for array, first touch on cuda host memory" << std::endl;
    HArray<double> data3( cudaContextPtr->getHostMemoryPtr() );
    doBench( data3, N );
}

