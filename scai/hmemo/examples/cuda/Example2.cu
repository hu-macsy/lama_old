/**
 * @file hmemo/examples/cuda/Example2.cu
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
 * @brief ToDo: Missing description in ./hmemo/examples/cuda/Example2.cu
 * @author Thomas Brandes
 * @date 13.07.2015
 */

#include <scai/hmemo.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <scai/logging.hpp>

#include <scai/common/cuda/CUDAError.hpp>

#include <iostream>

using namespace scai;
using namespace hmemo;
using common::ContextType;

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
}

int main()
{
    ContextPtr cudaContext1;
    ContextPtr cudaContext2;

    try
    {
        std::cout << "try to get " << ContextType::CUDA << " context from factory" << std::endl;
        cudaContext1 = Context::getContextPtr( ContextType::CUDA, 0 );
        std::cout << "cudaContext1 = " << *cudaContext1 << std::endl;
    }
    catch ( scai::common::Exception& ex )
    {
        std::cout << "could not get CUDA device 0" << std::endl;
        return 0;
    }

    try
    {
        std::cout << "try to get " << ContextType::CUDA << " context from factory" << std::endl;
        cudaContext2 = Context::getContextPtr( ContextType::CUDA, 1 );
        std::cout << "cudaContext2 = " << *cudaContext2 << std::endl;
    }
    catch ( scai::common::Exception& ex )
    {
        std::cout << "could not get CUDA device 1" << std::endl;
        return 0;
    }

    std::cout << "try to get " << ContextType::Host << " context from factory" << std::endl;
    ContextPtr hostContext = Context::getContextPtr( ContextType::Host, 1 );
    std::cout << "hostContext = " << *hostContext << std::endl;
    const IndexType N = 100;
    HArray<double> data;
    std::cout << "data = " << data << std::endl;
    {
        SCAI_LOG_INFO( logger, "write only on host" )
        WriteOnlyAccess<double> write( data, hostContext, N );
        double* v = write.get();

        for ( IndexType i = 0; i < N; ++i )
        {
            v[i] = 1.0;
        }
    }
    {
        SCAI_LOG_INFO( logger, "read on gpu 0" )
        ReadAccess<double> read( data, cudaContext1 );
        SCAI_CONTEXT_ACCESS( cudaContext1 )
        double s = sum( read.get(), data.size() );
        std::cout << "sum = " << s << ", should be " << N  << std::endl;
    }
    {
        SCAI_LOG_INFO( logger, "write on gpu1" )
        WriteAccess<double> write( data, cudaContext1 );
        SCAI_CONTEXT_ACCESS( cudaContext1 )
        add( write.get(), data.size() );
    }
    {
        SCAI_LOG_INFO( logger, "write on gpu2" )
        WriteAccess<double> write( data, cudaContext2 );
        SCAI_CONTEXT_ACCESS( cudaContext2 )
        add( write.get(), data.size() );
    }
    {
        SCAI_LOG_INFO( logger, "read on host" )
        ReadAccess<double> read( data, hostContext );
        const double* values = read.get();

        for ( IndexType i = 0; i < N; ++i )
        {
            SCAI_ASSERT_EQUAL( values[i], 3.0, "wrong value after add" )
        }
    }
}

