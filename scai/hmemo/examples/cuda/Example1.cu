/**
 * @file hmemo/examples/cuda/Example1.cu
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
 * @brief ToDo: Missing description in ./hmemo/examples/cuda/Example1.cu
 * @author Thomas Brandes
 * @date 10.07.2015
 */

#include <scai/hmemo.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <scai/logging.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Walltime.hpp>

#include <iostream>
#include <unistd.h>

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
void add( ValueType* array, const IndexType n )
{
    const int blockSize = 256;
    const int nblocks   = ( n + blockSize - 1 ) / blockSize;
    dim3 block( blockSize, 1, 1 );
    dim3 grid( nblocks, 1, 1 );
    add_kernel <<< grid, block>>>( array, n );
}

int main()
{
    std::cout << "try to get " << ContextType::CUDA << " context from factory" << std::endl;
    ContextPtr cudaContext = Context::getContextPtr( ContextType::CUDA );
    SCAI_ASSERT( cudaContext, "NULL context" )
    std::cout << "cudaContext = " << *cudaContext << std::endl;
    MemoryPtr cudaMemory = cudaContext->getMemoryPtr();
    SCAI_ASSERT( cudaMemory, "NULL memory" )
    std::cout << "cudaMemory = " << *cudaMemory << std::endl;
    std::cout << "try to get " << ContextType::Host << " context from factory" << std::endl;
    ContextPtr hostContext = Context::getContextPtr( ContextType::Host );
    SCAI_ASSERT( hostContext, "NULL context" )
    std::cout << "hostContext = " << *hostContext << std::endl;
    MemoryPtr hostMemory = hostContext->getMemoryPtr();
    SCAI_ASSERT( hostMemory, "NULL memory" )
    std::cout << "hostMemory = " << *hostMemory << std::endl;
    MemoryPtr cudaHostMemory = cudaContext->getHostMemoryPtr();
    SCAI_ASSERT( cudaHostMemory, "NULL memory" )
    std::cout << "cudaHostMemory = " << *cudaHostMemory << std::endl;
    const IndexType N = 100;
    HArray<double> data( cudaContext );
    std::cout << "data = " << data << std::endl;
    {
        SCAI_LOG_INFO( logger, "write only on cuda host" )
        // HostWriteOnlyAccess<double> write( data,  N );
        WriteOnlyAccess<double> writeData( data, N );
        double* dataHost = writeData;

        for ( IndexType i = 0; i < N; ++i )
        {
            writeData[i] = 1.0;
        }
    }
    std::cout << "After host write: data = " << data << std::endl;
    {
        SCAI_LOG_INFO( logger, "read on cuda" )
        ReadAccess<double> read( data, cudaContext );
        SCAI_CONTEXT_ACCESS( cudaContext )
        double s = sum( read.get(), data.size() );
        std::cout << "sum = " << s << ", should be " << N  << std::endl;
    }
    std::cout << "After cuda read: data = " << data << std::endl;
    {
        SCAI_LOG_INFO( logger, "write on cuda" )
        WriteAccess<double> write( data, cudaContext );
        SCAI_CONTEXT_ACCESS( cudaContext )
        add( static_cast<double*>( write ), data.size() );
    }
    std::cout << "After cuda write: data = " << data << std::endl;
    {
        SCAI_LOG_INFO( logger, "read on host" )
        ReadAccess<double> read( data, hostContext );
        common::Walltime::sleep( 1000 );

        for ( IndexType i = 0; i < N; ++i )
        {
            SCAI_ASSERT_EQUAL( read[i], 2 * 1.0, "wrong value after add, i = " << i )
        }
    }
    std::cout << "After host read: data = " << data << std::endl;
}

