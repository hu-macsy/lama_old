/**
 * @file CUDACSRUtils.cu
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of CSR utilities with CUDA
 * @author Bea Hornef, Thomas Brandes, Jiri Kraus
 * @date 04.07.2012
 */

// hpp
#include <scai/sparsekernel/cuda/CUDACSRUtils.hpp>

// local library
#include <scai/sparsekernel/cuda/CUDACSRUtils.hpp>
#include <scai/sparsekernel/cuda/CUDACOOUtils.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

// internal scai library
#include <scai/utilskernel/cuda/CUDAUtils.hpp>

#include <scai/hmemo/Memory.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAUtils.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/Constants.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>

#include <scai/common/macros/unused.hpp>

// CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

// thrust
#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/reduce.h>

// Parameters for Matrix Multiplication
#define NUM_HASH_RETRIES 16
#define NUM_ELEMENTS_PER_CHUNK 512
#define NUM_ELEMENTS_IN_SHARED 512
#define NUM_BLOCKS 9216
#define NUM_THREADS 32
#define NUM_WARPS NUM_THREADS/32
#define HASH_A 684
#define HASH_B 46165
#define HASH_P 88651
#define HASH_C0 1
#define HASH_C1 1
#define NUM_CHUNKS_PER_WARP 128

using namespace scai::common;
using namespace scai::hmemo;

namespace scai
{

using utilskernel::CUDAUtils;

using tasking::SyncToken;
using tasking::CUDAStreamSyncToken;

using common::CUDASettings;


namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( CUDACSRUtils::logger, "CUDA.CSRUtils" )

// not yet: __device__ const IndexType cudaNIndex = std::numeric_limits<IndexType>::max();

#define cudaNIndex static_cast<IndexType>( -1 )

IndexType CUDACSRUtils::sizes2offsets( IndexType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "sizes2offsets " << " #n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<IndexType> array_ptr( array );
    thrust::exclusive_scan( array_ptr, array_ptr + n + 1, array_ptr );
    thrust::host_vector<IndexType> numValues( array_ptr + n, array_ptr + n + 1 );
    return numValues[0];
}

/* --------------------------------------------------------------------------- */
/*     CUDA Kernels                                                            */
/* --------------------------------------------------------------------------- */

__global__
static void offsets2sizes_kernel( IndexType sizes[], const IndexType offsets[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        sizes[i] = offsets[i + 1] - offsets[i];
    }
}

/* --------------------------------------------------------------------------- */
/*     offsets2sizes                                                           */
/* --------------------------------------------------------------------------- */

void CUDACSRUtils::offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType n )
{
    SCAI_REGION( "CUDA.CSRUtils.offsets2sizes" )
    SCAI_LOG_INFO( logger, "offsets2sizes " << " #n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );
    offsets2sizes_kernel <<< dimGrid, dimBlock>>>( sizes, offsets, n );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "offsets2sizes" )
}

/* --------------------------------------------------------------------------- */
/*     hasDiagonalProperty                                                     */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
struct identic_functor
{
    __host__ __device__
    double operator()( thrust::tuple<ValueType, ValueType> x )
    {
        return thrust::get < 0 > ( x ) == thrust::get < 1 > ( x );
    }
};

//trivial kernel to check diagonal property
__global__ void hasDiagonalProperty_kernel(
    const IndexType numDiagonals,
    const IndexType ia[],
    const IndexType ja[],
    bool* hasProperty )
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= numDiagonals )
    {
        return;
    }

    if ( ! ( *hasProperty ) )
    {
        return;
    }

    if ( ia[i] == ia[i + 1] )
    {
        *hasProperty = false;
    }
    else if ( ja[ia[i]] != i )
    {
        *hasProperty = false;
    }
}

bool CUDACSRUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrIA[], const IndexType csrJA[] )
{
    SCAI_REGION( "CUDA.CSRUtils.hasDiagonalProperty" )

    if ( numDiagonals == 0 )
    {
        return true;
    }

    SCAI_CHECK_CUDA_ACCESS
    //make grid
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimGrid( ( numDiagonals - 1 ) / blockSize + 1, 1, 1 );// = makeGrid( numDiagonals, blockSize );
    dim3 dimBlock( blockSize, 1, 1 );
    bool* d_hasProperty;
    bool hasProperty;
    SCAI_CUDA_RT_CALL( cudaMalloc( ( void** ) &d_hasProperty, sizeof( bool ) ),
                       "allocate 4 bytes on the device for the result of hasDiagonalProperty_kernel" )
    SCAI_CUDA_RT_CALL( cudaMemset( d_hasProperty, 1, sizeof( bool ) ), "memset bool hasProperty = true" )
    hasDiagonalProperty_kernel <<< dimGrid, dimBlock>>>( numDiagonals, csrIA, csrJA, d_hasProperty );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "hasDiagonalProperty failed: are ia and ja correct?" )
    SCAI_CUDA_RT_CALL( cudaMemcpy( &hasProperty, d_hasProperty, sizeof( bool ), cudaMemcpyDeviceToHost ),
                       "copy the result of hasDiagonalProperty_kernel to host" )
    return hasProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDACSRUtils::convertCSR2CSC(
    IndexType cscIA[],
    IndexType cscJA[],
    ValueType cscValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    IndexType numRows,
    IndexType numColumns,
    IndexType numValues )
{
    SCAI_REGION( "CUDA.CSRUtils.CSR2CSC" )
    SCAI_LOG_INFO( logger, "convertCSR2CSC of " << numRows << " x " << numColumns << ", nnz = " << numValues )
    // Sort the csrJA ( same as cooJA ), apply it to cooIA and cooValues
    IndexType* cooIA;
    SCAI_CUDA_RT_CALL( cudaMalloc( &cooIA, sizeof( IndexType ) * numValues ),
                       "allocate temp for cooIA" )
    // Step 1 : build COO storage,  cooIA (to do), cooJA ( = csrJA ), cooValues ( = csrValues )
    //          -> translate the csrIA offset array to a cooIA array
    const IndexType numDiagonals = 0;// not supported yet
    CUDACOOUtils::offsets2ia( cscJA, numValues, csrIA, numRows, numDiagonals );
    // switch cooIA and cooJA, copy values and resort
    CUDAUtils::set( cooIA, csrJA, numValues, utilskernel::reduction::COPY );
    CUDAUtils::set( cscValues, csrValues, numValues, utilskernel::reduction::COPY );
    thrust::device_ptr<IndexType> ja_d( cooIA );
    thrust::device_ptr<ValueType> values_d( cscValues );
    thrust::device_ptr<IndexType> ia_d( cscJA );
    // sort by column indexes in ascending order
    // zip_iterator used to resort cscValues and cscJA in one step
    thrust::stable_sort_by_key( ja_d, ja_d + numValues,
                                thrust::make_zip_iterator( thrust::make_tuple( values_d, ia_d ) ) );
    // cscJA is now sorted, can become an offset array
    CUDACOOUtils::ia2offsets( cscIA, numColumns, 0, cooIA, numValues );
    SCAI_CUDA_RT_CALL( cudaFree( cooIA ), "free tmp cooIA" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void scale_kernel(
    ValueType* result,
    const ValueType* y_d,
    const ValueType beta,
    IndexType numRows )
{
    // result = beta * y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        result[i] = beta * y_d[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gemv_kernel_beta_zero(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType alpha,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows )
{
    // result = alpha * A * x_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        ValueType value = 0.0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetchVectorX<ValueType, useTexture>( x_d, csrJA[jj] );
        }

        result[i] = alpha * value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gemv_kernel_alpha_one(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType beta,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows )
{
    // result = A * x_d + beta * y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType summand = beta * y_d[i];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        ValueType value = 0.0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetchVectorX<ValueType, useTexture>( x_d, csrJA[jj] );
        }

        result[i] = value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gemv_kernel_beta_one(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType alpha,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows )
{
    // result = alpha * A * x_d + y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType summand = y_d[i];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        ValueType value = 0.0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetchVectorX<ValueType, useTexture>( x_d, csrJA[jj] );
        }

        result[i] = alpha * value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gemv_kernel_alpha_one_beta_zero(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows )
{
    // result = A * x_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        ValueType value = 0.0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetchVectorX<ValueType, useTexture>( x_d, csrJA[jj] );
        }

        result[i] = value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void assign_kernel(
    ValueType* result,
    const ValueType* y_d,
    IndexType numRows )
{
    // result = y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        result[i] = y_d[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gemv_kernel_alpha_one_beta_one(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows )
{
    // result = A * x_d + y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType summand = y_d[i];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        ValueType value = 0.0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetchVectorX<ValueType, useTexture>( x_d, csrJA[jj] );
        }

        result[i] = value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gemv_kernel(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType alpha,
    const ValueType beta,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows )
{
    // result = alpha * A * x_d + beta * y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType summand = beta * y_d[i];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        ValueType value = 0.0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetchVectorX<ValueType, useTexture>( x_d, csrJA[jj] );
        }

        result[i] = alpha * value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gevm_kernel_alpha_one_beta_one(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows,
    IndexType numColumns )
{
    // result = x_d * A + y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        ValueType summand = y_d[i];
        ValueType value = 0.0;

        for ( IndexType j = 0; j < numRows; ++j )
        {
            const IndexType rowStart = csrIA[j];
            const IndexType rowEnd = csrIA[j + 1];

            for ( IndexType k = rowStart; k < rowEnd; ++k )
            {
                if ( csrJA[k] == i )
                {
                    value += csrValues[k] * fetchVectorX<ValueType, useTexture>( x_d, j );
                }
            }
        }

        result[i] = value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gevm_kernel_alpha_one_beta_zero(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows,
    IndexType numColumns )
{
    // result = x_d * A
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        ValueType value = 0.0;

        for ( IndexType j = 0; j < numRows; ++j )
        {
            const IndexType rowStart = csrIA[j];
            const IndexType rowEnd = csrIA[j + 1];

            for ( IndexType k = rowStart; k < rowEnd; ++k )
            {
                if ( csrJA[k] == i )
                {
                    value += csrValues[k] * fetchVectorX<ValueType, useTexture>( x_d, j );
                }
            }
        }

        result[i] = value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gevm_kernel_alpha_one(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType beta,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows,
    IndexType numColumns )
{
    // result = x_d * A + beta * y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        ValueType summand = beta * y_d[i];
        ValueType value = 0.0;

        for ( IndexType j = 0; j < numRows; ++j )
        {
            const IndexType rowStart = csrIA[j];
            const IndexType rowEnd = csrIA[j + 1];

            for ( IndexType k = rowStart; k < rowEnd; ++k )
            {
                if ( csrJA[k] == i )
                {
                    value += csrValues[k] * fetchVectorX<ValueType, useTexture>( x_d, j );
                }
            }
        }

        result[i] = value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gevm_kernel_beta_one(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType alpha,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows,
    IndexType numColumns )
{
    // result = alpha * x_d * A + y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        ValueType value = 0.0;

        for ( IndexType j = 0; j < numRows; ++j )
        {
            const IndexType rowStart = csrIA[j];
            const IndexType rowEnd = csrIA[j + 1];

            for ( IndexType k = rowStart; k < rowEnd; ++k )
            {
                if ( csrJA[k] == i )
                {
                    value += csrValues[k] * fetchVectorX<ValueType, useTexture>( x_d, j );
                }
            }
        }

        result[i] = alpha * value + y_d[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gevm_kernel_beta_zero(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType alpha,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows,
    IndexType numColumns )
{
    // result = alpha * x_d * A
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        ValueType value = 0.0;

        for ( IndexType j = 0; j < numRows; ++j )
        {
            const IndexType rowStart = csrIA[j];
            const IndexType rowEnd = csrIA[j + 1];

            for ( IndexType k = rowStart; k < rowEnd; ++k )
            {
                if ( csrJA[k] == i )
                {
                    value += csrValues[k] * fetchVectorX<ValueType, useTexture>( x_d, j );
                }
            }
        }

        result[i] = alpha * value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void normal_gevm_kernel(
    ValueType* result,
    const ValueType* x_d,
    const ValueType* y_d,
    const ValueType alpha,
    const ValueType beta,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    IndexType numRows,
    IndexType numColumns )
{
    // result = alpha * x_d * A + beta * y_d
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        ValueType summand = beta * y_d[i];
        ValueType value = 0.0;

        for ( IndexType j = 0; j < numRows; ++j )
        {
            const IndexType rowStart = csrIA[j];
            const IndexType rowEnd = csrIA[j + 1];

            for ( IndexType k = rowStart; k < rowEnd; ++k )
            {
                if ( csrJA[k] == i )
                {
                    value += csrValues[k] * fetchVectorX<ValueType, useTexture>( x_d, j );
                }
            }
        }

        result[i] = alpha * value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void sparse_gemv_kernel_alpha_one(
    ValueType* result,
    const ValueType* x_d,
    const ValueType alpha,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    const IndexType* rowIndexes,
    IndexType numRows )
{
    // result = A * x_d
    const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( ii < numRows )
    {
        IndexType i = rowIndexes[ii];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        ValueType value = 0.0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetchVectorX<ValueType, useTexture>( x_d, csrJA[jj] );
        }

        result[i] += value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void sparse_gemv_kernel(
    ValueType* result,
    const ValueType* x_d,
    const ValueType alpha,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    const IndexType* rowIndexes,
    IndexType numRows )
{
    // result = alpha * A * x_d
    const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( ii < numRows )
    {
        IndexType i = rowIndexes[ii];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        ValueType value = 0.0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetchVectorX<ValueType, useTexture>( x_d, csrJA[jj] );
        }

        result[i] += alpha * value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void sparse_gevm_kernel(
    ValueType* result,
    const ValueType* x_d,
    const ValueType alpha,
    const ValueType* csrValues,
    const IndexType* csrIA,
    const IndexType* csrJA,
    const IndexType* rowIndexes,
    IndexType numColumns,
    IndexType numNonZeroRows )
{
    // result += alpha * x_d * A
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        ValueType value = 0.0;

        for ( IndexType jj = 0; jj < numNonZeroRows; ++jj )
        {
            IndexType j = rowIndexes[jj];
            const IndexType rowStart = csrIA[j];
            const IndexType rowEnd = csrIA[j + 1];

            for ( IndexType k = rowStart; k < rowEnd; ++k )
            {
                if ( csrJA[k] == i )
                {
                    value += csrValues[k] * fetchVectorX<ValueType, useTexture>( x_d, j );
                }
            }
        }

        result[i] += alpha * value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  scaleRows                                                         */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, typename OtherValueType>
__global__
void scaleRowsKernel(
    ValueType* values,
    const IndexType* ia,
    const IndexType numRows,
    const OtherValueType* diagonal )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType tmp = static_cast<OtherValueType>( diagonal[i] );

        for ( IndexType j = ia[i]; j < ia[i + 1]; ++j )
        {
            values[j] *= tmp;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void CUDACSRUtils::scaleRows(
    ValueType1 csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const ValueType2 values[] )
{
    SCAI_REGION( "CUDA.CSRUtils.scaleRows" )
    SCAI_LOG_INFO( logger, "scaleRows<" << TypeTraits<ValueType1>::id() << ","
                   << TypeTraits<ValueType2>::id() << ">"
                   << ", numrows= " << numRows )
    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    scaleRowsKernel <<< dimGrid, dimBlock>>>( csrValues, csrIA, numRows, values );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "CSRUtils:scaleRowsKernel FAILED" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDACSRUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType SCAI_UNUSED( numColumns ),
    const IndexType SCAI_UNUSED( nnz ),
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_REGION( "CUDA.CSRUtils.normalGEMV" )
    SCAI_LOG_INFO( logger, "normalGEMV<" << TypeTraits<ValueType>::id() << ">" <<
                   " result[ " << numRows << "] = " << alpha << " * A(csr) * x + " << beta << " * y " )
    SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0; // default stream if no syncToken is given
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    bool useTexture = CUDASettings::useTexture();
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // asynchronous execution takes other stream and will not synchronize later
        stream = syncToken->getCUDAStream();
    }

    SCAI_LOG_INFO( logger, "Start normal_gemv_kernel<" << TypeTraits<ValueType>::id()
                   << ", useTexture = " << useTexture << ">" );

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == constants::ONE && beta == constants::ONE )
        {
            // result = A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == constants::ONE && beta == constants::ZERO )
        {
            // result = A * x_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == constants::ZERO && beta == constants::ONE )
        {
            // result = y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( assign_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            assign_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>( result, y, numRows );
        }
        else if ( alpha == constants::ONE )
        {
            // result = A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, beta, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == constants::ZERO )
        {
            // result = A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( scale_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            scale_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>( result, y, beta, numRows );
        }
        else if ( beta == constants::ONE )
        {
            // result = alpha * A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows );
        }
        else if ( beta == constants::ZERO )
        {
            // result = alpha * A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows );
        }
        else
        {
            // result = alpha * A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, beta, csrValues, csrIA, csrJA, numRows );
        }
    }
    else
    {
        if ( alpha == constants::ONE && beta == constants::ONE )
        {
            // result = A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == constants::ONE && beta == constants::ZERO )
        {
            // result = A * x_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == constants::ZERO && beta == constants::ONE )
        {
            // result = y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( assign_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            assign_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>( result, y, numRows );
        }
        else if ( alpha == constants::ONE )
        {
            // result = A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, beta, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == constants::ZERO )
        {
            // result = beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( scale_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            scale_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>( result, y, beta, numRows );
        }
        else if ( beta == constants::ONE )
        {
            // result = alpha * A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows );
        }
        else if ( beta == constants::ZERO )
        {
            // result = alpha * A * x_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows );
        }
        else
        {
            // result = alpha * A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, beta, csrValues, csrIA, csrJA, numRows );
        }
    }

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "normalGEMV, stream = " << stream )
        SCAI_LOG_DEBUG( logger, "normalGEMV<" << TypeTraits<ValueType>::id() << "> synchronized" )
    }

    if ( useTexture )
    {
        if ( !syncToken )
        {
            vectorUnbindTexture( x );
        }
        else
        {
            // get routine with the right signature
            void ( *unbind ) ( const ValueType* ) = &vectorUnbindTexture;
            // delay unbind until synchroniziaton
            syncToken->pushRoutine( common::bind( unbind, x ) );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDACSRUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger, "normalGEVM<" << TypeTraits<ValueType>::id() << ">" <<
                   " result[ " << numColumns << "] = " << alpha << " * A(csr) * x + " << beta << " * y " )
    SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0; // default stream if no syncToken is given
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numColumns, dimBlock.x );
    bool useTexture = CUDASettings::useTexture();
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    SCAI_LOG_INFO( logger, "Start normal_gevm_kernel<" << TypeTraits<ValueType>::id()
                   << ", useTexture = " << useTexture << ">" );

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == constants::ONE && beta == constants::ONE )
        {
            // result = A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_alpha_one_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else if ( alpha == constants::ONE && beta == constants::ZERO )
        {
            // result = A * x_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_alpha_one_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else if ( alpha == constants::ZERO && beta == constants::ONE )
        {
            // result = y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( assign_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            assign_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>( result, y, numColumns );
        }
        else if ( alpha == constants::ONE )
        {
            // result = A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, beta, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else if ( alpha == constants::ZERO )
        {
            // result = beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( scale_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            scale_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>( result, y, beta, numColumns );
        }
        else if ( beta == constants::ONE )
        {
            // result = alpha * A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else if ( beta == constants::ZERO )
        {
            // result = alpha * A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else
        {
            // result = alpha * A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, beta, csrValues, csrIA, csrJA, numRows, numColumns );
        }
    }
    else
    {
        if ( alpha == constants::ONE && beta == constants::ONE )
        {
            // result = A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_alpha_one_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else if ( alpha == constants::ONE && beta == constants::ZERO )
        {
            // result = A * x_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_alpha_one_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else if ( alpha == constants::ZERO && beta == constants::ONE )
        {
            // result = y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( assign_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            assign_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>( result, y, numColumns );
        }
        else if ( alpha == constants::ONE )
        {
            // result = A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, beta, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else if ( alpha == constants::ZERO )
        {
            // result = beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( scale_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            scale_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>( result, y, beta, numColumns );
        }
        else if ( beta == constants::ONE )
        {
            // result = alpha * A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else if ( beta == constants::ZERO )
        {
            // result = alpha * A * x_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows, numColumns );
        }
        else
        {
            // result = alpha * A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gevm_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, beta, csrValues, csrIA, csrJA, numRows, numColumns );
        }
    }

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "normalGEVM, stream = " << stream )
        SCAI_LOG_DEBUG( logger, "normalGEVM<" << TypeTraits<ValueType>::id() << "> synchronized" )
    }

    if ( useTexture )
    {
        if ( !syncToken )
        {
            vectorUnbindTexture( x );
        }
        else
        {
            // get routine with the right signature
            void ( *unbind ) ( const ValueType* ) = &vectorUnbindTexture;
            // delay unbind until synchroniziaton
            syncToken->pushRoutine( common::bind( unbind, x ) );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDACSRUtils::sparseGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_REGION( "CUDA.CSRUtils.sparseGEMV" )
    SCAI_LOG_INFO( logger,
                   "sparseGEMV<" << TypeTraits<ValueType>::id() << ">" << ", #non-zero rows = " << numNonZeroRows )
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    const int blockSize = CUDASettings::getBlockSize( numNonZeroRows );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numNonZeroRows, dimBlock.x );
    bool useTexture = CUDASettings::useTexture();

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == constants::ONE )
        {
            sparse_gemv_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numNonZeroRows );
        }
        else
        {
            sparse_gemv_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numNonZeroRows );
        }
    }
    else
    {
        if ( alpha == constants::ONE )
        {
            sparse_gemv_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numNonZeroRows );
        }
        else
        {
            sparse_gemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numNonZeroRows );
        }
    }

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "sparseGEMV, stream = " << stream )
        SCAI_LOG_INFO( logger, "sparseGEMV<" << TypeTraits<ValueType>::id() << "> synchronized" )
    }

    if ( useTexture )
    {
        if ( !syncToken )
        {
            vectorUnbindTexture( x );
        }
        else
        {
            // get routine with the right signature
            void ( *unbind ) ( const ValueType* ) = &vectorUnbindTexture;
            // delay unbind until synchroniziaton
            syncToken->pushRoutine( common::bind( unbind, x ) );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDACSRUtils::sparseGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numColumns,
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "sparseGEVM<" << TypeTraits<ValueType>::id() << ">" << ", #non-zero rows = " << numNonZeroRows )
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;
    // check if asynchronous execution is wanted
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    const int blockSize = CUDASettings::getBlockSize( numNonZeroRows );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numNonZeroRows, dimBlock.x );
    bool useTexture = CUDASettings::useTexture();

    if ( useTexture )
    {
        vectorBindTexture( x );
        SCAI_LOG_DEBUG( logger, "sparse_gevm_kernel<useTexture=true>" )
        sparse_gevm_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
        ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numColumns, numNonZeroRows );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "sparse_gevm_kernel<useTexture=false>" )
        sparse_gevm_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
        ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numColumns, numNonZeroRows );
    }

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "sparseGEVM, stream = " << stream )
        SCAI_LOG_INFO( logger, "sparseGEVM<" << TypeTraits<ValueType>::id() << "> synchronized" )
    }

    if ( useTexture )
    {
        if ( !syncToken )
        {
            vectorUnbindTexture( x );
        }
        else
        {
            // get routine with the right signature
            void ( *unbind ) ( const ValueType* ) = &vectorUnbindTexture;
            // delay unbind until synchroniziaton
            syncToken->pushRoutine( common::bind( unbind, x ) );
        }
    }
}

/* --------------------------------------------------------------------------- */
/*                          Jacobi                                             */
/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void csr_jacobi_kernel(
    const IndexType* const csrIA,
    const IndexType* const csrJA,
    const ValueType* const csrValues,
    const IndexType numRows,
    const ValueType* const rhs,
    ValueType* const solution,
    const ValueType* const oldSolution,
    const ValueType omega )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType temp = rhs[i];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        const ValueType diag = csrValues[rowStart];

        for ( IndexType jj = rowStart + 1; jj < rowEnd; ++jj )
        {
            temp -= csrValues[jj] * fetchVectorX<ValueType, useTexture>( oldSolution, csrJA[jj] );
        }

        if ( omega == 0.5 )
        {
            solution[i] = omega * ( fetchVectorX<ValueType, useTexture>( oldSolution, i ) + temp / diag );
        }
        else if ( omega == 1.0 )
        {
            solution[i] = temp / diag;
        }
        else
        {
            solution[i] = omega * ( temp / diag ) + ( 1.0 - omega ) * fetchVectorX<ValueType, useTexture>( oldSolution, i );
        }
    }
}

template<typename ValueType>
__inline__ __device__ ValueType getSharedValue( ValueType* shared, const ValueType* const value, const IndexType index )
{
    if ( index / blockDim.x == blockIdx.x )
    {
        return shared[index % blockDim.x];
    }
    else
    {
        return value[index];
    }
}

//these templates allow to combine dynamic shared memory with templates
template<typename ValueType>
struct SharedMemory
{
    //! @brief Return a pointer to the runtime-sized shared memory array.
    //! @returns Pointer to runtime-sized shared memory array
    __device__
    ValueType* getPointer()
    {
        extern __device__ void Error_UnsupportedType(); // Ensure that we won't compile any un-specialized types
        Error_UnsupportedType();
        return ( ValueType* ) 0;
    }

};

template<>
struct SharedMemory<float>
{
    __device__
    float* getPointer()
    {
        extern __shared__ float s_float[];
        return s_float;
    }
};

template<>
struct SharedMemory<double>
{
    __device__
    double* getPointer()
    {
        extern __shared__ double s_double[];
        return s_double;
    }
};

//this is just like the other jacobi kernel, but it performs a coalesced prefetch of the old solution
//instead of using the texture memory
template<typename ValueType>
__global__ void csr_alternate_jacobi_kernel(
    const IndexType* const csrIA,
    const IndexType* const csrJA,
    const ValueType* const csrValues,
    const IndexType numRows,
    const ValueType* const rhs,
    ValueType* const solution,
    const ValueType* const oldSolution,
    const ValueType omega )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    SharedMemory<ValueType> smem;
    ValueType* shared = smem.getPointer();

    if ( i < numRows )
    {
        //this is the prefetch
        shared[threadIdx.x] = oldSolution[i];
        __syncthreads();
        ValueType temp = rhs[i];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd = csrIA[i + 1];
        const ValueType diag = csrValues[rowStart];

        for ( IndexType jj = rowStart + 1; jj < rowEnd; ++jj )
        {
            temp -= csrValues[jj] * getSharedValue<ValueType>( shared, oldSolution, csrJA[jj] );
        }

        if ( omega == 0.5 )
        {
            solution[i] = omega * ( getSharedValue<ValueType>( shared, oldSolution, i ) + temp / diag );
        }
        else if ( omega == 1.0 )
        {
            solution[i] = temp / diag;
        }
        else
        {
            solution[i] = omega * ( temp / diag ) + ( 1.0 - omega ) * getSharedValue<ValueType>( shared, oldSolution, i );
        }
    }
}

template<typename ValueType>
void CUDACSRUtils::jacobi(
    ValueType* const solution,
    const IndexType* const csrIA,
    const IndexType* const csrJA,
    const ValueType* const csrValues,
    const ValueType* const oldSolution,
    const ValueType* const rhs,
    const ValueType omega,
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger, "jacobi, #rows = " << numRows )
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;
    bool useTexture = CUDASettings::useTexture();
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
        useTexture = false; // not yet supported
    }

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    SCAI_LOG_INFO( logger, "Start csr_jacobi_kernel<" << TypeTraits<ValueType>::id()
                   << ", useTexture = " << useTexture << ">" );

    if ( useTexture )
    {
        vectorBindTexture( oldSolution );
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobi_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
        csr_jacobi_kernel <ValueType, true> <<< dimGrid, dimBlock, 0, stream>>>( csrIA, csrJA, csrValues, numRows,
                rhs, solution, oldSolution, omega );
    }
    else
    {
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobi_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
        csr_jacobi_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>>( csrIA, csrJA, csrValues, numRows, rhs,
                solution, oldSolution, omega );
    }

    if ( !syncToken )
    {
        cudaStreamSynchronize( stream );
    }

    if ( useTexture )
    {
        vectorUnbindTexture( oldSolution );
    }
}

/* --------------------------------------------------------------------------- */
/*                          Jacobi halo                                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void csr_jacobiHalo_kernel(
    ValueType* const solution,
    const IndexType* const localIA,
    const ValueType* const localValues,
    const IndexType* const haloIA,
    const IndexType* const haloJA,
    const ValueType* const haloValues,
    const IndexType* const rowIndexes,
    const IndexType numNonEmptyRows,
    const ValueType* const oldSolution,
    const ValueType omega )
{
    const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( ii < numNonEmptyRows )
    {
        IndexType i = ii; // default: rowIndexes is identity

        if ( rowIndexes )
        {
            i = rowIndexes[ii];
        }

        ValueType temp = 0.0;
        const IndexType rowStart = haloIA[i];
        const IndexType rowEnd = haloIA[i + 1];

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            temp += haloValues[jj] * fetchVectorX<ValueType, useTexture>( oldSolution, haloJA[jj] );
        }

        const ValueType diag = localValues[localIA[i]];
        solution[i] -= temp * ( omega / diag );
    }
}

template<typename ValueType>
void CUDACSRUtils::jacobiHalo(
    ValueType solution[],
    const IndexType localIA[],
    const ValueType localValues[],
    const IndexType haloIA[],
    const IndexType haloJA[],
    const ValueType haloValues[],
    const IndexType haloRowIndexes[],
    const ValueType oldSolution[],
    const ValueType omega,
    const IndexType numNonEmptyRows )
{
    SCAI_LOG_INFO( logger, "jacobiHalo, #non-empty rows = " << numNonEmptyRows )
    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numNonEmptyRows, dimBlock.x );
    bool useTexture = CUDASettings::useTexture();
    useTexture = false;

    if ( useTexture )
    {
        vectorBindTexture( oldSolution );
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHalo_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
        csr_jacobiHalo_kernel <ValueType, true> <<< dimGrid, dimBlock>>>( solution, localIA, localValues, haloIA,
                haloJA, haloValues, haloRowIndexes,
                numNonEmptyRows, oldSolution, omega );
    }
    else
    {
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHalo_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
        csr_jacobiHalo_kernel<ValueType, false> <<< dimGrid, dimBlock>>>( solution, localIA, localValues, haloIA,
                haloJA, haloValues, haloRowIndexes, numNonEmptyRows,
                oldSolution, omega );
    }

    SCAI_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED" )
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED" )

    if ( useTexture )
    {
        vectorUnbindTexture( oldSolution );
    }
}

/* --------------------------------------------------------------------------- */
/*                          Jacobi halo with diagonal array                    */
/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void csr_jacobiHaloWithDiag_kernel(
    ValueType* const solution,
    const ValueType* const localDiagValues,
    const IndexType* const haloIA,
    const IndexType* const haloJA,
    const ValueType* const haloValues,
    const IndexType* const rowIndexes,
    const IndexType numNonEmptyRows,
    const ValueType* const oldSolution,
    const ValueType omega )
{
    const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( ii < numNonEmptyRows )
    {
        IndexType i = ii; // default: rowIndexes is identity

        if ( rowIndexes )
        {
            i = rowIndexes[ii];
        }

        ValueType temp = 0.0;
        const IndexType rowStart = haloIA[i];
        const IndexType rowEnd = haloIA[i + 1];

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            temp += haloValues[jj] * fetchVectorX<ValueType, useTexture>( oldSolution, haloJA[jj] );
        }

        const ValueType diag = localDiagValues[i];
        solution[i] -= temp * ( omega / diag );
    }
}

template<typename ValueType>
void CUDACSRUtils::jacobiHaloWithDiag(
    ValueType solution[],
    const ValueType localDiagValues[],
    const IndexType haloIA[],
    const IndexType haloJA[],
    const ValueType haloValues[],
    const IndexType haloRowIndexes[],
    const ValueType oldSolution[],
    const ValueType omega,
    const IndexType numNonEmptyRows )
{
    SCAI_LOG_INFO( logger, "jacobiHaloWithDiag, #non-empty rows = " << numNonEmptyRows )
    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numNonEmptyRows, dimBlock.x );
    bool useTexture = CUDASettings::useTexture();
    useTexture = false;

    if ( useTexture )
    {
        vectorBindTexture( oldSolution );
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHaloWithDiag_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }
    else
    {
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHaloWithDiag_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }

    if ( useTexture )
    {
        csr_jacobiHaloWithDiag_kernel <ValueType, true> <<< dimGrid, dimBlock>>>( solution, localDiagValues, haloIA,
                haloJA, haloValues, haloRowIndexes,
                numNonEmptyRows, oldSolution, omega );
    }
    else
    {
        csr_jacobiHaloWithDiag_kernel<ValueType, false> <<< dimGrid, dimBlock>>>( solution, localDiagValues, haloIA,
                haloJA, haloValues, haloRowIndexes, numNonEmptyRows,
                oldSolution, omega );
    }

    SCAI_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_CSRJACOBIHALOWITHDIAG_CUDAKERNEL_FAILED" )
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "LAMA_STATUS_CSRJACOBIHALOWITHDIAG_CUDAKERNEL_FAILED" )

    if ( useTexture )
    {
        vectorUnbindTexture( oldSolution );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             helper                                                                 */
/* ------------------------------------------------------------------------------------------------------------------ */

__device__ __inline__ IndexType multHlp_getNumActiveThreads(
    IndexType aColIt,
    IndexType aColEnd,
    const IndexType* aIA,
    IndexType aRowIt,
    IndexType offset )
{
#ifdef CUDA_CAP_20
    IndexType end = __popc ( __ballot ( aColIt < aColEnd ) );
#else
    IndexType aColStart = aIA[aRowIt] + offset;
    IndexType end = ( aColEnd - aColStart );
#endif
    return end;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixAddSizes                                                         */
/* ------------------------------------------------------------------------------------------------------------------ */

template<int nWarps>
__global__ void matrixAddSizesKernel(
    IndexType* cIa,
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const IndexType* aIa,
    const IndexType* aJa,
    const IndexType* bIa,
    const IndexType* bJa )
{
    __shared__ volatile IndexType sColA[nWarps];
    __shared__ volatile bool sFound[nWarps];
    IndexType localWarpId = threadIdx.x / warpSize;
    IndexType globalWarpId = ( blockIdx.x * blockDim.x + threadIdx.x ) / warpSize;
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;
//IndexType numWarpsLocal  = blockDim.x / warpSize;
    IndexType numWarpsGlobal = ( blockDim.x * gridDim.x ) / warpSize;
    IndexType rowIt = globalWarpId;

    for ( ; __any( rowIt < numRows ); rowIt += numWarpsGlobal )
    {
        if ( rowIt < numRows )
        {
            if ( diagonalProperty && rowIt >= numColumns )
            {
                diagonalProperty = false;
            }

            IndexType aColIt = aIa[rowIt] + laneId;
            IndexType aColEnd = aIa[rowIt + 1];
            IndexType bColIt = bIa[rowIt] + laneId;
            IndexType bColEnd = bIa[rowIt + 1];

            if ( laneId == 0 )
            {
                cIa[rowIt] = bColEnd - bColIt;
            }

            for ( IndexType aColItOffset = 0; __any( aColIt < aColEnd ); aColIt += warpSize, aColItOffset += warpSize )
            {
                IndexType colA = aColIt < aColEnd ? aJa[aColIt] : cudaNIndex;
                IndexType end = multHlp_getNumActiveThreads( aColIt, aColEnd, aIa, rowIt, aColItOffset );

                for ( IndexType k = 0; k < end && k < warpSize; k++ )
                {
                    if ( laneId == k )
                    {
                        sColA[localWarpId] = colA;
                    }

                    sFound[localWarpId] = false;

                    for ( IndexType bColItOffset = 0; !sFound[localWarpId] && __any( ( bColIt + bColItOffset ) < bColEnd );
                            bColItOffset += warpSize )
                    {
                        IndexType colB = ( bColIt + bColItOffset ) < bColEnd ? bJa[bColIt + bColItOffset] : cudaNIndex;

                        if ( sColA[localWarpId] == colB )
                        {
                            sFound[localWarpId] = true;
                        }
                    }

                    if ( laneId == 0 && !sFound[localWarpId] )
                    {
                        cIa[rowIt]++;
                    }
                }
            }
        }
    }
}

IndexType CUDACSRUtils::matrixAddSizes(
    IndexType cIa[],
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const IndexType aIa[],
    const IndexType aJa[],
    const IndexType bIa[],
    const IndexType bJa[] )
{
    SCAI_REGION( "CUDA.CSRUtils.matrixAddSizes" )
    SCAI_LOG_INFO(
        logger,
        "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )
    SCAI_CHECK_CUDA_ACCESS
// Reset cIa
    thrust::device_ptr<IndexType> cIaPtr( cIa );
    thrust::fill( cIaPtr, cIaPtr + numRows, 0 );
// TODO: Check if diagonal property needs special attention
    matrixAddSizesKernel<NUM_WARPS> <<< NUM_BLOCKS, NUM_THREADS>>>( cIa, numRows, numColumns, diagonalProperty,
            aIa, aJa, bIa, bJa );
    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
// Convert sizes array to offset array
    thrust::exclusive_scan( cIaPtr, cIaPtr + numRows + 1, cIaPtr );
// Copy numValues from cIa to Host
// TODO: use cuMem cpy
    thrust::device_ptr<IndexType> iaPtr( cIa );
    thrust::host_vector<IndexType> numValues( iaPtr + numRows, iaPtr + numRows + 1 );
    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
// TODO: write it!
    return numValues[0];
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             hashTable Methods                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

__device__
inline bool multHlp_insertIndexex( IndexType colB,
                                   IndexType sHashTableIndexes[],
                                   IndexType aRowIt,
                                   IndexType* chunkPtr,
                                   volatile IndexType chunkList[],
                                   IndexType numReservedChunks,
                                   IndexType* cIA )
{
    const IndexType one = 1;

    IndexType fx = HASH_A * colB;
    IndexType gx = ( fx + HASH_B ) % HASH_P;

    if ( numReservedChunks == 0 )
    {
        for ( IndexType i = 0; i < NUM_HASH_RETRIES; i++ )
        {
            IndexType hash = ( gx + HASH_C0 * i + HASH_C1 *  i * i ) % NUM_ELEMENTS_IN_SHARED;
            IndexType val = common::CUDAUtils::atomicCAS( &sHashTableIndexes[hash], cudaNIndex, colB );

            if ( val == cudaNIndex )
            {
                common::CUDAUtils::atomicAdd( &cIA[aRowIt], one );
                return true;
            }

            if ( val == colB )
            {
                return true;
            }
        }

        return false;
    }

    for ( IndexType i = 0; i < NUM_HASH_RETRIES; i++ )
    {
        IndexType globalHash = ( gx + HASH_C0 * i + HASH_C1 * ( IndexType ) i * i ) % ( NUM_ELEMENTS_PER_CHUNK * numReservedChunks );
        IndexType localHash = globalHash % NUM_ELEMENTS_PER_CHUNK;
        IndexType chunk = globalHash / NUM_ELEMENTS_PER_CHUNK;
        IndexType val = common::CUDAUtils::atomicCAS( &chunkPtr[chunkList[chunk] * NUM_ELEMENTS_PER_CHUNK + localHash], cudaNIndex, colB );

        if ( val == cudaNIndex )
        {
            common::CUDAUtils::atomicAdd( &cIA[aRowIt], one );
            return true;
        }

        if ( val == colB )
        {
            return true;
        }
    }

    return false;
}

template <typename ValueType>
__device__
inline bool multHlp_insertValues( IndexType colB,
                                  IndexType* sHashTableIndexes,
                                  ValueType* sHashTableValues,
                                  IndexType* indexChunks,
                                  ValueType* valueChunks,
                                  volatile IndexType chunkList[],
                                  IndexType numReservedChunks,
                                  ValueType valB,
                                  ValueType sValA )
{
    IndexType fx = HASH_A * colB;
    IndexType gx = ( fx + HASH_B ) % HASH_P;

    if ( numReservedChunks == 0 )
    {
        for ( IndexType i = 0; i < NUM_HASH_RETRIES; i++ )
        {
            IndexType hash = ( gx + HASH_C0 * i + HASH_C1 * i * i ) % NUM_ELEMENTS_IN_SHARED;
            IndexType val = common::CUDAUtils::atomicCAS( &sHashTableIndexes[hash], cudaNIndex, colB );

            if ( val == cudaNIndex )
            {
                sHashTableValues[hash] = valB * sValA;
                return true;
            }

            if ( val == colB )
            {
                sHashTableValues[hash] += valB * sValA;
                return true;
            }
        }

        return false;
    }

    for ( IndexType i = 0; i < NUM_HASH_RETRIES; i++ )
    {
        IndexType globalHash = ( gx + HASH_C0 * i + HASH_C1 * ( IndexType ) i * i ) % ( NUM_ELEMENTS_PER_CHUNK * numReservedChunks );
        IndexType localHash = globalHash % NUM_ELEMENTS_PER_CHUNK;
        IndexType chunk = globalHash / NUM_ELEMENTS_PER_CHUNK;
        IndexType val = common::CUDAUtils::atomicCAS( &indexChunks[chunkList[chunk] * NUM_ELEMENTS_PER_CHUNK + localHash], cudaNIndex, colB );

        if ( val == cudaNIndex )
        {
            valueChunks[chunkList[chunk] * NUM_ELEMENTS_PER_CHUNK + localHash] = sValA * valB;
            return true;
        }

        if ( val == colB )
        {
            valueChunks[chunkList[chunk] * NUM_ELEMENTS_PER_CHUNK + localHash] += sValA * valB;
            return true;
        }
    }

    return false;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixMultiplySizes                                                    */
/* ------------------------------------------------------------------------------------------------------------------ */

__device__
inline bool multHlp_nextRow( IndexType* row,
                             IndexType numRows
#ifdef USE_LOAD_BALANCING
                             , IndexType* rowCounter
#endif
                           )
{
#ifdef USE_LOAD_BALANCING
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;
    IndexType localWarpId = threadIdx.x / warpSize;
    __shared__ volatile IndexType sRowIt[NUM_WARPS];

    if ( laneId == 0 )
    {
        IndexType one = 1;
        sRowIt[localWarpId] = common::CUDAUtils::atomicAdd( rowCounter, one );
    }

    *row = sRowIt[localWarpId];

    if ( *row < numRows )
    {
        return true;
    }
    else
    {
        return false;
    }

#else
    IndexType numWarpsGlobal = ( blockDim.x * gridDim.x ) / warpSize;
    *row += numWarpsGlobal;

    if ( *row < numRows )
    {
        return true;
    }
    else
    {
        return false;
    }

#endif
}

__device__
inline void multHlp_releaseChunks ( IndexType* chunkList,
                                    volatile IndexType* sChunkList,
                                    volatile IndexType* sReservedChunks,
                                    IndexType chunkCount )
{
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;

    if ( laneId == 0 )
    {
        for ( IndexType i = *sReservedChunks - 1; i >= *sReservedChunks - chunkCount; --i )
        {
            IndexType headItem;
            IndexType old;

            do
            {
                headItem = chunkList[0];
                chunkList[sChunkList[i] + 1] = headItem;
                old = common::CUDAUtils::atomicCAS( const_cast<IndexType*>( &chunkList[0] ), headItem, sChunkList[i] );
            }
            while ( old != headItem );
        }
    }

    *sReservedChunks = *sReservedChunks - chunkCount;
}

__device__
inline bool multHlp_reserveChunks( IndexType* chunkList,
                                   volatile IndexType* sChunkList,
                                   volatile IndexType* sReservedChunks,
                                   IndexType chunkCount )
{
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;

    if ( chunkCount > NUM_CHUNKS_PER_WARP )
    {
//        printf("to many chunks %i\n", chunkCount);
        return false;
    }

    if ( laneId == 0 && chunkCount > 0 && *sReservedChunks != chunkCount )
    {
        if ( *sReservedChunks < chunkCount )
        {
            for ( IndexType i = *sReservedChunks; i < chunkCount; ++i )
            {
                IndexType headItem;
                IndexType nextItem;
                IndexType old;

                do
                {
                    headItem = chunkList[0];

                    if ( headItem != cudaNIndex )
                    {
                        __threadfence();
                        nextItem = chunkList[headItem + 1];

                        old = common::CUDAUtils::atomicCAS( const_cast<IndexType*>( &chunkList[0] ), headItem, nextItem );

                        if ( old == headItem )
                        {
                            sChunkList[i] = headItem;
                        }
                    }
                    else
                    {
//                        printf("no more chunks!\n");
                        return false;
                    }
                }
                while ( old != headItem );
            }

            *sReservedChunks = chunkCount;
            return true;
        }
        else
        {
            multHlp_releaseChunks ( chunkList, sChunkList, sReservedChunks, *sReservedChunks - chunkCount );
            return true;
        }
    }
    else
    {
        return true;
    }
}

__device__
inline void multHlp_initializeChunks ( IndexType* sHashTable,
                                       IndexType* chunks,
                                       const IndexType numElementsPerChunk,
                                       volatile IndexType* sChunkList,
                                       volatile IndexType sReservedChunks )
{
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;

    if ( sReservedChunks == 0 )
    {
        for ( IndexType i = 0; i < NUM_ELEMENTS_IN_SHARED; i += warpSize )
        {
            if ( i + laneId < NUM_ELEMENTS_IN_SHARED )
            {
                sHashTable[i + laneId] = cudaNIndex;
            }
        }

        return;
    }

    for ( IndexType i = 0; i < sReservedChunks; ++i )
    {
        IndexType chunkId = sChunkList[i];

        for ( IndexType j = laneId; j < numElementsPerChunk; j += warpSize )
        {
            chunks[chunkId * numElementsPerChunk + j] = cudaNIndex;
        }
    }
}

__device__
inline IndexType multHlp_growth ( IndexType numChunks )
{
    if ( numChunks == 0 )
    {
        return 2;
    }
    else
    {
        return numChunks * 2;
    }
}

__device__
inline IndexType multHlp_calcOptChunkCount ( IndexType row,
        const IndexType* cIA,
        const IndexType numElementsPerChunk )
{
    IndexType numElements = cIA[row + 1] - cIA[row];

    if ( numElements * 2 < NUM_ELEMENTS_IN_SHARED )
    {
        return 0;
    }
    else
    {
        return ( ( ( cIA[row + 1] - cIA[row] ) * 2 ) / numElementsPerChunk ) + 1;
    }
}

__global__
void matrixMultiplySizesKernel(
    const IndexType* aIA,
    const IndexType* aJA,
    const IndexType* bIA,
    const IndexType* bJA,
    IndexType* cIA,
    const IndexType numRows,
    const IndexType numColumns,
    IndexType* chunkPtr,
    IndexType* chunkList,
    IndexType numChunks,
    bool* hashError,
    bool diagonalProperty )
{
    __shared__ IndexType sHashTable[NUM_ELEMENTS_IN_SHARED];
    __shared__ volatile IndexType sReservedChunks;
    __shared__ volatile IndexType sChunkList[NUM_CHUNKS_PER_WARP];
    __shared__ volatile IndexType sColA;
    __shared__ volatile IndexType sRowIt;
    __shared__ volatile bool sInsertMiss;
    IndexType globalWarpId = ( blockIdx.x * blockDim.x + threadIdx.x ) / warpSize;
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;
    IndexType colB;
    IndexType aRowIt = globalWarpId;
    bool localSystemError = false;
    sReservedChunks = 0;

    if ( aRowIt < numRows )
    {
        do
        {
            do
            {
                sInsertMiss = false;
                IndexType aColIt = aIA[aRowIt] + laneId;
                IndexType aColEnd = aIA[aRowIt + 1];

                if ( laneId == 0 && diagonalProperty )
                {
                    cIA[aRowIt]++;
                }

                multHlp_initializeChunks( sHashTable,
                                          chunkPtr,
                                          NUM_ELEMENTS_PER_CHUNK,
                                          sChunkList,
                                          sReservedChunks );

                for ( IndexType offset = 0; __any( aColIt < aColEnd ); aColIt += warpSize, offset += warpSize )
                {
                    IndexType colA = aColIt < aColEnd ? aJA[aColIt] : cudaNIndex;
                    IndexType end = multHlp_getNumActiveThreads( aColIt, aColEnd, aIA, aRowIt, offset );

                    for ( IndexType k = 0; k < end && k < warpSize; k++ )
                    {
                        if ( laneId == k )
                        {
                            sColA = colA;
                        }

                        IndexType bColIt = bIA[sColA] + laneId;
                        IndexType bColEnd = bIA[sColA + 1];

                        for ( ; __any( bColIt < bColEnd ); bColIt += warpSize )
                        {
                            colB = bColIt < bColEnd ? bJA[bColIt] : cudaNIndex;

                            if ( colB != cudaNIndex && ( !diagonalProperty || colB != aRowIt ) )
                            {
                                bool inserted = multHlp_insertIndexex( colB,
                                                                       sHashTable,
                                                                       aRowIt,
                                                                       chunkPtr,
                                                                       sChunkList,
                                                                       sReservedChunks,
                                                                       cIA );

                                if ( !inserted )
                                {
                                    sInsertMiss = true;
                                }
                            }
                        }
                    }
                }

                // only release if insertion was ok, otherwire reserve some more
                // STEP x: release reserved chunks
                if ( laneId == 0 )
                {
                    if ( sInsertMiss )
                    {
                        cIA[aRowIt] = 0;

                        if ( !multHlp_reserveChunks( chunkList, sChunkList, &sReservedChunks, multHlp_growth( sReservedChunks ) ) )
                        {
                            // ABORT KERNEL HERE;
                            localSystemError = true;
                        }
                    }
                }

                if ( __any( localSystemError ) )
                {
                    *hashError = true;
                    return;
                }
            }
            while ( sInsertMiss );
        }
        while ( multHlp_nextRow( &aRowIt, numRows ) );
    }

    // release all remaining chunks
    multHlp_releaseChunks( chunkList, sChunkList, &sReservedChunks, sReservedChunks );
}

struct multHlp_chunkFill
{
    const IndexType n;
    multHlp_chunkFill( IndexType _n )
        : n( _n )
    {
    }
    __device__
    IndexType operator()( IndexType i )
    {
        if ( i == ( n - 1 ) )
        {
            return cudaNIndex;
        }

        return i;
    }
};

IndexType CUDACSRUtils::matrixMultiplySizes(
    IndexType cIa[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType /* k */,
    bool diagonalProperty,
    const IndexType aIa[],
    const IndexType aJa[],
    const IndexType bIa[],
    const IndexType bJa[] )
{
    SCAI_REGION( "CUDA.CSR.matrixMultiplySizes" )
    SCAI_LOG_INFO(
        logger,
        "matrixMutliplySizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )
    SCAI_CHECK_CUDA_ACCESS
    // Reset cIa
    thrust::device_ptr<IndexType> cIaPtr( cIa );
    thrust::fill( cIaPtr, cIaPtr + numRows, 0 );
    ContextPtr loc = Context::getContextPtr( context::CUDA );
    MemoryPtr mem = loc->getMemoryPtr();
    bool hashErrorHost = false;
    bool* hashError = ( bool* ) mem->allocate( sizeof( bool ) );
    cudaMemcpy( hashError, &hashErrorHost, sizeof( bool ), cudaMemcpyHostToDevice );
    size_t free;
    size_t total;
    cuMemGetInfo( &free, &total );
    int nnz_a;
    int nnz_b;
    cudaMemcpy( &nnz_a, &aIa[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );
    cudaMemcpy( &nnz_b, &bIa[numColumns], sizeof( IndexType ), cudaMemcpyDeviceToHost );
    int avgDensity = ( nnz_a / numRows + nnz_b / numColumns ) / 2;
    int numChunks;
    int maxNumChunks = ( free - ( 100 * 1024 * 1024 ) ) / ( NUM_ELEMENTS_PER_CHUNK * sizeof ( IndexType ) * 2 );
    int chunksPerWarp = NUM_BLOCKS * ( ( avgDensity * 8 ) / NUM_ELEMENTS_PER_CHUNK + 1 );

    if ( chunksPerWarp > maxNumChunks )
    {
        numChunks = maxNumChunks;
    }
    else
    {
        numChunks = chunksPerWarp;
    }

    unsigned int hashTableAllocatedBytes = numChunks * NUM_ELEMENTS_PER_CHUNK * sizeof( IndexType );
    IndexType* hashTable = ( IndexType* ) mem->allocate( hashTableAllocatedBytes );
    // chunkList table needs one integers per chunk plus 1 start pointer
    unsigned int chunkListAllocatedBytes = numChunks * sizeof( IndexType ) + sizeof( IndexType );
    IndexType* chunkList = ( IndexType* ) mem->allocate( chunkListAllocatedBytes );
    thrust::device_ptr<IndexType> chunkListPtr( chunkList );
    thrust::transform( thrust::make_counting_iterator( 0 ),
                       thrust::make_counting_iterator( numChunks + 1 ),
                       chunkListPtr,
                       multHlp_chunkFill( numChunks + 1 ) );
    matrixMultiplySizesKernel <<< NUM_BLOCKS, NUM_THREADS>>>( aIa,
            aJa,
            bIa,
            bJa,
            cIa,
            numRows,
            numColumns,
            hashTable,
            chunkList,
            numChunks,
            hashError,
            diagonalProperty );
    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
    cudaMemcpy( &hashErrorHost, hashError, sizeof( bool ), cudaMemcpyDeviceToHost );

    if ( hashErrorHost )
    {
        COMMON_THROWEXCEPTION( "Multiplication failed!" );
    }

    // Free hashTable and hashError
    mem->free( ( void* ) hashError, sizeof( bool ) );
    mem->free( ( void* ) hashTable, hashTableAllocatedBytes );
    mem->free( ( void* ) chunkList, chunkListAllocatedBytes );
    // Convert sizes array to offset array
    thrust::exclusive_scan( cIaPtr, cIaPtr + numRows + 1, cIaPtr );
    IndexType numValues;
    cudaMemcpy( &numValues, &cIa[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );
    return numValues;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixAdd                                                              */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, int nWarps>
__global__
void matrixAddKernel(
    IndexType* cJA,
    ValueType* cValues,
    const IndexType* cIA,
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const ValueType alpha,
    const IndexType* aIA,
    const IndexType* aJA,
    const ValueType* aValues,
    const ValueType beta,
    const IndexType* bIA,
    const IndexType* bJA,
    const ValueType* bValues )
{
// TODO: Just naive implementation, could be done faster, but works!
// TODO: Check if diagonal property needs special attention
    __shared__ volatile IndexType sColA[nWarps];
    __shared__ volatile ValueType sValA[nWarps];
    __shared__ volatile IndexType sFoundJa[nWarps];
    IndexType localWarpId = threadIdx.x / warpSize;
    IndexType globalWarpId = ( blockIdx.x * blockDim.x + threadIdx.x ) / warpSize;
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;
//IndexType numWarpsLocal  = blockDim.x / warpSize;
    IndexType numWarpsGlobal = ( blockDim.x * gridDim.x ) / warpSize;
    IndexType rowIt = globalWarpId;

    for ( ; __any( rowIt < numRows ); rowIt += numWarpsGlobal )
    {
        if ( rowIt < numRows )
        {
            if ( diagonalProperty && rowIt >= numColumns )
            {
                diagonalProperty = false;
            }

            IndexType aColIt = aIA[rowIt] + laneId;
            IndexType aColEnd = aIA[rowIt + 1];
            IndexType bColIt = bIA[rowIt] + laneId;
            IndexType bColEnd = bIA[rowIt + 1];
            IndexType cColIt = cIA[rowIt] + laneId;

// Copy values of b to C
            for ( IndexType bColOffset = 0; __any( ( bColIt + bColOffset ) < bColEnd ); bColOffset += warpSize )
            {
                IndexType colB = ( bColIt + bColOffset ) < bColEnd ? bJA[bColIt + bColOffset] : cudaNIndex;
                ValueType valB = ( bColIt + bColOffset ) < bColEnd ? bValues[bColIt + bColOffset] : static_cast<ValueType>( 0 );

                if ( colB != cudaNIndex )
                {
                    cJA[cColIt + bColOffset] = colB;
                    cValues[cColIt + bColOffset] = valB * beta;
                }
            }

// Offset in c after coping b to c
            IndexType cColOffset = bIA[rowIt + 1] - bIA[rowIt];

// Add values of a to c
            for ( IndexType aColItOffset = 0; __any( aColIt < aColEnd ); aColIt += warpSize, aColItOffset += warpSize )
            {
                IndexType colA = aColIt < aColEnd ? aJA[aColIt] : cudaNIndex;
                ValueType valA = aColIt < aColEnd ? aValues[aColIt] : static_cast<ValueType>( 0 );
                IndexType end = multHlp_getNumActiveThreads( aColIt, aColEnd, aIA, rowIt, aColItOffset );

                for ( IndexType k = 0; k < end && k < warpSize; k++ )
                {
                    if ( laneId == k )
                    {
                        sColA[localWarpId] = colA;
                        sValA[localWarpId] = valA;
                        sFoundJa[localWarpId] = cudaNIndex;
                    }

                    for ( IndexType bColItOffset = 0; ( sFoundJa[localWarpId] == cudaNIndex ) && __any( ( bColIt + bColItOffset ) < bColEnd );
                            bColItOffset += warpSize )
                    {
                        IndexType colB = ( bColIt + bColItOffset ) < bColEnd ? bJA[bColIt + bColItOffset] : cudaNIndex;

                        if ( sColA[localWarpId] == colB )
                        {
                            sFoundJa[localWarpId] = laneId + bColItOffset;
                        }
                    }

                    if ( laneId == 0 )
                    {
                        if ( sFoundJa[localWarpId] == cudaNIndex )
                        {
                            // Element is new element, add new element
                            cJA[cColIt + cColOffset] = colA;
                            cValues[cColIt + cColOffset] = sValA[localWarpId] * alpha;
                            cColOffset++;
                        }
                        else
                        {
                            // Element exists, add values
                            // We can use cColIt, because this is thread with laneId = 0!
                            cValues[cColIt + sFoundJa[localWarpId]] += sValA[localWarpId] * alpha;
                        }
                    }
                }
            }
        }
    }
}

template<typename ValueType>
void CUDACSRUtils::matrixAdd(
    IndexType cJA[],
    ValueType cValues[],
    const IndexType cIA[],
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const ValueType alpha,
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const ValueType beta,
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[] )
{
    SCAI_REGION( "CUDA.CSRUtils.matrixAdd" )
    SCAI_LOG_INFO( logger, "matrixAdd for " << numRows << "x" << numColumns << " matrix" )
    SCAI_CHECK_CUDA_ACCESS
    matrixAddKernel<ValueType, NUM_WARPS> <<< NUM_BLOCKS, NUM_THREADS>>>( cJA, cValues, cIA, numRows, numColumns,
            diagonalProperty, alpha, aIA, aJA, aValues, beta, bIA, bJA, bValues );
    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
}
/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixMultiply                                                         */
/* ------------------------------------------------------------------------------------------------------------------ */

template <typename ValueType>
__device__
inline void multHlp_copyHashtable ( volatile IndexType* sColA,
                                    const IndexType* cIA,
                                    IndexType laneId,
                                    IndexType aRowIt,
                                    const ValueType alpha,
                                    IndexType* cJA,
                                    ValueType* cValues,
                                    IndexType* sHashTableIndexes,
                                    ValueType* sHashTableValues,
                                    IndexType* indexChunks,
                                    ValueType* valueChunks,
                                    volatile IndexType chunkList[],
                                    IndexType numReservedChunks,
                                    bool diagonalProperty,
                                    ValueType diagonalElement )

{
    // TODO: rename sColA => destinationOffset!

    *sColA = 0;
    IndexType rowOffset = cIA[aRowIt];
    IndexType one = 1;
    IndexType hashCol;
    ValueType hashVal;

    if ( diagonalProperty && laneId == 0 )
    {
        cJA[rowOffset] = aRowIt;
        cValues[rowOffset] = diagonalElement * alpha;
        *sColA = 1;
    }

    if ( numReservedChunks == 0 )
    {
        for ( int j = laneId; j < NUM_ELEMENTS_IN_SHARED; j += warpSize )
        {
            hashCol = sHashTableIndexes[j];
            hashVal = sHashTableValues[j];
#if SCAI_CUDA_COMPUTE_CAPABILITY >= 20
            IndexType localOffset;
            // TODO: be carefull here, ballot is warpsize Bit's long!
            IndexType ballot = __ballot ( hashCol != cudaNIndex );
            localOffset = __popc( ballot << ( warpSize - laneId ) );

            if ( hashCol != cudaNIndex )
            {
                cJA[rowOffset + *sColA + localOffset] = hashCol;
                cValues[rowOffset + *sColA + localOffset] = hashVal * alpha;
            }

            *sColA += __popc( ballot );
#else

            if ( hashCol != cudaNIndex )
            {
                // the volatile attribute must be cast away
                IndexType offset = common::CUDAUtils::atomicAdd( const_cast<IndexType*>( sColA ), one );
                cJA[rowOffset + offset] = hashCol;
                cValues[rowOffset + offset] = hashVal * alpha;
            }

#endif
        }

        return;
    }

    for ( int i = 0; i < numReservedChunks; ++i )
    {
        for ( int j = laneId; j < NUM_ELEMENTS_PER_CHUNK; j += warpSize )
        {
            hashCol = indexChunks[chunkList[i] * NUM_ELEMENTS_PER_CHUNK + j];
            hashVal = valueChunks[chunkList[i] * NUM_ELEMENTS_PER_CHUNK + j];
#if SCAI_CUDA_COMPUTE_CAPABILITY >= 20
            IndexType localOffset;
            // TODO: be carefull here, ballot is warpsize Bit's long!
            IndexType ballot = __ballot ( hashCol != cudaNIndex );
            localOffset = __popc( ballot << ( warpSize - laneId ) );

            if ( hashCol != cudaNIndex )
            {
                cJA[rowOffset + *sColA + localOffset] = hashCol;
                cValues[rowOffset + *sColA + localOffset] = hashVal * alpha;
            }

            if ( laneId == 0 )
            {
                *sColA += __popc( ballot );
            }

#else

            if ( hashCol != cudaNIndex )
            {
                IndexType offset = common::CUDAUtils::atomicAdd( const_cast<IndexType*>( sColA ), IndexType( 1 ) );
                cJA[rowOffset + offset] = hashCol;
                cValues[rowOffset + offset] = hashVal * alpha;
            }

#endif
        }
    }
}

template<typename ValueType>
__global__
void matrixMultiplyKernel(
    const IndexType* aIA,
    const IndexType* aJA,
    const ValueType* aValues,
    const IndexType* bIA,
    const IndexType* bJA,
    const ValueType* bValues,
    const IndexType* cIA,
    const ValueType alpha,
    IndexType* cJA,
    ValueType* cValues,
    const IndexType numRows,
    const IndexType numColumns,
    IndexType* indexChunks,
    ValueType* valueChunks,
    IndexType* chunkList,
    const IndexType numChunks,
    bool* hashError,
    bool diagonalProperty )
{
    __shared__ IndexType sHashTableIndexes[NUM_ELEMENTS_IN_SHARED];
    __shared__ ValueType sHashTableValues[NUM_ELEMENTS_IN_SHARED];
    __shared__ volatile IndexType sReservedChunks;
    __shared__ volatile IndexType sChunkList[NUM_CHUNKS_PER_WARP];
    __shared__ volatile IndexType sColA;
    __shared__ volatile ValueType sValA;
    __shared__ volatile IndexType sRowIt;
    __shared__ volatile bool sInsertMiss;
    __shared__ volatile ValueType diagonalElement;
    IndexType globalWarpId = ( blockIdx.x * blockDim.x + threadIdx.x ) / warpSize;
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;
    IndexType colB;
    IndexType aRowIt = globalWarpId;
    bool localSystemError = false;
    sReservedChunks = 0;

    if ( aRowIt < numRows )
    {
        do
        {
            IndexType optimalChunkCount = multHlp_calcOptChunkCount ( aRowIt, cIA, NUM_ELEMENTS_PER_CHUNK );

            // reserve Chunks
            if ( !multHlp_reserveChunks( chunkList, sChunkList, &sReservedChunks, optimalChunkCount ) )
            {
                // ABORT KERNEL HERE;
                localSystemError = true;
            }

            if ( __any( localSystemError ) )
            {
                *hashError = true;
                return;
            }

            do
            {
                sInsertMiss = false;
                IndexType aColIt = aIA[aRowIt] + laneId;
                IndexType aColEnd = aIA[aRowIt + 1];

                if ( laneId == 0 && diagonalProperty )
                {
                    diagonalElement = 0.0;
                }

                multHlp_initializeChunks( sHashTableIndexes,
                                          indexChunks,
                                          NUM_ELEMENTS_PER_CHUNK,
                                          sChunkList,
                                          sReservedChunks );

                for ( IndexType offset = 0; __any( aColIt < aColEnd ); aColIt += warpSize, offset += warpSize )
                {
                    IndexType colA = aColIt < aColEnd ? aJA[aColIt] : cudaNIndex;
                    ValueType valA = aColIt < aColEnd ? aValues[aColIt] : static_cast<ValueType>( 0 );
                    IndexType end = multHlp_getNumActiveThreads( aColIt, aColEnd, aIA, aRowIt, offset );

                    for ( IndexType k = 0; k < end && k < warpSize; k++ )
                    {
                        if ( laneId == k )
                        {
                            sColA = colA;
                            sValA = valA;
                        }

                        IndexType bColIt = bIA[sColA] + laneId;
                        IndexType bColEnd = bIA[sColA + 1];

                        for ( ; __any( bColIt < bColEnd ); bColIt += warpSize )
                        {
                            colB = bColIt < bColEnd ? bJA[bColIt] : cudaNIndex;
                            ValueType valB = bColIt < bColEnd ? bValues[bColIt] : static_cast<ValueType>( 0 );

                            if ( diagonalProperty && colB == aRowIt )
                            {
                                diagonalElement += sValA * valB;
                            }
                            else
                            {
                                if ( colB != cudaNIndex && ( !diagonalProperty || colB != aRowIt ) )
                                {
                                    bool inserted = multHlp_insertValues( colB,
                                                                          sHashTableIndexes,
                                                                          sHashTableValues,
                                                                          indexChunks,
                                                                          valueChunks,
                                                                          sChunkList,
                                                                          sReservedChunks,
                                                                          valB,
                                                                          sValA );

                                    if ( !inserted )
                                    {
                                        sInsertMiss = true;
                                    }
                                }
                            }
                        }
                    }
                }

                if ( !sInsertMiss )
                {
                    multHlp_copyHashtable ( &sColA,
                                            cIA,
                                            laneId,
                                            aRowIt,
                                            alpha,
                                            cJA,
                                            cValues,
                                            sHashTableIndexes,
                                            sHashTableValues,
                                            indexChunks,
                                            valueChunks,
                                            sChunkList,
                                            sReservedChunks,
                                            diagonalProperty,
                                            diagonalElement );
                }
                else
                {
                    if ( !multHlp_reserveChunks( chunkList, sChunkList, &sReservedChunks, multHlp_growth( sReservedChunks ) ) )
                    {
                        // ABORT KERNEL HERE;
                        localSystemError = true;
                    }

                    if ( __any( localSystemError ) )
                    {
                        *hashError = true;
                        return;
                    }
                }
            }
            while ( sInsertMiss );
        }
        while ( multHlp_nextRow( &aRowIt, numRows ) );
    }

    // release all remaining chunks
    multHlp_releaseChunks( chunkList, sChunkList, &sReservedChunks, sReservedChunks );
}

template<typename ValueType>
void CUDACSRUtils::matrixMultiply(
    const IndexType cIa[],
    IndexType cJa[],
    ValueType cValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType /* k */,
    const ValueType alpha,
    bool diagonalProperty,
    const IndexType aIa[],
    const IndexType aJa[],
    const ValueType aValues[],
    const IndexType bIa[],
    const IndexType bJa[],
    const ValueType bValues[] )
{
    SCAI_REGION( "CUDA.CSRUtils.matrixMultiply" )
    SCAI_LOG_INFO( logger, "matrixMultiply for " << numRows << "x" << numColumns << " matrix" )
    SCAI_CHECK_CUDA_ACCESS
    ContextPtr loc = Context::getContextPtr( context::CUDA );
    MemoryPtr mem = loc->getMemoryPtr();
    bool hashErrorHost = false;
    bool* hashError = ( bool* ) mem->allocate( sizeof( bool ) );
    cudaMemcpy( hashError, &hashErrorHost, sizeof( bool ), cudaMemcpyHostToDevice );
    size_t free;
    size_t total;
    cuMemGetInfo( &free, &total );
    IndexType nnz_a;
    IndexType nnz_b;
    cudaMemcpy( &nnz_a, &aIa[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );
    cudaMemcpy( &nnz_b, &bIa[numColumns], sizeof( IndexType ), cudaMemcpyDeviceToHost );
    IndexType avgDensity = ( nnz_a / numRows + nnz_b / numColumns ) / 2;
    IndexType numChunks;
    IndexType maxNumChunks = ( free - ( 100 * 1024 * 1024 ) ) / ( NUM_ELEMENTS_PER_CHUNK * sizeof ( IndexType ) * 2 );
    IndexType chunksPerWarp = NUM_BLOCKS * ( ( avgDensity * 8 ) / NUM_ELEMENTS_PER_CHUNK + 1 );

    if ( chunksPerWarp > maxNumChunks )
    {
        numChunks = maxNumChunks;
    }
    else
    {
        numChunks = chunksPerWarp;
    }

    unsigned int hashTableAllocatedBytes = numChunks * NUM_ELEMENTS_PER_CHUNK * ( sizeof( IndexType ) + sizeof( ValueType ) );
    void* chunks = ( void* ) mem->allocate( hashTableAllocatedBytes );
    IndexType* indexChunks = ( IndexType* ) chunks;
    ValueType* valueChunks = ( ValueType* ) ( indexChunks + numChunks * NUM_ELEMENTS_PER_CHUNK );
    // chunkList table needs one integers per chunk plus 1 start pointer
    unsigned int chunkListAllocatedBytes = numChunks * sizeof( IndexType ) + sizeof( IndexType );
    IndexType* chunkList = ( IndexType* ) mem->allocate( chunkListAllocatedBytes );
    thrust::device_ptr<IndexType> chunkListPtr( chunkList );
    IndexType zero = 0;
    thrust::transform( thrust::make_counting_iterator( zero ),
                       thrust::make_counting_iterator( numChunks + 1 ),
                       chunkListPtr,
                       multHlp_chunkFill( numChunks + 1 ) );
    matrixMultiplyKernel <<< NUM_BLOCKS, NUM_THREADS>>>( aIa,
            aJa,
            aValues,
            bIa,
            bJa,
            bValues,
            cIa,
            alpha,
            cJa,
            cValues,
            numRows,
            numColumns,
            indexChunks,
            valueChunks,
            chunkList,
            numChunks,
            hashError,
            diagonalProperty );
    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
    cudaMemcpy( &hashErrorHost, hashError, sizeof( bool ), cudaMemcpyDeviceToHost );

    if ( hashErrorHost )
    {
        COMMON_THROWEXCEPTION( "Multiplication failed!" );
    }

    // Free hashTable and hashError
    mem->free( ( void* ) hashError, sizeof( bool ) );
    mem->free( ( void* ) chunks, hashTableAllocatedBytes );
    mem->free( ( void* ) chunkList, chunkListAllocatedBytes );
    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDACSRUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_DEBUG( logger, "set CSR routines for CUDA in Interface" )
    KernelRegistry::set<CSRKernelTrait::sizes2offsets>( sizes2offsets, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::offsets2sizes>( offsets2sizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAddSizes>( matrixAddSizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, ctx, flag );
}

template<typename ValueType>
void CUDACSRUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_DEBUG( logger, "register CSRUtils CUDA-routines for CUDA at kernel registry [" << flag
                     << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<ValueType> >( convertCSR2CSC, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::normalGEVM<ValueType> >( normalGEVM, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sparseGEVM<ValueType> >( sparseGEVM, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAdd<ValueType> >( matrixAdd, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiply<ValueType> >( matrixMultiply, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobiHaloWithDiag<ValueType> >( jacobiHaloWithDiag, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void CUDACSRUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_DEBUG( logger, "register CSRUtils CUDA-routines for CUDA at kernel registry [" << flag
                   << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::scaleRows<ValueType, OtherValueType> >( scaleRows, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDACSRUtils::CUDACSRUtils()
{
    SCAI_LOG_INFO( logger, "register CSRUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_CUDA_LIST, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDACSRUtils::~CUDACSRUtils()
{
    SCAI_LOG_INFO( logger, "unregister CSRUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_CUDA_LIST, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDACSRUtils CUDACSRUtils::guard;    // guard variable for registration

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

unsigned int CUDACSRUtils::lastHashTableSize = 1024;

} /* end namespace sparsekernel */

} /* end namespace scai */
