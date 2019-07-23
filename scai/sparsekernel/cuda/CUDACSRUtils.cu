/**
 * @file CUDACSRUtils.cu
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
#include <scai/utilskernel/cuda/CUDASparseUtils.hpp>

#include <scai/hmemo/Memory.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAUtils.hpp>
#include <scai/common/SCAITypes.hpp>
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
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/reduce.h>

#include <functional>

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

#if __CUDACC_VER_MAJOR__ >= 9
    #define FULL_MASK 0xffffffff
    #define CUDA_ANY( cond ) __any_sync( FULL_MASK, ( cond ) )
#else
    #define CUDA_ANY( cond ) __any( cond )
#endif

using namespace scai::common;
using namespace scai::hmemo;

namespace scai
{

using utilskernel::CUDAUtils;
using utilskernel::CUDASparseUtils;

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
/*     getColumnPositions                                                      */
/* --------------------------------------------------------------------------- */

struct notEqual
{
    const IndexType mOutOfRange;

    notEqual( const IndexType val ) : mOutOfRange( val )
    {
    }

    __host__ __device__
    bool operator()( const IndexType x )
    {
        return x != mOutOfRange;
    }
};

__global__
static void get_col_pos_kernel( IndexType row[], IndexType pos[], const IndexType j,
                                const IndexType csrIA[], const IndexType numRows,
                                const IndexType csrJA[], const IndexType numValues )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        row[i] = numRows;     // out of range value indicates not found
        pos[i] = numValues;   // out of range value indicates not found

        for ( IndexType k = csrIA[i]; k < csrIA[i + 1]; ++k )
        {
            if ( csrJA[k] == j )
            {
                row[i] = i;
                pos[i] = k;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType CUDACSRUtils::getColumnPositions( 
    IndexType row[], 
    IndexType pos[], 
    const IndexType j,
    const IndexType csrIA[], 
    const IndexType numRows,
    const IndexType csrJA[], 
    const IndexType numValues )
{
    SCAI_REGION( "CUDA.CSRUtils.getColumnPositions" )

    SCAI_LOG_INFO( logger, "getColumnPositions: j = " << j << ", #rows = " << numRows << ", #nnz = " << numValues )

    SCAI_CHECK_CUDA_ACCESS

    // compute 'full' row, pos arrays

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    get_col_pos_kernel <<< dimGrid, dimBlock>>>( row, pos, j, csrIA, numRows, csrJA, numValues );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "get_row_kernel" )

    thrust::device_ptr<IndexType> d_pos( pos );
    thrust::device_ptr<IndexType> d_row( row );

    IndexType cnt1 = thrust::copy_if( d_pos,
                                      d_pos + numRows,
                                      d_pos,
                                      notEqual( numValues ) ) - d_pos;

    IndexType cnt2 = thrust::copy_if( d_row,
                                      d_row + numRows,
                                      d_row,
                                      notEqual( numRows ) ) - d_row;

    SCAI_ASSERT_EQ_ERROR( cnt1, cnt2, "serious size mismatch of row/pos arrays" )

    return cnt1;
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

// trivial kernel to check diagonal property
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

    bool found = false;

    for ( IndexType jj = ia[i]; jj < ia[i+1]; ++jj )
    {
        if ( ja[jj] == i )
        {
            found = true;
            break;
         }
    }
 
    if ( !found )
    {
        *hasProperty = false;
    }
}

bool CUDACSRUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrIA[], const IndexType csrJA[], const bool )
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

    bool* deviceFlag;
    bool hostFlag;     // will contain the flag copied from device

    SCAI_CUDA_DRV_CALL( cuMemAlloc( ( CUdeviceptr* ) &deviceFlag, sizeof( bool ) ),
                       "allocate bool flag on the device for the result of hasDiagonalProperty_kernel" )

    SCAI_CUDA_DRV_CALL( cuMemsetD8( ( CUdeviceptr ) deviceFlag, ( unsigned char ) 1, sizeof( bool ) ), "memset bool hasSortedRows = true" )

    hasDiagonalProperty_kernel <<< dimGrid, dimBlock>>>( numDiagonals, csrIA, csrJA, deviceFlag );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "hasDiagonalProperty::kernel failed, most likely arrays csrIA and/or csrJA are invalid" )

    SCAI_CUDA_DRV_CALL( cuMemcpyDtoH( &hostFlag, ( CUdeviceptr ) deviceFlag, sizeof( bool ) ),
                       "copy the result of hasDiagonalProperty_kernel to host" )

    SCAI_CUDA_DRV_CALL( cuMemFree( ( CUdeviceptr ) deviceFlag ), "cuMemFree( " << deviceFlag << " ) failed" )

    return hostFlag;
}

/* --------------------------------------------------------------------------- */

__inline__ __device__ bool isSorted( const IndexType* data, const IndexType n )
{
    bool sortedFlag = true;

    for ( IndexType i = 1; i < n; ++i )
    {
        if ( data[i] < data[i - 1] )
        {
            sortedFlag = false;
            break;
        }
    }

    return sortedFlag;
}

__global__ void hasSortedRows_kernel(
    bool* hasSortedRows,
    const IndexType* csrIA,
    const IndexType* csrJA,
    const IndexType numRows )
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= numRows )
    {
        return;
    }

    if ( ! *hasSortedRows )
    {
        return;
    }

    IndexType start = csrIA[i];

    const IndexType end = csrIA[i + 1];

    if ( ! isSorted( &csrJA[start], end - start ) )
    {
        *hasSortedRows = false;
    }
}

bool CUDACSRUtils::hasSortedRows(
    const IndexType csrIA[],
    const IndexType csrJA[],
    const IndexType numRows,
    const IndexType,
    const IndexType )
{
    SCAI_REGION( "CUDA.CSRUtils.hasSortedRows" )

    if ( numRows == 0 )
    {
        return true;
    }

    SCAI_LOG_INFO( logger, "hasSortedRows, numRows = " << numRows )

    SCAI_CHECK_CUDA_ACCESS

    // make grid

    const int blockSize = CUDASettings::getBlockSize(); 
    dim3 dimGrid( ( numRows - 1 ) / blockSize + 1, 1, 1 );// = makeGrid( numDiagonals, blockSize );
    dim3 dimBlock( blockSize, 1, 1 );

    bool* deviceFlag;
    bool hostFlag;     // will contain the flag copied from device

    SCAI_CUDA_DRV_CALL( cuMemAlloc( ( CUdeviceptr* ) &deviceFlag, sizeof( bool ) ),
                       "allocate global bool variable on the device for the result of hasSortedRows_kernel" )

    SCAI_CUDA_DRV_CALL( cuMemsetD8( ( CUdeviceptr ) deviceFlag, ( unsigned char ) 1, sizeof( bool ) ), "memset bool hasSortedRows = true" )

    hasSortedRows_kernel <<< dimGrid, dimBlock>>>( deviceFlag, csrIA, csrJA, numRows );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "hasSortedKernel failed: most likely arrays csrIA and/or csrJA are invalid" )

    SCAI_CUDA_DRV_CALL( cuMemcpyDtoH( &hostFlag, ( CUdeviceptr ) deviceFlag, sizeof( bool ) ),
                       "copy the result of hasSortedRows_kernel to host" )

    SCAI_CUDA_DRV_CALL( cuMemFree( ( CUdeviceptr ) deviceFlag ), "cuMemFree( " << deviceFlag << " ) failed" )

    return hostFlag;
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
    CUDACOOUtils::offsets2ia( cscJA, numValues, csrIA, numRows );
    // switch cooIA and cooJA, copy values and resort
    CUDASparseUtils::set( cooIA, csrJA, numValues, common::BinaryOp::COPY );
    CUDASparseUtils::set( cscValues, csrValues, numValues, common::BinaryOp::COPY );
    thrust::device_ptr<IndexType> ja_d( cooIA );
    thrust::device_ptr<ValueType> values_d( cscValues );
    thrust::device_ptr<IndexType> ia_d( cscJA );
    // sort by column indexes in ascending order
    // zip_iterator used to resort cscValues and cscJA in one step
    thrust::stable_sort_by_key( ja_d, ja_d + numValues,
                                thrust::make_zip_iterator( thrust::make_tuple( values_d, ia_d ) ) );
    // cscJA is now sorted, can become an offset array
    CUDACOOUtils::ia2offsets( cscIA, numColumns, cooIA, numValues );
    SCAI_CUDA_RT_CALL( cudaFree( cooIA ), "free tmp cooIA" )
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

template<typename ValueType>
__global__
void normal_gevm_kernel(
    ValueType result[],
    const ValueType x[],
    const ValueType alpha,
    const ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    IndexType numRows )
{
    // Note: atomicAdd dominates performance
    // result += alpha * x_d * A
    // result[j] += alpha * x_d[i] * A[i,j] for all (i,j) non-zero entries

    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd   = csrIA[i + 1];
        const ValueType xi       = x[i];

        for ( IndexType k = rowStart; k < rowEnd; ++k )
        {
            IndexType j = csrJA[k];
            ValueType v = alpha * csrValues[k] * xi;
            common::CUDAUtils::atomicAdd( &result[j], v );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void sparse_gevm_kernel(
    ValueType result[],
    const ValueType x[],
    const ValueType alpha,
    const ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const IndexType rowIndexes[],
    IndexType numRows )
{
    // Note: atomicAdd dominates performance
    // result += alpha * x_d * A
    // result[j] += alpha * x_d[i] * A[i,j] for all (i,j) non-zero entries

    const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( ii < numRows )
    {
        const IndexType i        = rowIndexes[ii];
        const IndexType rowStart = csrIA[i];
        const IndexType rowEnd   = csrIA[i + 1];
        const ValueType xi       = x[i];

        for ( IndexType k = rowStart; k < rowEnd; ++k )
        {
            IndexType j = csrJA[k];
            ValueType v = alpha * csrValues[k] * xi;
            common::CUDAUtils::atomicAdd( &result[j], v );
        }
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

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setRows                                                           */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void setRowsKernel(
    ValueType* values,
    const IndexType* ia,
    const IndexType numRows,
    const ValueType* diagonal,
    const common::BinaryOp op )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType tmp = diagonal[i];

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            values[jj] = common::applyBinary( values[jj], op, tmp );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDACSRUtils::setRows(
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const ValueType values[],
    const common::BinaryOp op  )
{
    SCAI_REGION( "CUDA.CSRUtils.setRows" )
    SCAI_LOG_INFO( logger, "setRows<" << TypeTraits<ValueType>::id() << ">"
                   << ", numrows= " << numRows )
    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    setRowsKernel <<< dimGrid, dimBlock>>>( csrValues, csrIA, numRows, values, op );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "CSRUtils:setRowsKernel FAILED" )
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
    const IndexType numColumns,
    const IndexType SCAI_UNUSED( nnz ),
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const common::MatrixOp op )
{
    SCAI_REGION( "CUDA.CSRUtils.normalGEMV" )

    if ( alpha == Constants::ZERO )
    {
        // result = beta * y 

        IndexType nTarget = common::isTranspose( op ) ? numColumns : numRows;

        CUDAUtils::binaryOpScalar( result, y, beta, nTarget, common::BinaryOp::MULT, false );

        return;
    }

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

    if ( common::isTranspose( op ) )
    {
        useTexture = false;

        // set result = beta * y, not needed if beta == 1 and y == result

        CUDAUtils::binaryOpScalar( result, y, beta, numColumns, common::BinaryOp::MULT, false );

        SCAI_LOG_DEBUG( logger, "Launch normal_gevm_kernel<" << TypeTraits<ValueType>::id() << ">" );

        normal_gevm_kernel<ValueType> <<< dimGrid, dimBlock, 0, stream >>>
        ( result, x, alpha, csrValues, csrIA, csrJA, numRows );
    }
    else if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == Constants::ONE && beta == Constants::ONE )
        {
            // result = A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == Constants::ONE && beta == Constants::ZERO )
        {
            // result = A * x_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == Constants::ONE )
        {
            // result = A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, beta, csrValues, csrIA, csrJA, numRows );
        }
        else if ( beta == Constants::ONE )
        {
            // result = alpha * A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows );
        }
        else if ( beta == Constants::ZERO )
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
        if ( alpha == Constants::ONE && beta == Constants::ONE )
        {
            // result = A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == Constants::ONE && beta == Constants::ZERO )
        {
            // result = A * x_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, csrValues, csrIA, csrJA, numRows );
        }
        else if ( alpha == Constants::ONE )
        {
            // result = A * x_d + beta * y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, beta, csrValues, csrIA, csrJA, numRows );
        }
        else if ( beta == Constants::ONE )
        {
            // result = alpha * A * x_d + y_d
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, y, alpha, csrValues, csrIA, csrJA, numRows );
        }
        else if ( beta == Constants::ZERO )
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
            syncToken->pushRoutine( std::bind( unbind, x ) );
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
    const ValueType csrValues[],
    const common::MatrixOp op )
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

    // Note: use of texture for x is only helpful for normal GEMV

    bool useTexture = CUDASettings::useTexture() && ( !common::isTranspose( op ) );

    if ( common::isTranspose( op ) )
    {
        sparse_gevm_kernel<ValueType> <<< dimGrid, dimBlock, 0, stream >>>
        ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numNonZeroRows );
    }
    else if ( useTexture )
    {
        vectorBindTexture( x );

        sparse_gemv_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
        ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numNonZeroRows );
    }
    else
    {
        sparse_gemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
        ( result, x, alpha, csrValues, csrIA, csrJA, rowIndexes, numNonZeroRows );
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
            syncToken->pushRoutine( std::bind( unbind, x ) );
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
        ValueType diag = 0;

        for ( IndexType jj = rowStart; jj < rowEnd; ++jj )
        {
            const IndexType j = csrJA[jj];

            if ( j == i )
            {
                diag = csrValues[jj];
            }
            else
            {
                temp -= csrValues[jj] * fetchVectorX<ValueType, useTexture>( oldSolution, j );
            }
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
    SCAI_REGION( "CUDA.CSR.jacobi" )

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
void CUDACSRUtils::jacobiHalo(
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
    SCAI_REGION( "CUDA.CSR.jacobiHalo" )

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
    }
    else
    {
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHalo_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }

    if ( useTexture )
    {
        csr_jacobiHalo_kernel <ValueType, true> <<< dimGrid, dimBlock>>>( solution, localDiagValues, haloIA,
                haloJA, haloValues, haloRowIndexes,
                numNonEmptyRows, oldSolution, omega );
    }
    else
    {
        csr_jacobiHalo_kernel<ValueType, false> <<< dimGrid, dimBlock>>>( solution, localDiagValues, haloIA,
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
    IndexType numWarpsGlobal = ( blockDim.x * gridDim.x ) / warpSize;
    IndexType rowIt = globalWarpId;

    for ( ; CUDA_ANY( rowIt < numRows ); rowIt += numWarpsGlobal )
    {
        if ( rowIt < numRows )
        {
            IndexType aColIt = aIa[rowIt] + laneId;
            IndexType aColEnd = aIa[rowIt + 1];
            IndexType bColIt = bIa[rowIt] + laneId;
            IndexType bColEnd = bIa[rowIt + 1];

            if ( laneId == 0 )
            {
                cIa[rowIt] = bColEnd - bColIt;
            }

            for ( IndexType aColItOffset = 0; CUDA_ANY( aColIt < aColEnd ); aColIt += warpSize, aColItOffset += warpSize )
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

                    for ( IndexType bColItOffset = 0; !sFound[localWarpId] &&  CUDA_ANY( ( bColIt + bColItOffset ) < bColEnd );
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
    const IndexType aIa[],
    const IndexType aJa[],
    const IndexType bIa[],
    const IndexType bJa[] )
{
    SCAI_REGION( "CUDA.CSRUtils.matrixAddSizes" )
    SCAI_LOG_INFO(
        logger,
        "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" )
    SCAI_CHECK_CUDA_ACCESS
// Reset cIa
    thrust::device_ptr<IndexType> cIaPtr( cIa );
    thrust::fill( cIaPtr, cIaPtr + numRows, 0 );
// TODO: Check if diagonal property needs special attention
    matrixAddSizesKernel<NUM_WARPS> <<< NUM_BLOCKS, NUM_THREADS>>>( cIa, numRows, numColumns, aIa, aJa, bIa, bJa );
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

    if ( laneId == 0 && chunkCount > 0 )
    {
        // This loop should also work for unsigned index type

        for ( IndexType i = *sReservedChunks; --i > ( *sReservedChunks - chunkCount ); )
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
    bool* hashError )
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

                multHlp_initializeChunks( sHashTable,
                                          chunkPtr,
                                          NUM_ELEMENTS_PER_CHUNK,
                                          sChunkList,
                                          sReservedChunks );

                for ( IndexType offset = 0; CUDA_ANY( aColIt < aColEnd ); aColIt += warpSize, offset += warpSize )
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

                        for ( ; CUDA_ANY( bColIt < bColEnd ); bColIt += warpSize )
                        {
                            colB = bColIt < bColEnd ? bJA[bColIt] : cudaNIndex;

                            if ( colB != cudaNIndex )
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

                if ( CUDA_ANY( localSystemError ) )
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
    const IndexType k,
    const IndexType aIa[],
    const IndexType aJa[],
    const IndexType bIa[],
    const IndexType bJa[] )
{
    SCAI_REGION( "CUDA.CSR.matrixMultiplySizes" )
    SCAI_LOG_INFO( logger, "matrixMultiplySizes for " << numRows << " x " << numColumns << " matrix" )
    SCAI_CHECK_CUDA_ACCESS

    // Reset cIa
    thrust::device_ptr<IndexType> cIaPtr( cIa );
    thrust::fill( cIaPtr, cIaPtr + numRows, 0 );

    ContextPtr loc = Context::getContextPtr( ContextType::CUDA );
    MemoryPtr mem = loc->getMemoryPtr();

    bool hashErrorHost = false;
    bool* hashError = ( bool* ) mem->allocate( sizeof( bool ) );
    SCAI_CUDA_RT_CALL( cudaMemcpy( hashError, &hashErrorHost, sizeof( bool ), cudaMemcpyHostToDevice ), "memcpy of hashError" );

    size_t free;
    size_t total;
    cuMemGetInfo( &free, &total );
    SCAI_LOG_DEBUG( logger, "free = " << free << ", total = " << total )

    IndexType nnz_a;
    IndexType nnz_b;
    SCAI_CUDA_RT_CALL( cudaMemcpy( &nnz_a, &aIa[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost ), "memcpy of nnz_a" );
    SCAI_CUDA_RT_CALL( cudaMemcpy( &nnz_b, &bIa[k], sizeof( IndexType ), cudaMemcpyDeviceToHost ), "memcpy of nnz_b" );

    IndexType avgDensity = ( nnz_a / numRows + nnz_b / numColumns ) / 2;
    IndexType numChunks;
    SCAI_ASSERT_GT_ERROR ( free, static_cast<IndexType>( 100 * 1024 * 1024 ), "insufficient free memory" );
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

    SCAI_LOG_DEBUG( logger, "numChunks = " << numChunks << ", max = " << maxNumChunks << ", per warp = " << chunksPerWarp )

    size_t hashTableAllocatedBytes = static_cast<size_t>( numChunks ) * NUM_ELEMENTS_PER_CHUNK * sizeof( IndexType );

    SCAI_LOG_DEBUG( logger, "hashTableAllcoatedBytes= " << hashTableAllocatedBytes )

    IndexType* hashTable = reinterpret_cast<IndexType*>( mem->allocate( hashTableAllocatedBytes ) );

    // chunkList table needs one integer per chunk plus 1 start pointer

    size_t chunkListAllocatedBytes = static_cast<size_t>( numChunks ) * sizeof( IndexType ) + sizeof( IndexType );

    IndexType* chunkList = reinterpret_cast<IndexType*>( mem->allocate( chunkListAllocatedBytes ) );

    thrust::device_ptr<IndexType> chunkListPtr( chunkList );
    thrust::transform( thrust::make_counting_iterator( IndexType( 0 ) ),
                       thrust::make_counting_iterator( numChunks + 1 ),
                       chunkListPtr,
                       multHlp_chunkFill( numChunks + 1 ) );

    matrixMultiplySizesKernel <<< NUM_BLOCKS, NUM_THREADS>>>(
        aIa,
        aJa,
        bIa,
        bJa,
        cIa,
        numRows,
        numColumns,
        hashTable,
        chunkList,
        numChunks,
        hashError );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "snyc after matrixMultiplySizesKernel" );

    SCAI_CUDA_RT_CALL( cudaMemcpy( &hashErrorHost, hashError, sizeof( bool ), cudaMemcpyDeviceToHost ), "memcpy hashError" );

    if ( hashErrorHost )
    {
        COMMON_THROWEXCEPTION( "Multiplication for Sizes failed!" );
    }

    // Free hashTable and hashError

    mem->free( ( void* ) hashError, sizeof( bool ) );
    mem->free( ( void* ) hashTable, hashTableAllocatedBytes );
    mem->free( ( void* ) chunkList, chunkListAllocatedBytes );

    // Convert sizes array to offset array
    thrust::exclusive_scan( cIaPtr, cIaPtr + numRows + 1, cIaPtr );
    IndexType numValues;
    cudaMemcpy( &numValues, &cIa[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost );

    SCAI_CHECK_CUDA_ERROR

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
    IndexType numWarpsGlobal = ( blockDim.x * gridDim.x ) / warpSize;
    IndexType rowIt = globalWarpId;

    for ( ; CUDA_ANY( rowIt < numRows ); rowIt += numWarpsGlobal )
    {
        if ( rowIt < numRows )
        {
            IndexType aColIt = aIA[rowIt] + laneId;
            IndexType aColEnd = aIA[rowIt + 1];
            IndexType bColIt = bIA[rowIt] + laneId;
            IndexType bColEnd = bIA[rowIt + 1];
            IndexType cColIt = cIA[rowIt] + laneId;

            // Copy values of b to C

            for ( IndexType bColOffset = 0; CUDA_ANY( ( bColIt + bColOffset ) < bColEnd ); bColOffset += warpSize )
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
            for ( IndexType aColItOffset = 0; CUDA_ANY( aColIt < aColEnd ); aColIt += warpSize, aColItOffset += warpSize )
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

                    for ( IndexType bColItOffset = 0; ( sFoundJa[localWarpId] == cudaNIndex ) && CUDA_ANY( ( bColIt + bColItOffset ) < bColEnd );
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
                            cJA[cColIt + cColOffset] = sColA[localWarpId];
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
            alpha, aIA, aJA, aValues, beta, bIA, bJA, bValues );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "sync after matrixAdd kernel" )
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
                                    IndexType numReservedChunks )
{
    // TODO: rename sColA => destinationOffset!

    *sColA = 0;
    IndexType rowOffset = cIA[aRowIt];
    IndexType one = 1;
    IndexType hashCol;
    ValueType hashVal;

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
    bool* hashError )
{
    __shared__ IndexType sHashTableIndexes[NUM_ELEMENTS_IN_SHARED];
    __shared__ ValueType sHashTableValues[NUM_ELEMENTS_IN_SHARED];
    __shared__ volatile IndexType sReservedChunks;
    __shared__ volatile IndexType sChunkList[NUM_CHUNKS_PER_WARP];
    __shared__ volatile IndexType sColA;
    __shared__ volatile ValueType sValA;
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
            IndexType optimalChunkCount = multHlp_calcOptChunkCount ( aRowIt, cIA, NUM_ELEMENTS_PER_CHUNK );

            // reserve Chunks
            if ( !multHlp_reserveChunks( chunkList, sChunkList, &sReservedChunks, optimalChunkCount ) )
            {
                // ABORT KERNEL HERE;
                localSystemError = true;
            }

            if ( CUDA_ANY( localSystemError ) )
            {
                *hashError = true;
                return;
            }

            do
            {
                sInsertMiss = false;
                IndexType aColIt = aIA[aRowIt] + laneId;
                IndexType aColEnd = aIA[aRowIt + 1];

                multHlp_initializeChunks( sHashTableIndexes,
                                          indexChunks,
                                          NUM_ELEMENTS_PER_CHUNK,
                                          sChunkList,
                                          sReservedChunks );

                for ( IndexType offset = 0; CUDA_ANY( aColIt < aColEnd ); aColIt += warpSize, offset += warpSize )
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

                        for ( ; CUDA_ANY( bColIt < bColEnd ); bColIt += warpSize )
                        {
                            colB = bColIt < bColEnd ? bJA[bColIt] : cudaNIndex;
                            ValueType valB = bColIt < bColEnd ? bValues[bColIt] : static_cast<ValueType>( 0 );

                            if ( colB != cudaNIndex  )
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
                                            sReservedChunks );
                }
                else
                {
                    if ( !multHlp_reserveChunks( chunkList, sChunkList, &sReservedChunks, multHlp_growth( sReservedChunks ) ) )
                    {
                        // ABORT KERNEL HERE;
                        localSystemError = true;
                    }

                    if ( CUDA_ANY( localSystemError ) )
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
    const IndexType k,
    const ValueType alpha,
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

    ContextPtr loc = Context::getContextPtr( ContextType::CUDA );
    MemoryPtr mem = loc->getMemoryPtr();

    bool hashErrorHost = false;
    bool* hashError = ( bool* ) mem->allocate( sizeof( bool ) );
    SCAI_CUDA_RT_CALL( cudaMemcpy( hashError, &hashErrorHost, sizeof( bool ), cudaMemcpyHostToDevice ), "memcpy of hashError" );

    size_t free;
    size_t total;
    cuMemGetInfo( &free, &total );

    IndexType nnz_a;
    IndexType nnz_b;
    SCAI_CUDA_RT_CALL( cudaMemcpy( &nnz_a, &aIa[numRows], sizeof( IndexType ), cudaMemcpyDeviceToHost ), "memcpy of nnz_a" );
    SCAI_CUDA_RT_CALL( cudaMemcpy( &nnz_b, &bIa[k], sizeof( IndexType ), cudaMemcpyDeviceToHost ), "memcpy of nnz_b" );

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

    SCAI_LOG_DEBUG( logger, "numChunks = " << numChunks << ", max = " << maxNumChunks << ", per warp = " << chunksPerWarp )

    unsigned int hashTableAllocatedBytes = numChunks * NUM_ELEMENTS_PER_CHUNK *
                                           ( sizeof( IndexType ) + sizeof( ValueType ) );

    SCAI_LOG_DEBUG( logger, "hashTableAllcoatedBytes= " << hashTableAllocatedBytes )

    void* chunks = ( void* ) mem->allocate( hashTableAllocatedBytes );
    IndexType* indexChunks = ( IndexType* ) chunks;
    ValueType* valueChunks = ( ValueType* ) ( indexChunks + numChunks * NUM_ELEMENTS_PER_CHUNK );

    // chunkList table needs one integers per chunk plus 1 start pointer
    unsigned int chunkListAllocatedBytes = numChunks * sizeof( IndexType ) + sizeof( IndexType );

    IndexType* chunkList = ( IndexType* ) mem->allocate( chunkListAllocatedBytes );

    thrust::device_ptr<IndexType> chunkListPtr( chunkList );
    thrust::transform( thrust::make_counting_iterator( IndexType( 0 ) ),
                       thrust::make_counting_iterator( numChunks + 1 ),
                       chunkListPtr,
                       multHlp_chunkFill( numChunks + 1 ) );

    matrixMultiplyKernel <<< NUM_BLOCKS, NUM_THREADS>>>(
        aIa,
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
        hashError );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "sync after matrixMultiply kernel" )

    SCAI_CUDA_RT_CALL( cudaMemcpy( &hashErrorHost, hashError, sizeof( bool ), cudaMemcpyDeviceToHost ), "memcpy hashError" );

    if ( hashErrorHost )
    {
        COMMON_THROWEXCEPTION( "Multiplication failed!" );
    }

    // Free hashTable and hashError
    mem->free( ( void* ) hashError, sizeof( bool ) );
    mem->free( ( void* ) chunks, hashTableAllocatedBytes );
    mem->free( ( void* ) chunkList, chunkListAllocatedBytes );

    SCAI_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */

namespace gpu
{
template <class T>
__device__ void swap ( T& a, T& b )
{
    T c( a );
    a = b;
    b = c;
}
}

template<typename ValueType>
__global__
void sortRowKernel(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )

    {
        // use serial bubble sort as sort algorithm for one row

        IndexType start = csrIA[i];
        IndexType end   = csrIA[i + 1] - 1;

        bool sorted = false;

        while ( !sorted )
        {
            sorted = true; // will be reset if any wrong order appears

            for ( IndexType jj = start; jj < end; ++jj )
            {
                bool swapIt = csrJA[jj] > csrJA[jj + 1];

                if ( swapIt )
                {
                    sorted = false;
                    gpu::swap( csrJA[jj], csrJA[jj + 1] );
                    gpu::swap( csrValues[jj], csrValues[jj + 1] );
                }
            }

            --end;
        }
    }
}

template<typename ValueType>
void CUDACSRUtils::sortRows(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType,
    const IndexType )
{
    SCAI_REGION( "CUDA.CSR.sortRow" )

    SCAI_LOG_INFO( logger, "sort elements in each of " << numRows << " rows" )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    sortRowKernel <<< dimGrid, dimBlock>>>( csrJA, csrValues, csrIA, numRows );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "sortRows" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void shiftDiagKernel(
    IndexType count[],
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType diagonalIndexes[],
    const IndexType numDiagonals )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i >= numDiagonals )
    {
        return;
    }

    IndexType diagonalIndex = i;
 
    if ( diagonalIndexes != NULL )
    {
        diagonalIndex = diagonalIndexes[i];
    }

    IndexType start = csrIA[i];
    IndexType end   = csrIA[i+1] - 1;

    if ( end < start )
    {
        count[i] = 0;  // empty row, no diagonal element
        return;
    }

    if ( csrJA[start] == diagonalIndex )
    {
        count[i] = 1;  // diagonal element is already first
        return;
    }

    bool found = false;
    count[i]   = 0;      // count this row only if diagonal element found

    ValueType diagonalValue;  // will be set when found

    // traverse reverse

    while ( end > start )
    {
        // check if it is the diagonal element, save the diagonal value

        if ( not found && csrJA[end] == diagonalIndex )
        {
            found = true;
            diagonalValue = csrValues[end];
        }

        // move up elements to fill the gap of diagonal element
        if ( found )
        {
            csrJA[end] = csrJA[end - 1];
            csrValues[end] = csrValues[end - 1];
        }

        end--;
    }

    if ( found )
    {
        // now set the first row element as the diagonal element
        csrValues[start] = diagonalValue;
        csrJA[start] = diagonalIndex;
        count[i] = 1;
    }
}

template<typename ValueType>
IndexType CUDACSRUtils::shiftDiagonal(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType numDiagonals,
    const IndexType csrIA[],
    const IndexType diagonalIndexes[] )
{
    SCAI_REGION( "CUDA.CSR.shiftDiag" )

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_vector<IndexType> count( numDiagonals );

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numDiagonals, dimBlock.x );

    shiftDiagKernel <<< dimGrid, dimBlock>>>( count.data().get(), csrJA, csrValues, csrIA, diagonalIndexes, numDiagonals );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "shiftDiagonal" )

    IndexType numFirstDiagonals = thrust::reduce( count.begin(), count.end(), 0, thrust::plus<IndexType>() );

    SCAI_LOG_INFO( logger, "shiftDiagonal for CSR data, " << numFirstDiagonals << " of " << numDiagonals << " are now first entry" )

    return numFirstDiagonals;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void countNonZerosKernel(
    IndexType sizes[],
    const IndexType ia[],
    const IndexType ja[],
    const ValueType values[],
    const IndexType numRows,
    const RealType<ValueType> eps )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType cnt = 0;

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            bool nonZero = common::Math::abs( values[jj] ) > common::Math::abs( eps );

            if ( nonZero )
            {
                ++cnt;
            }
        }

        sizes[i] = cnt;
    }
}

template<typename ValueType>
void CUDACSRUtils::countNonZeros(
    IndexType sizes[],
    const IndexType ia[],
    const IndexType ja[],
    const ValueType values[],
    const IndexType numRows,
    const RealType<ValueType> eps )
{
    SCAI_REGION( "CUDA.CSRUtils.countNonZeros" )

    SCAI_LOG_INFO( logger, "countNonZeros of CSR<" << TypeTraits<ValueType>::id() << ">( " << numRows
                   << "), eps = " << eps )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    countNonZerosKernel <<< dimGrid, dimBlock>>>( sizes, ia, ja, values, numRows, eps );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void compressKernel(
    IndexType newJA[],
    ValueType newValues[],
    const IndexType newIA[],
    const IndexType ia[],
    const IndexType ja[],
    const ValueType values[],
    const IndexType numRows,
    const RealType<ValueType> eps )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType offs = newIA[i];

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            if ( common::Math::abs( values[jj] ) <= eps )
            {
                continue;  // skip this zero value
            }

            newJA[ offs ]     = ja[jj];
            newValues[ offs ] = values[jj];
            ++offs;
        }
    }
}

template<typename ValueType>
void CUDACSRUtils::compress(
    IndexType newJA[],
    ValueType newValues[],
    const IndexType newIA[],
    const IndexType ia[],
    const IndexType ja[],
    const ValueType values[],
    const IndexType numRows,
    const RealType<ValueType> eps )
{
    SCAI_REGION( "CUDA.CSR.compress" )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    compressKernel <<< dimGrid, dimBlock>>>( newJA, newValues, newIA, ia, ja, values, numRows, eps );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "compress" )
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDACSRUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "set CSR routines for CUDA in Interface" )
    KernelRegistry::set<CSRKernelTrait::sizes2offsets>( sizes2offsets, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::offsets2sizes>( offsets2sizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::getColumnPositions>( getColumnPositions, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::hasSortedRows>( hasSortedRows, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAddSizes>( matrixAddSizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, ctx, flag );
}

template<typename ValueType>
void CUDACSRUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register CSRUtils CUDA-routines for CUDA at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<ValueType> >( convertCSR2CSC, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sortRows<ValueType> >( sortRows, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::shiftDiagonal<ValueType> >( shiftDiagonal, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::compress<ValueType> >( compress, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::countNonZeros<ValueType> >( countNonZeros, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAdd<ValueType> >( matrixAdd, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiply<ValueType> >( matrixMultiply, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::setRows<ValueType> >( setRows, ctx, flag );
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
}

CUDACSRUtils::~CUDACSRUtils()
{
    SCAI_LOG_INFO( logger, "unregister CSRUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDACSRUtils CUDACSRUtils::guard;    // guard variable for registration

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

unsigned int CUDACSRUtils::lastHashTableSize = 1024;

} /* end namespace sparsekernel */

} /* end namespace scai */
