/**
 * @file CUDADIAUtils.cu
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
 * @brief Implementation of DIA utilities with CUDA
 * @author Bea Hornef, Thomas Brandes
 * @date 04.07.2012
 */

// hpp
#include <scai/sparsekernel/cuda/CUDADIAUtils.hpp>

// local library
#include <scai/sparsekernel/DIAKernelTrait.hpp>


// internal scai library
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>

// thrust
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

using namespace scai::tasking;

namespace scai
{

using common::TypeTraits;
using common::CUDASettings;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( CUDADIAUtils::logger, "CUDA.DIAUtils" )

/* --------------------------------------------------------------------------- */

template<bool useTexture, bool useSharedMemory>
__inline__ __device__
IndexType fetchOffset( const IndexType* const offset_d, IndexType[], const IndexType i )
{
    return offset_d[i];
}

template<>
__inline__ __device__
IndexType fetchOffset<true, false>( const IndexType* const offset_d, IndexType[], const IndexType i )
{
    return fetchVectorX<IndexType, true>( offset_d, i );
}

template<>
__inline__ __device__
IndexType fetchOffset<true, true>( const IndexType* const, IndexType offset_sm[], const IndexType i )
{
    return offset_sm[i];
}

template<>
__inline__ __device__
IndexType fetchOffset<false, true>( const IndexType* const, IndexType offset_sm[], const IndexType i )
{
    return offset_sm[i];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void normal_gemv_kernel(
    ValueType* result,
    const ValueType* x,
    const ValueType* y,
    const ValueType alpha,
    const ValueType beta,
    const ValueType* diagonalValues,
    const IndexType* offsets_d,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals )
{
    extern __shared__ IndexType offsets_sm[];

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType summand = beta * y[i];
        ValueType temp = 0.0;

        for ( IndexType idiag = 0; idiag < numDiagonals; idiag++ )
        {
            IndexType j = i + fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, idiag );

            if ( common::Utils::validIndex( j, numColumns ) )
            {
                ValueType val = diagonalValues[ numRows * idiag + i ];
                temp += val * fetchVectorX<ValueType, useTexture>( x, j );
            }
        }

        result[i] = alpha * temp + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void normal_gemv_kernel_alpha_one_beta_one(
    ValueType* result,
    const ValueType* x,
    const ValueType* y,
    const ValueType* diagonalValues,
    const IndexType* offsets_d,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals )
{
    extern __shared__ IndexType offsets_sm[];

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType summand = y[i];
        ValueType temp = 0.0;

        for ( IndexType idiag = 0; idiag < numDiagonals; idiag++ )
        {
            IndexType j = i + fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, idiag );

            if ( common::Utils::validIndex( j, numColumns ) )
            {
                ValueType val = diagonalValues[ numRows * idiag + i ];
                temp += val * fetchVectorX<ValueType, useTexture>( x, j );
            }
        }

        result[i] = temp + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void normal_gemv_kernel_alpha_one_beta_zero(
    ValueType* result,
    const ValueType* x,
    const ValueType* y,
    const ValueType* diagonalValues,
    const IndexType* offsets_d,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals )
{
    extern __shared__ IndexType offsets_sm[];

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType temp = 0.0;

        for ( IndexType idiag = 0; idiag < numDiagonals; idiag++ )
        {
            IndexType j = i + fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, idiag );

            if ( common::Utils::validIndex( j, numColumns ) )
            {
                ValueType val = diagonalValues[ numRows * idiag + i ];
                temp += val * fetchVectorX<ValueType, useTexture>( x, j );
            }
        }

        result[i] = temp;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__
void assign_kernel(
    ValueType* result,
    const ValueType* y,
    const IndexType numRows )
{
    IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        result[i] = y[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void normal_gemv_kernel_alpha_one(
    ValueType* result,
    const ValueType* x,
    const ValueType* y,
    const ValueType beta,
    const ValueType* diagonalValues,
    const IndexType* offsets_d,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals )
{
    extern __shared__ IndexType offsets_sm[];

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType summand = 0.0;

        if ( beta != 0.0 )
        {
            summand = beta * y[i];
        }

        ValueType temp = 0.0;

        for ( IndexType idiag = 0; idiag < numDiagonals; idiag++ )
        {
            IndexType j = i + fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, idiag );

            if ( common::Utils::validIndex( j, numColumns ) )
            {
                ValueType val = diagonalValues[ numRows * idiag + i ];
                temp += val * fetchVectorX<ValueType, useTexture>( x, j );
            }
        }

        result[i] = temp + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void normal_gemv_kernel_alpha_zero(
    ValueType* result,
    const ValueType* x,
    const ValueType* y,
    const ValueType beta,
    const ValueType* diagonalValues,
    const IndexType* offsets_d,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals )
{
    extern __shared__ IndexType offsets_sm[];

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        result[i] = beta * y[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void normal_gemv_kernel_beta_one(
    ValueType* result,
    const ValueType* x,
    const ValueType* y,
    const ValueType alpha,
    const ValueType* diagonalValues,
    const IndexType* offsets_d,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals )
{
    extern __shared__ IndexType offsets_sm[];

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType summand = y[i];
        ValueType temp = 0.0;

        for ( IndexType idiag = 0; idiag < numDiagonals; idiag++ )
        {
            IndexType j = i + fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, idiag );

            if ( common::Utils::validIndex( j, numColumns ) )
            {
                ValueType val = diagonalValues[ numRows * idiag + i ];
                temp += val * fetchVectorX<ValueType, useTexture>( x, j );
            }
        }

        result[i] = alpha * temp + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void normal_gemv_kernel_beta_zero(
    ValueType* result,
    const ValueType* x,
    const ValueType* y,
    const ValueType alpha,
    const ValueType* diagonalValues,
    const IndexType* offsets_d,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals )
{
    extern __shared__ IndexType offsets_sm[];

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType temp = 0.0;

        for ( IndexType idiag = 0; idiag < numDiagonals; idiag++ )
        {
            IndexType j = i + fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, idiag );

            if ( common::Utils::validIndex( j, numColumns ) )
            {
                ValueType val = diagonalValues[ numRows * idiag + i ];
                temp += val * fetchVectorX<ValueType, useTexture>( x, j );
            }
        }

        result[i] = alpha * temp;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDADIAUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    SCAI_REGION( "CUDA.DIA.normalGEMV" )
    SCAI_LOG_INFO( logger, "normalGEMV<" << TypeTraits<ValueType>::id() << ">"
                   << " result[ " << numRows << "] = " << alpha
                   << " * A( #diags = " << numDiagonals << " ) * x + " << beta << " * y " )
    const IndexType blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    const bool useSharedMem = CUDASettings::useSharedMem();

    const bool useTexture = CUDASettings::useTexture();

    SCAI_LOG_INFO( logger, "Start normal_gemv_kernel<" << TypeTraits<ValueType>::id()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

    int sharedMemSize = 0;

    if ( useSharedMem )
    {
        sharedMemSize = numDiagonals * sizeof( int );
    }

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( !useSharedMem )
        {
            vectorBindTexture( diaOffsets );
        }

        if ( useSharedMem )
        {
            if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ONE )
            {
                normal_gemv_kernel_alpha_one_beta_one<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_alpha_one_beta_zero<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ZERO && beta == scai::common::constants::ONE )
            {
                assign_kernel<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, y, numRows );
            }
            else if ( alpha == scai::common::constants::ONE )
            {
                normal_gemv_kernel_alpha_one<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_alpha_zero<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( beta == scai::common::constants::ONE )
            {
                normal_gemv_kernel_beta_one<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( beta == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_beta_zero<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else
            {
                normal_gemv_kernel<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
        }
        else
        {
            if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ONE )
            {
                normal_gemv_kernel_alpha_one_beta_one<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_alpha_one_beta_zero<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ZERO && beta == scai::common::constants::ONE )
            {
                assign_kernel<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, y, numRows );
            }
            else if ( alpha == scai::common::constants::ONE )
            {
                normal_gemv_kernel_alpha_one<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_alpha_zero<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( beta == scai::common::constants::ONE )
            {
                normal_gemv_kernel_beta_one<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( beta == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_beta_zero<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else
            {
                normal_gemv_kernel<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
        }
    }
    else
    {
        if ( useSharedMem )
        {
            if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ONE )
            {
                normal_gemv_kernel_alpha_one_beta_one<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_alpha_one_beta_zero<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ZERO && beta == scai::common::constants::ONE )
            {
                assign_kernel<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, y, numRows );
            }
            else if ( alpha == scai::common::constants::ONE )
            {
                normal_gemv_kernel_alpha_one<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_alpha_zero<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( beta == scai::common::constants::ONE )
            {
                normal_gemv_kernel_beta_one<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( beta == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_beta_zero<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else
            {
                normal_gemv_kernel<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
        }
        else
        {
            if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ONE )
            {
                normal_gemv_kernel_alpha_one_beta_one<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_alpha_one_beta_zero<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ZERO && beta == scai::common::constants::ONE )
            {
                assign_kernel<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, y, numRows );
            }
            else if ( alpha == scai::common::constants::ONE )
            {
                normal_gemv_kernel_alpha_one<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( alpha == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_alpha_zero<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( beta == scai::common::constants::ONE )
            {
                normal_gemv_kernel_beta_one<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else if ( beta == scai::common::constants::ZERO )
            {
                normal_gemv_kernel_beta_zero<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else
            {
                normal_gemv_kernel<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                    result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
        }
    }

    if ( !syncToken )
    {
        // synchronize now, unbind used textures
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "normalGEMV for DIA" )

        if ( useTexture )
        {
            vectorUnbindTexture( x );

            if ( !useSharedMem )
            {
                vectorUnbindTexture( diaOffsets );
            }
        }
    }
    else
    {
        // synchronize by syncToken, delay unbind texture
        if ( useTexture )
        {
            void ( *unbindV ) ( const ValueType* ) = &vectorUnbindTexture;
            void ( *unbindI ) ( const IndexType* ) = &vectorUnbindTexture;
            syncToken->pushRoutine( common::bind( unbindV, x ) );

            if ( !useSharedMem )
            {
                syncToken->pushRoutine( common::bind( unbindI, diaOffsets ) );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void normal_gevm_kernel(
    ValueType* result,
    const ValueType* x,
    const ValueType* y,
    const ValueType alpha,
    const ValueType beta,
    const ValueType* diagonalValues,
    const IndexType* offsets_d,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals )
{
    extern __shared__ IndexType offsets_sm[];

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    IndexType k = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( k < numColumns )
    {
        ValueType summand = 0;

        if ( beta != 0 )
        {
            summand = beta * y[k];
        }

        ValueType temp = 0.0;

        for ( IndexType ii = 0; ii < numDiagonals; ii++ )
        {
            IndexType i = k - fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, ii );

            if ( common::Utils::validIndex( i, numRows ) )
            {
                temp += diagonalValues[ numRows * ii + i ] * fetchVectorX<ValueType, useTexture>( x, i );
            }
        }

        result[k] = alpha * temp + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDADIAUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    SCAI_REGION( "CUDA.DIA.normalGEVM" )
    SCAI_LOG_INFO( logger, "normalGEVM<" << TypeTraits<ValueType>::id() << ">"
                   << " result[ " << numRows << "] = " << alpha
                   << " * A( #diags = " << numDiagonals << " ) * x + " << beta << " * y " )
    const IndexType blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numColumns, dimBlock.x );
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    const bool useSharedMem = CUDASettings::useSharedMem();

    const bool useTexture = CUDASettings::useTexture();

    SCAI_LOG_INFO( logger, "Start normal_gevm_kernel<" << TypeTraits<ValueType>::id()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

    int sharedMemSize = 0;

    if ( useSharedMem )
    {
        sharedMemSize = numDiagonals * sizeof( int );
    }

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( !useSharedMem )
        {
            // @ToDo: be careful, some CUDA devices do not support multiple bind textures, e.g. GeForce 460
            vectorBindTexture( diaOffsets );
        }

        if ( useSharedMem )
        {
            normal_gevm_kernel<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
        }
        else
        {
            normal_gevm_kernel<ValueType, true, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
        }
    }
    else
    {
        if ( useSharedMem )
        {
            normal_gevm_kernel<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
        }
        else
        {
            normal_gevm_kernel<ValueType, false, false> <<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
        }
    }

    if ( !syncToken )
    {
        // synchronize now, unbind used textures
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "normalGEMV for DIA" )

        if ( useTexture )
        {
            vectorUnbindTexture( x );

            if ( !useSharedMem )
            {
                vectorUnbindTexture( diaOffsets );
            }
        }
    }
    else
    {
        // synchronize by syncToken, delay unbind texture
        if ( useTexture )
        {
            void ( *unbindV ) ( const ValueType* ) = &vectorUnbindTexture;
            void ( *unbindI ) ( const IndexType* ) = &vectorUnbindTexture;
            syncToken->pushRoutine( common::bind( unbindV, x ) );

            if ( !useSharedMem )
            {
                syncToken->pushRoutine( common::bind( unbindI, diaOffsets ) );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDADIAUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    SCAI_LOG_DEBUG( logger, "register DIAUtils CUDA-routines for CUDA at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    const common::context::ContextType ctx = common::context::CUDA;
    KernelRegistry::set<DIAKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::normalGEVM<ValueType> >( normalGEVM, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDADIAUtils::CUDADIAUtils()
{
    SCAI_LOG_INFO( logger, "register DIAUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDADIAUtils::~CUDADIAUtils()
{
    SCAI_LOG_INFO( logger, "unregister DIAUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDADIAUtils CUDADIAUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
