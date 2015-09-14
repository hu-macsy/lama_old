/**
 * @file CUDADIAUtils.cu
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Implementation of DIA utilities with CUDA
 * @author Bea Hornef, Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

#include <scai/common/Assert.hpp>
#include <scai/common/bind.hpp>

#include <scai/common/cuda/CUDAError.hpp>

#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>

#include <scai/lama/cuda/utils.cu.h>
#include <scai/lama/cuda/CUDADIAUtils.hpp>
#include <scai/lama/cuda/CUDASettings.hpp>

// thrust
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

#include <boost/preprocessor.hpp>

using namespace scai::hmemo;
using namespace scai::tasking;

using scai::common::getScalarType;

using scai::common::ComplexFloat;
using scai::common::ComplexDouble;

namespace scai
{

namespace lama
{

    SCAI_LOG_DEF_LOGGER( CUDADIAUtils::logger, "CUDA.DIAUtils" )

    /* --------------------------------------------------------------------------- */

#include <scai/lama/cuda/CUDATexVector.hpp>

    /* --------------------------------------------------------------------------- */

    template<bool useTexture, bool useSharedMemory>
    __inline__ __device__
    int fetchOffset( const int* const offset_d, int[], const int i )
    {
        return offset_d[i];
    }

    template<>
    __inline__ __device__
    int fetchOffset<true, false>( const int* const offset_d, int[], const int i )
    {
        return fetchVectorX<int, true>( offset_d, i );
    }

    template<>
    __inline__ __device__
    int fetchOffset<true, true>( const int* const, int offset_sm[], const int i )
    {
        return offset_sm[i];
    }

    template<>
    __inline__ __device__
    int fetchOffset<false, true>( const int* const, int offset_sm[], const int i )
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
            int k = threadIdx.x;
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

                if ( j >= 0 && j < numColumns )
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
            int k = threadIdx.x;
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

                if ( j >= 0 && j < numColumns )
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
            int k = threadIdx.x;
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

                if ( j >= 0 && j < numColumns )
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
            int k = threadIdx.x;
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

                if ( j >= 0 && j < numColumns )
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
            int k = threadIdx.x;
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
            int k = threadIdx.x;
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

                if ( j >= 0 && j < numColumns )
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
            int k = threadIdx.x;
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

                if ( j >= 0 && j < numColumns )
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
                    const ValueType diaValues[],
                    SyncToken* syncToken )
    {
        SCAI_REGION( "CUDA.DIA.normalGEMV" )

        SCAI_LOG_INFO( logger, "normalGEMV<" << getScalarType<ValueType>() << ">"
                        << " result[ " << numRows << "] = " << alpha
                        << " * A( #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

        const IndexType blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        SCAI_CHECK_CUDA_ACCESS

        cudaStream_t stream = 0;

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
        }

        const bool useSharedMem = CUDASettings::useSharedMem();
        const bool useTexture = CUDASettings::useTexture();

        SCAI_LOG_INFO( logger, "Start normal_gemv_kernel<" << getScalarType<ValueType>()
                        << "> <<< blockSize = " << blockSize << ", stream = " << stream
                        << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

        int sharedMemSize = 0;
        if ( useSharedMem )
        {
            sharedMemSize = numDiagonals * sizeof(int);
        }

        if ( useTexture )
        {
            vectorBindTexture( x );

            if ( !useSharedMem )
            {
                vectorBindTexture( diaOffsets );
            }

            if( useSharedMem )
            {
                if( alpha == 1 && beta == 1 )
                {
                    normal_gemv_kernel_alpha_one_beta_one<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 1 && beta == 0 )
                {
                    normal_gemv_kernel_alpha_one_beta_zero<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 0 && beta == 1 )
                {
                    assign_kernel<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, y, numRows );
                }
                else if ( alpha == 1 )
                {
                    normal_gemv_kernel_alpha_one<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 0 )
                {
                    normal_gemv_kernel_alpha_zero<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( beta == 1 )
                {
                    normal_gemv_kernel_beta_one<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( beta == 0 )
                {
                    normal_gemv_kernel_beta_zero<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else
                {
                    normal_gemv_kernel<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
            }
            else
            {
                if( alpha == 1 && beta == 1 )
                {
                    normal_gemv_kernel_alpha_one_beta_one<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 1 && beta == 0 )
                {
                    normal_gemv_kernel_alpha_one_beta_zero<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 0 && beta == 1 )
                {
                    assign_kernel<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, y, numRows );
                }
                else if ( alpha == 1 )
                {
                    normal_gemv_kernel_alpha_one<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 0 )
                {
                    normal_gemv_kernel_alpha_zero<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( beta == 1 )
                {
                    normal_gemv_kernel_beta_one<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( beta == 0 )
                {
                    normal_gemv_kernel_beta_zero<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else
                {
                    normal_gemv_kernel<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
            }
        }
        else
        {
            if( useSharedMem )
            {
                if( alpha == 1 && beta == 1 )
                {
                    normal_gemv_kernel_alpha_one_beta_one<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 1 && beta == 0 )
                {
                    normal_gemv_kernel_alpha_one_beta_zero<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 0 && beta == 1 )
                {
                    assign_kernel<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, y, numRows );
                }
                else if ( alpha == 1 )
                {
                    normal_gemv_kernel_alpha_one<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 0 )
                {
                    normal_gemv_kernel_alpha_zero<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( beta == 1 )
                {
                    normal_gemv_kernel_beta_one<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( beta == 0 )
                {
                    normal_gemv_kernel_beta_zero<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else
                {
                    normal_gemv_kernel<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
            }
            else
            {
                if( alpha == 1 && beta == 1 )
                {
                    normal_gemv_kernel_alpha_one_beta_one<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 1 && beta == 0 )
                {
                    normal_gemv_kernel_alpha_one_beta_zero<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 0 && beta == 1 )
                {
                    assign_kernel<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, y, numRows );
                }
                else if ( alpha == 1 )
                {
                    normal_gemv_kernel_alpha_one<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( alpha == 0 )
                {
                    normal_gemv_kernel_alpha_zero<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( beta == 1 )
                {
                    normal_gemv_kernel_beta_one<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else if ( beta == 0 )
                {
                    normal_gemv_kernel_beta_zero<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                    result, x, y, alpha, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
                }
                else
                {
                    normal_gemv_kernel<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
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
            int k = threadIdx.x;
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
            ValueType summand = beta * y[k];

            ValueType temp = 0.0;

            for ( IndexType ii = 0; ii < numDiagonals; ii++ )
            {
                IndexType i = k - fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, ii );

                if ( i >= 0 && i < numRows )
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
                    const ValueType diaValues[],
                    SyncToken* syncToken )
    {
        SCAI_REGION( "CUDA.DIA.normalGEVM" )

        SCAI_LOG_INFO( logger, "normalGEVM<" << getScalarType<ValueType>() << ">"
                        << " result[ " << numRows << "] = " << alpha
                        << " * A( #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

        const IndexType blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numColumns, dimBlock.x );

        SCAI_CHECK_CUDA_ACCESS

        cudaStream_t stream = 0;

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
        }

        const bool useSharedMem = CUDASettings::useSharedMem();
        const bool useTexture = CUDASettings::useTexture();

        SCAI_LOG_INFO( logger, "Start normal_gevm_kernel<" << getScalarType<ValueType>()
                        << "> <<< blockSize = " << blockSize << ", stream = " << stream
                        << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

        int sharedMemSize = 0;
        if ( useSharedMem )
        {
            sharedMemSize = numDiagonals * sizeof(int);
        }

        if ( useTexture )
        {
            vectorBindTexture( x );

            if ( !useSharedMem )
            {
                // @ToDo: be careful, some CUDA devices do not support multiple bind textures, e.g. GeForce 460
                vectorBindTexture( diaOffsets );
            }

            if( useSharedMem )
            {
                normal_gevm_kernel<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else
            {
                normal_gevm_kernel<ValueType, true, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
        }
        else
        {
            if( useSharedMem )
            {
                normal_gevm_kernel<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                                result, x, y, alpha, beta, diaValues, diaOffsets, numRows, numColumns, numDiagonals );
            }
            else
            {
                normal_gevm_kernel<ValueType, false, false><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
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

    void CUDADIAUtils::setInterface( DIAUtilsInterface& DIAUtils )
    {
        SCAI_LOG_INFO( logger, "set DIA routines for CUDA in Interface" )

#define LAMA_DIA_UTILS_REGISTER(z, I, _)                                                       \
    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEMV, ARITHMETIC_CUDA_TYPE_##I )                \
    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEVM, ARITHMETIC_CUDA_TYPE_##I )                \
                                                                                         
    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT, LAMA_DIA_UTILS_REGISTER, _ )

#undef LAMA_DIA_UTILS_REGISTER

    }

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDADIAUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( context::CUDA );
    setInterface( interface.DIAUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDADIAUtils::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
