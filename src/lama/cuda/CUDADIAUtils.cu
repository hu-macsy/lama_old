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

#include <lama/exception/LAMAAssert.hpp>

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDADIAUtils.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>
#include <lama/cuda/CUDASettings.hpp>
#include <lama/tracing.hpp>

// thrust
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

#include <boost/bind.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDADIAUtils::logger, "CUDA.DIAUtils" )

/* --------------------------------------------------------------------------- */

texture<float, 1> texDIAVectorSXref;

texture<int2, 1> texDIAVectorDXref;

texture<int, 1> texDIAVectorIref;

__inline__ void vectorDIABindTexture( const float* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texDIAVectorSXref, vector ), "bind float vector x to texture" )
}

__inline__ void vectorDIABindTexture( const double* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texDIAVectorDXref, vector ), "bind double vector x to texture" )
}

__inline__ void vectorDIABindTexture( const int* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texDIAVectorIref, vector ), "bind int vector x to texture" )
}

__inline__ void vectorDIAUnbindTexture( const float* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texDIAVectorSXref ), "unbind float vector x from texture" )
}

__inline__ void vectorDIAUnbindTexture( const double* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texDIAVectorDXref ), "unbind double vector x from texture" )
}

__inline__ void vectorDIAUnbindTexture( const int* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( texDIAVectorIref ), "unbind int vector x from texture" )
}

template<typename ValueType, bool useTexture>
__inline__ __device__ 
ValueType fetchDIAVectorX( const ValueType* const x, const int i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetchDIAVectorX<float, true>( const float* const, const int i )
{
    return tex1Dfetch( texDIAVectorSXref, i );
}

template<>
__inline__ __device__
double fetchDIAVectorX<double, true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texDIAVectorDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<>
__inline__ __device__
int fetchDIAVectorX<int, true>( const int* const, const int i )
{
    return tex1Dfetch( texDIAVectorIref, i );
}

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
    return fetchDIAVectorX<int, true>( offset_d, i );
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
__global__ void diagemvKernel(
    ValueType* result,
    const IndexType numRows,
    const ValueType alpha,
    const ValueType* x,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType* offsets_d,
    const ValueType* diagonalValues,
    const ValueType beta,
    const ValueType* y )
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
                temp += val * fetchDIAVectorX<ValueType, useTexture>( x, j );
            }
        }

        result[i] = alpha * temp + summand;
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
    class SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.DIA.normalGEMV" )

    LAMA_LOG_INFO( logger, "normalGEMV<" << Scalar::getType<ValueType>() << ">"
                            << " result[ " << numRows << "] = " << alpha
                            << " * A( #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

    const IndexType blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    const bool useSharedMem = CUDASettings::useSharedMem();
    const bool useTexture   = CUDASettings::useTexture();

    LAMA_LOG_INFO( logger, "Start diaGemvKernel<" << Scalar::getType<ValueType>()
                           << "> <<< blockSize = " << blockSize << ", stream = " << stream
                           << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

    if ( useTexture )
    {
        vectorDIABindTexture( x );

        if ( !useSharedMem )
        {
            vectorDIABindTexture( diaOffsets );

            diagemvKernel<ValueType, true, false><<< dimGrid, dimBlock, 0, stream >>>( 
                result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );
        }
        else
        {
            const int sharedMemSize = numDiagonals * sizeof(int);

            diagemvKernel<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>( 
                result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );
        }
    }
    else
    {
        if ( !useSharedMem )
        {
            diagemvKernel<ValueType, false, false><<< dimGrid, dimBlock, 0, stream >>>( 
                result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );
        }
        else
        {
            const int sharedMemSize = numDiagonals * sizeof(int);

            diagemvKernel<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>( 
                result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );
        }
    }

    if ( !syncToken )
    {
        // synchronize now, unbind used textures

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "normalGEMV for DIA" )

        if ( useTexture )
        {
            vectorDIAUnbindTexture( x );

            if ( !useSharedMem )
            {
                vectorDIAUnbindTexture( diaOffsets );
            }
        }
    }
    else
    {
        // synchronize by syncToken, delay unbind texture 

        if ( useTexture )
        {
            void ( *unbindV ) ( const ValueType* ) = &vectorDIAUnbindTexture;
            void ( *unbindI ) ( const IndexType* ) = &vectorDIAUnbindTexture;
    
            syncToken->pushRoutine( boost::bind( unbindV, x ) );

            if ( !useSharedMem )
            {
                syncToken->pushRoutine( boost::bind( unbindI, diaOffsets ) );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__ void diagevmKernel(
    ValueType* result,
    const IndexType numRows,
    const ValueType alpha,
    const ValueType* x,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType* offsets_d,
    const ValueType* diagonalValues,
    const ValueType beta,
    const ValueType* y )
{
    extern __shared__ IndexType offsets_sm[];

    /*if ( useSharedMem )
    {
        int k = threadIdx.x;
        while ( k < numDiagonals )
        {
            offsets_sm[k] = offsets_d[k];
            k += blockDim.x;
        }
        __syncthreads();
    }*/

    IndexType k = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( k < numColumns )
    {
        ValueType summand = 0.0;

        if ( beta != 0.0 )
        {
            summand = beta * y[k];
        }

        ValueType temp = 0.0;

        for( IndexType i = 0; i < numRows; ++k )
        {
            for ( IndexType ii = 0; ii < numDiagonals; ii++ )
            {
                //IndexType j = k + fetchOffset<useTexture, useSharedMem>( offsets_d, offsets_sm, ii );
                IndexType j = k + offsets_d[ ii ];

                if ( j >= 0 && j < numColumns )
                {
                    if( j == k )
                    {
                        temp += diagonalValues[ numRows * ii + k ] * fetchDIAVectorX<ValueType, useTexture>( x, i );
                    }
                }
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
    class SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.DIA.normalGEVM" )

    LAMA_LOG_INFO( logger, "normalGEVM<" << Scalar::getType<ValueType>() << ">"
                            << " result[ " << numRows << "] = " << alpha
                            << " * A( #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

    const IndexType blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    const bool useSharedMem = CUDASettings::useSharedMem();
    const bool useTexture   = CUDASettings::useTexture();

    LAMA_LOG_INFO( logger, "Start diaGevmKernel<" << Scalar::getType<ValueType>()
                           << "> <<< blockSize = " << blockSize << ", stream = " << stream
                           << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

    if ( useTexture )
    {
        vectorDIABindTexture( x );

        if ( !useSharedMem )
        {
            vectorDIABindTexture( diaOffsets );

            diagevmKernel<ValueType, true, false><<< dimGrid, dimBlock, 0, stream >>>(
                result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );
        }
        else
        {
            const int sharedMemSize = numDiagonals * sizeof(int);

            diagevmKernel<ValueType, true, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );
        }
    }
    else
    {
        if ( !useSharedMem )
        {
            diagevmKernel<ValueType, false, false><<< dimGrid, dimBlock, 0, stream >>>(
                result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );
        }
        else
        {
            const int sharedMemSize = numDiagonals * sizeof(int);

            diagevmKernel<ValueType, false, true><<< dimGrid, dimBlock, sharedMemSize, stream >>>(
                result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );
        }
    }

    if ( !syncToken )
    {
        // synchronize now, unbind used textures

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "normalGEMV for DIA" )

        if ( useTexture )
        {
            vectorDIAUnbindTexture( x );

            if ( !useSharedMem )
            {
                vectorDIAUnbindTexture( diaOffsets );
            }
        }
    }
    else
    {
        // synchronize by syncToken, delay unbind texture

        if ( useTexture )
        {
            void ( *unbindV ) ( const ValueType* ) = &vectorDIAUnbindTexture;
            void ( *unbindI ) ( const IndexType* ) = &vectorDIAUnbindTexture;

            syncToken->pushRoutine( boost::bind( unbindV, x ) );

            if ( !useSharedMem )
            {
                syncToken->pushRoutine( boost::bind( unbindI, diaOffsets ) );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

void CUDADIAUtils::setInterface( DIAUtilsInterface& DIAUtils )
{
    LAMA_LOG_INFO( logger, "set DIA routines for CUDA in Interface" )

    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEVM, float )
    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEVM, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDADIAUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.DIAUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDADIAUtils::initialized = registerInterface();

} // namespace lama
