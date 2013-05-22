/**
 * @file CUDADIAUtils.cpp
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
 * $Id$
 */

#include <lama/exception/LAMAAssert.hpp>

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDADIAUtils.hpp>
#include <lama/tracing.hpp>

// thrust
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDADIAUtils::logger, "CUDA.DIAUtils" )

/* --------------------------------------------------------------------------- */

template<typename T,bool useTexture>
__inline__    __device__ T fetch_DIAx( const T* const x, const IndexType i )
{
    return x[i];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__ void diagemvpbvKernelSmall(
    ValueType* result,
    const IndexType numRows,
    const ValueType alpha,
    const ValueType* x,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType* diagonalOffsets,
    const ValueType* diagonalValues,
    const ValueType beta,
    const ValueType* y )
{
    extern __shared__ IndexType diagonalOffsetsShared[];

    // fill the shared array for diagonals

    if ( threadIdx.x < numDiagonals )
    {
        diagonalOffsetsShared[threadIdx.x] = diagonalOffsets[threadIdx.x];
    }

    __syncthreads();
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
            IndexType j = i + diagonalOffsetsShared[idiag];

            if ( j >= 0 && j < numColumns )
            {
                ValueType val = diagonalValues[numRows * idiag + i];
                temp += val * fetch_DIAx<ValueType,true>( x, j );
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
    class SyncToken* /* syncToken */)
{
    LAMA_REGION( "CUDA.DIA.normalGEMV" )

    const IndexType block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    LAMA_CHECK_CUDA_ACCESS
    ;
    diagemvpbvKernelSmall <<< dimGrid, dimBlock, block_size * sizeof( IndexType )>>>
    ( result, numRows, alpha, x, numColumns, numDiagonals, diaOffsets, diaValues, beta, y );

    LAMA_CHECK_CUDA_ERROR
    ;

    cudaStreamSynchronize( 0 );
}

/* --------------------------------------------------------------------------- */

void CUDADIAUtils::setInterface( DIAUtilsInterface& DIAUtils )
{
    LAMA_LOG_INFO( logger, "set DIA routines for CUDA in Interface" )

    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEMV, double )
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
