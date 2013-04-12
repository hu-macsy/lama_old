/**
 * @file CUDACOOUtils.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Implementation of COO utilities with CUDA
 * @author Bea Hornef, Thomas Brandes
 * @date 04.07.2012
 * $Id$
 */

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDACOOUtils.hpp>

// thrust
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDACOOUtils::logger, "CUDA.COOUtils" )

/* --------------------------------------------------------------------------- */

template<typename T,bool useTexture>
__inline__     __device__ T fetch_COOx( const T* const x, const IndexType i )
{
    return x[i];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__ void cooInitKernel( ValueType* result, const IndexType numRows, const ValueType beta, const ValueType* y )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        result[i] = beta * y[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__ void cooGemvKernel(
    ValueType* result,
    const IndexType numRows,
    const ValueType alpha,
    const ValueType* x,
    const IndexType numValues,
    const IndexType* cooIA,
    const IndexType* cooJA,
    const ValueType* cooValues )
{
    const int k = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( k < numValues )
    {
        IndexType i = cooIA[k];
        IndexType j = cooJA[k];

        // we must use atomic updates as different threads might update same row i

        const ValueType resultUpdate = alpha * cooValues[k] * x[j];

        // ToDo: make it atomic
        result[i] += resultUpdate;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDACOOUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const IndexType numValues,
    class SyncToken* /* syncToken */)
{
    const IndexType block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    LAMA_CHECK_CUDA_ACCESS

    cooInitKernel<<< dimGrid, dimBlock, block_size * sizeof( IndexType )>>>
    ( result, numRows, beta, y );

    cooGemvKernel<<< dimGrid, dimBlock, block_size * sizeof( IndexType )>>>
    ( result, numRows, alpha, x, numValues, cooIA, cooJA, cooValues );

    LAMA_CHECK_CUDA_ERROR

    cudaStreamSynchronize( 0 );
}

/* --------------------------------------------------------------------------- */

void CUDACOOUtils::setInterface( COOUtilsInterface& COOUtils )
{
    LAMA_LOG_INFO( logger, "set COO routines for CUDA in Interface" )

    LAMA_INTERFACE_REGISTER_T( COOUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( COOUtils, normalGEMV, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDACOOUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.COOUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDACOOUtils::initialized = registerInterface();

} // namespace lama
