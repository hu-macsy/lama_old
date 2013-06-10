/**
 * @file CUDACOOUtils.cpp
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
 * @brief Implementation of COO utilities with CUDA
 * @author Bea Hornef, Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDACOOUtils.hpp>
#include <lama/tracing.hpp>

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

__device__ inline void atomicAddDouble( double* address, double val )
{
    unsigned long long int* address_as_ull =
       (unsigned long long int*) address;

    unsigned long long int old = *address_as_ull, assumed;

    do 
    {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));
    } while (assumed != old);
}

__device__ inline void atomicAddFloat( float* address, float val)

{
    int i_val = __float_as_int(val);

    int tmp0 = 0;

    int tmp1;

    while( ( tmp1 = atomicCAS( ( int * ) address, tmp0, i_val ) ) != tmp0 )

    {
        tmp0 = tmp1;
        i_val = __float_as_int(val + __int_as_float(tmp1));
    }
}

template<typename ValueType>
__global__ void cooGemvKernel(
    ValueType* result,
    const IndexType numRows,
    const ValueType alpha,
    const ValueType* x,
    const IndexType numValues,
    const IndexType* cooIA,
    const IndexType* cooJA,
    const ValueType* cooValues );

template<>
__global__ void cooGemvKernel(
    double* result,
    const IndexType numRows,
    const double alpha,
    const double* x,
    const IndexType numValues,
    const IndexType* cooIA,
    const IndexType* cooJA,
    const double* cooValues )
{
    const int k = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( k < numValues )
    {
        IndexType i = cooIA[k];
        IndexType j = cooJA[k];

        // we must use atomic updates as different threads might update same row i

        const double resultUpdate = alpha * cooValues[k] * x[j];

        // This solution is very slow, but works

        atomicAddDouble( &result[i], resultUpdate );
    }
}

template<>
__global__ void cooGemvKernel(
    float* result,
    const IndexType numRows,
    const float alpha,
    const float* x,
    const IndexType numValues,
    const IndexType* cooIA,
    const IndexType* cooJA,
    const float* cooValues )
{
    const int k = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( k < numValues )
    {
        IndexType i = cooIA[k];
        IndexType j = cooJA[k];

        // we must use atomic updates as different threads might update same row i

        const float resultUpdate = alpha * cooValues[k] * x[j];

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 200

        // new fast solution

        atomicAdd( &result[i], resultUpdate );
#else
        // old slow solution

        atomicAddFloat( &result[i], resultUpdate );
#endif

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
    LAMA_REGION( "CUDA.COO.normalGEMV" )

    LAMA_LOG_INFO( logger, "normalGEMV, #rows = " << numRows << ", #vals = " << numValues )

    const IndexType block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numValues, dimBlock.x );

    LAMA_CHECK_CUDA_ACCESS

    cooInitKernel<<< dimGrid, dimBlock>>>
    ( result, numRows, beta, y );

    cudaStreamSynchronize( 0 );

    dim3 dimGrid1 = makeGrid( numValues, dimBlock.x );

    cooGemvKernel<<< dimGrid1, dimBlock>>>
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
