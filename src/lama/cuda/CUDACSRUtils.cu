/**
 * @file CUDACSRUtils.cu
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
 * @brief Implementation of CSR utilities with CUDA
 * @author Bea Hornef, Thomas Brandes, Jiri Kraus
 * @date 04.07.2012
 * @since 1.0.0
 */

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDACSRUtils.hpp>
#include <lama/cuda/CUDASettings.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>

// macros
#include <lama/macros/unused.hpp>

#include <cuda.h>
#include <cusparse.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <lama/tracing.hpp>

#include <lama/ContextFactory.hpp>

// thrust
#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>

//#include "cuPrintf.cuh"
//#include "cuPrintf.cu"

//#define CUDA_CAP_20

// TODO: defines for matrix multiplication, should be removed later!
#define STEP 1
#define NUM_BLOCKS 64
#define HASH_TABLE_SIZE 2048

#define NUM_THREADS (STEP*32)
#define NUM_WARPS NUM_THREADS/32
#define MAX_HASH_TRIES 96
#define HASH_A 684
#define HASH_B 46165
#define HASH_P 88651
#define HASH_C0 1
#define HASH_C1 1
#define SIZE_LOCAL_HASHTABLE 1024

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDACSRUtils::logger, "CUDA.CSRUtils" )

IndexType CUDACSRUtils::sizes2offsets( IndexType array[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "sizes2offsets " << " #n = " << n )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<IndexType> array_ptr( const_cast<IndexType*>( array ) );
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
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

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
    LAMA_LOG_INFO( logger, "offsets2sizes " << " #n = " << n )

    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    offsets2sizes_kernel<<<dimGrid, dimBlock>>>( sizes, offsets, n );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "offsets2sizes" )
}

/** Correct handle for cusparse is set due to the context. */

extern cusparseHandle_t CUDAContext_cusparseHandle;

/* --------------------------------------------------------------------------- */
/*     hasDiagonalProperty                                                     */
/* --------------------------------------------------------------------------- */

template<typename T>
struct identic_functor
{
    __host__ __device__
    double operator()( thrust::tuple<T,T> x )
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

    //this is the actual check
    if ( ja[ia[i]] != i )
    {
        *hasProperty = false;
    }
}

bool CUDACSRUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrIA[], const IndexType csrJA[] )
{
    if ( numDiagonals == 0 )
    {
        return false;
    }

    //make grid
    const int blockSize = 256;
    dim3 dimGrid( ( numDiagonals - 1 ) / blockSize + 1, 1, 1 ); // = makeGrid( numDiagonals, blockSize );
    dim3 dimBlock( blockSize, 1, 1 );

    LAMA_CHECK_CUDA_ACCESS

    bool* d_hasProperty;
    bool hasProperty;

    LAMA_CUDA_RT_CALL( cudaMalloc( (void**)&d_hasProperty, sizeof(bool) ),
                       "allocate 4 bytes on the device for the result of hasDiagonalProperty_kernel" )
    LAMA_CUDA_RT_CALL( cudaMemset( d_hasProperty, 1, sizeof(bool) ), "memset bool hasProperty = true" )

    hasDiagonalProperty_kernel<<<dimGrid, dimBlock>>>( numDiagonals, csrIA, csrJA, d_hasProperty );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "hasDiagonalProperty failed: are ia and ja correct?" )

    LAMA_CUDA_RT_CALL( cudaMemcpy( &hasProperty, d_hasProperty, sizeof(bool), cudaMemcpyDeviceToHost ),
                       "copy the result of hasDiagonalProperty_kernel to host" )

    return hasProperty;
}

/* --------------------------------------------------------------------------- */
/*     Template specialization convertCSR2CSC<float>                           */
/* --------------------------------------------------------------------------- */

template<>
void CUDACSRUtils::convertCSR2CSC(
    IndexType cscIA[],
    IndexType cscJA[],
    float cscValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const float csrValues[],
    int numRows,
    int numColumns,
    int /* numValues */)
{
    LAMA_LOG_INFO( logger,
                   "convertCSR2CSC<float> -> cusparseScsr2csc" << ", matrix size = " << numRows << " x " << numColumns )

    LAMA_CUSPARSE_CALL(
        cusparseScsr2csc( CUDAContext_cusparseHandle, numRows, numColumns, csrValues, csrIA, csrJA, cscValues, cscJA, cscIA, 1, CUSPARSE_INDEX_BASE_ZERO ),
        "convertCSR2SCC<float>" )

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "convertCSR2CSC" )
}

/* --------------------------------------------------------------------------- */
/*     Template specialization convertCSR2CSC<double>                          */
/* --------------------------------------------------------------------------- */

template<>
void CUDACSRUtils::convertCSR2CSC(
    IndexType cscIA[],
    IndexType cscJA[],
    double cscValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const double csrValues[],
    int numRows,
    int numColumns,
    int /* numValues */)
{
    LAMA_LOG_INFO( logger,
                   "convertCSR2CSC<double> -> cusparseDcsr2csc" << ", matrix size = " << numRows << " x " << numColumns )

    LAMA_CUSPARSE_CALL( cusparseDcsr2csc( CUDAContext_cusparseHandle, numRows, numColumns, 
                                          csrValues, csrIA, csrJA, cscValues, cscJA, cscIA, 1, 
                                          CUSPARSE_INDEX_BASE_ZERO ),
                        "convertCSR2SCC<double>" )

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "convertCSR2CSC" )
}

/* --------------------------------------------------------------------------- */
/*     Using indirectly accessed vector in texture                             */
/* --------------------------------------------------------------------------- */

texture<float,1> textureFloatXRef;

texture<int2,1> textureDoubleXRef;

__inline__ void csrBindTexture( const float* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, textureFloatXRef, vector ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
}

__inline__ void csrBindTexture( const double* vector )
{
    LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, textureDoubleXRef, vector ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
}

__inline__ void csrUnbindTexture( const float* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( textureFloatXRef ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" )
}

__inline__ void csrUnbindTexture( const double* )
{
    LAMA_CUDA_RT_CALL( cudaUnbindTexture( textureDoubleXRef ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" )
}

template<typename T, bool useTexture>
__inline__  __device__ T fetch_x_i( const T* const x, const int i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetch_x_i<float, true>( const float* const, const int i )
{
    return tex1Dfetch( textureFloatXRef, i );
}

template<>
__inline__ __device__
double fetch_x_i<double, true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( textureDoubleXRef, i );
    return __hiloint2double( v.y, v.x );
}

/* --------------------------------------------------------------------------- */

template<typename T, bool useTexture>
__global__
void normal_gemv_kernel(
    T* result,
    int n,
    const T alpha,
    const T* csrValues,
    const int* csrIA,
    const int* csrJA,
    const T* x_d,
    const T beta,
    const T* y_d )
{
    // result = alpha * A * x_d + beta * y_d

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        T summand = beta * y_d[i];
        const int rowStart = csrIA[i];
        const int rowEnd = csrIA[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetch_x_i<T, useTexture>( x_d, csrJA[jj] );
        }

        result[i] = alpha * value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename T, bool useTexture>
__global__
void sparse_gemv_kernel(
    int n,
    const T alpha,
    const T* csrValues,
    const int* csrIA,
    const int* csrJA,
    const T* x_d,
    const IndexType* rows,
    T* result )
{
    const int ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( ii < n )
    {
        IndexType i = rows[ii];
        const int rowStart = csrIA[i];
        const int rowEnd = csrIA[i + 1];
        T value = 0.0;

        for ( int jj = rowStart; jj < rowEnd; ++jj )
        {
            value += csrValues[jj] * fetch_x_i<T, useTexture>( x_d, csrJA[jj] );
        }

        result[i] += alpha * value;
    }
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
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "normalGEMV<" << typeid(ValueType).name() << ">" << ", #rows = " << numRows )

    LAMA_LOG_INFO( logger,
                   "alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", result = " << result )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    const int blockSize = 256;

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
        useTexture = false;
    }

    LAMA_LOG_INFO( logger, "Start csr_normal_gemv_kernel<" << typeid( ValueType ).name()
                           << ", useTexture = " << useTexture << ">" );

    if ( useTexture )
    {
        csrBindTexture( x );

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        normal_gemv_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
                    ( result, numRows, alpha, csrValues, csrIA, csrJA, x, beta, y );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        normal_gemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
                    ( result, numRows, alpha, csrValues, csrIA, csrJA, x, beta, y );
    }

    // ToDo: implement and choose more efficient variants for special
    //       cases: alpha = 1.0, alpha = -1.0, beta = 1.0, beta = 0.0, beta = -1.0

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "normalGEMV, stream = " << stream )
        LAMA_LOG_INFO( logger, "normalGEMV<" << typeid(ValueType).name() << "> synchronized" )
    }

    if ( useTexture )
    {
        csrUnbindTexture( x );
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
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "sparseGEMV<" << typeid(ValueType).name() << ">" << ", #non-zero rows = " << numNonZeroRows )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    const int blockSize = 256;

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numNonZeroRows, dimBlock.x );

    sparse_gemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
                    ( numNonZeroRows, alpha, csrValues, csrIA, csrJA, x, rowIndexes, result );

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "sparseGEMV, stream = " << stream )
        LAMA_LOG_INFO( logger, "sparseGEMV<" << typeid(ValueType).name() << "> synchronized" )
    }
}

/* --------------------------------------------------------------------------- */
/*                          Jacobi                                             */
/* --------------------------------------------------------------------------- */

template<typename T,bool useTexture>
__global__
void csr_jacobi_kernel(
    const int* const csrIA,
    const int* const csrJA,
    const T* const csrValues,
    const int numRows,
    const T* const rhs,
    T* const solution,
    const T* const oldSolution,
    const T omega )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( i < numRows )
    {
        T temp = rhs[i];
        const int rowStart = csrIA[i];
        const int rowEnd = csrIA[i + 1];
        const T diag = csrValues[rowStart];
        for ( int jj = rowStart + 1; jj < rowEnd; ++jj )
        {
            temp -= csrValues[jj] * fetch_x_i<T,useTexture>( oldSolution, csrJA[jj] );
        }
        if ( omega == 0.5 )
        {
            solution[i] = omega * ( fetch_x_i<T,useTexture>( oldSolution, i ) + temp / diag );
        }
        else if ( omega == 1.0 )
        {
            solution[i] = temp / diag;
        }
        else
        {
            solution[i] = omega * ( temp / diag ) + ( 1.0 - omega ) * fetch_x_i<T,useTexture>( oldSolution, i );
        }
    }
}

template<typename T>
__inline__  __device__ T getSharedValue( T* shared, const T* const value, const int index )
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
template<typename T>
struct SharedMemory
{
    //! @brief Return a pointer to the runtime-sized shared memory array.
    //! @returns Pointer to runtime-sized shared memory array
    __device__
    T* getPointer()
    {
        extern __device__ void Error_UnsupportedType(); // Ensure that we won't compile any un-specialized types
        Error_UnsupportedType();
        return (T*) 0;
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
template<typename T>
__global__ void csr_alternate_jacobi_kernel(
    const int* const csrIA,
    const int* const csrJA,
    const T* const csrValues,
    const int numRows,
    const T* const rhs,
    T* const solution,
    const T* const oldSolution,
    const T omega )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    SharedMemory<T> smem;
    T* shared = smem.getPointer();
    if ( i < numRows )
    {
        //this is the prefetch
        shared[threadIdx.x] = oldSolution[i];
        __syncthreads();

        T temp = rhs[i];
        const int rowStart = csrIA[i];
        const int rowEnd = csrIA[i + 1];
        const T diag = csrValues[rowStart];
        for ( int jj = rowStart + 1; jj < rowEnd; ++jj )
        {
            temp -= csrValues[jj] * getSharedValue<T>( shared, oldSolution, csrJA[jj] );
        }
        if ( omega == 0.5 )
        {
            solution[i] = omega * ( getSharedValue<T>( shared, oldSolution, i ) + temp / diag );
        }
        else if ( omega == 1.0 )
        {
            solution[i] = temp / diag;
        }
        else
        {
            solution[i] = omega * ( temp / diag ) + ( 1.0 - omega ) * getSharedValue<T>( shared, oldSolution, i );
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
    const IndexType numRows,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "jacobi, #rows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    bool useTexture = CUDASettings::useTexture();

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();

        useTexture = false;  // not yet supported
    }

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    LAMA_LOG_INFO( logger, "Start csr_jacobi_kernel<" << typeid( ValueType ).name()
                           << ", useTexture = " << useTexture << ">" );

    if ( useTexture )
    {
        csrBindTexture( oldSolution);

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobi_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        csr_jacobi_kernel <ValueType, true> <<<dimGrid, dimBlock, 0, stream>>>( csrIA, csrJA, csrValues, numRows,
                rhs, solution, oldSolution, omega );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobi_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        csr_jacobi_kernel<ValueType, false> <<<dimGrid, dimBlock, 0, stream>>>( csrIA, csrJA, csrValues, numRows, rhs,
                solution, oldSolution, omega );
    }

    if ( !syncToken )
    {
        cudaStreamSynchronize( stream );
    }

    if ( useTexture )
    {
        csrUnbindTexture( oldSolution );
    }
}

/* --------------------------------------------------------------------------- */
/*                          Jacobi halo                                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType,bool useTexture>
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
            temp += haloValues[jj] * fetch_x_i<ValueType, useTexture>( oldSolution, haloJA[jj] );
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
    LAMA_LOG_INFO( logger, "jacobiHalo, #non-empty rows = " << numNonEmptyRows )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 128;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numNonEmptyRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    useTexture = false;

    if ( useTexture )
    {
        csrBindTexture( oldSolution );

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHalo_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        csr_jacobiHalo_kernel <ValueType, true> <<<dimGrid, dimBlock>>>( solution, localIA, localValues, haloIA,
                haloJA, haloValues, haloRowIndexes,
                numNonEmptyRows, oldSolution, omega );
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHalo_kernel<ValueType, false>, cudaFuncCachePreferL1),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        csr_jacobiHalo_kernel<ValueType, false> <<<dimGrid, dimBlock>>>( solution, localIA, localValues, haloIA,
                haloJA, haloValues, haloRowIndexes, numNonEmptyRows,
                oldSolution, omega );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED" )
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), "LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED" )

    if ( useTexture )
    {
        csrUnbindTexture( oldSolution );
    }
}

/* --------------------------------------------------------------------------- */
/*                          Jacobi halo with diagonal array                    */
/* --------------------------------------------------------------------------- */


template<typename ValueType,bool useTexture>
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
            temp += haloValues[jj] * fetch_x_i<ValueType,useTexture>( oldSolution, haloJA[jj] );
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
    LAMA_LOG_INFO( logger, "jacobiHaloWithDiag, #non-empty rows = " << numNonEmptyRows )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 128;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numNonEmptyRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    useTexture = false;

    if ( useTexture )
    {
        csrBindTexture( oldSolution );

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHaloWithDiag_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr_jacobiHaloWithDiag_kernel<ValueType, false>, cudaFuncCachePreferL1),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

    }

    if ( useTexture )
    {
        csr_jacobiHaloWithDiag_kernel <ValueType, true> <<<dimGrid, dimBlock>>>( solution, localDiagValues, haloIA,
                haloJA, haloValues, haloRowIndexes,
                numNonEmptyRows, oldSolution, omega );
    }
    else
    {
        csr_jacobiHaloWithDiag_kernel<ValueType, false> <<<dimGrid, dimBlock>>>( solution, localDiagValues, haloIA,
                haloJA, haloValues, haloRowIndexes, numNonEmptyRows,
                oldSolution, omega );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_CSRJACOBIHALOWITHDIAG_CUDAKERNEL_FAILED" )
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), "LAMA_STATUS_CSRJACOBIHALOWITHDIAG_CUDAKERNEL_FAILED" )

    if ( useTexture )
    {
        csrUnbindTexture( oldSolution );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             helper                                                                 */
/* ------------------------------------------------------------------------------------------------------------------ */

__device__ __inline__ IndexType getNumActiveThreads(
    IndexType aColIt,
    IndexType aColEnd,
    const IndexType *aIA,
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
    IndexType *cIa,
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const IndexType *aIa,
    const IndexType *aJa,
    const IndexType *bIa,
    const IndexType *bJa )
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
                IndexType colA = aColIt < aColEnd ? aJa[aColIt] : -1;
                IndexType end = getNumActiveThreads( aColIt, aColEnd, aIa, rowIt, aColItOffset );

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
                        IndexType colB = ( bColIt + bColItOffset ) < bColEnd ? bJa[bColIt + bColItOffset] : -1;
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
    LAMA_REGION( "CUDA.CSR.matrixAddSizes" )

    LAMA_LOG_INFO(
        logger,
        "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

    LAMA_CHECK_CUDA_ACCESS

// Reset cIa
    thrust::device_ptr<IndexType> cIaPtr( cIa );
    thrust::fill( cIaPtr, cIaPtr + numRows, 0 );

// TODO: Check if diagonal property needs special attention
    matrixAddSizesKernel<NUM_WARPS><<<NUM_BLOCKS, NUM_THREADS>>>( cIa, numRows, numColumns, diagonalProperty,
            aIa, aJa, bIa, bJa );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR

// Convert sizes array to offset array
    thrust::exclusive_scan( cIaPtr, cIaPtr + numRows + 1, cIaPtr );

// Copy numValues from cIa to Host
// TODO: use cuMem cpy
    thrust::device_ptr<IndexType> iaPtr( cIa );
    thrust::host_vector<IndexType> numValues( iaPtr + numRows, iaPtr + numRows + 1 );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR

// TODO: write it!
    return numValues[0];
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixMultiplySizes                                                    */
/* ------------------------------------------------------------------------------------------------------------------ */

template<int nWarps>
__global__
void matrixMultiplySizesKernel(
    const IndexType *aIA,
    const IndexType *aJA,
    const IndexType *bIA,
    const IndexType *bJA,
    IndexType *cIA,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType rowOffset,
    IndexType* hashTables,
    IndexType numElementsHashTable,
    bool* hashError,
    bool diagonalProperty )
{
//TODO: We need the numColumns of the matrix in the kernel to use diagonalProperty correct in cases of not square matrices

    __shared__ volatile IndexType sColA[nWarps];
    __shared__ IndexType sHashTableJa[nWarps * SIZE_LOCAL_HASHTABLE];
    __shared__ volatile bool sGlobalHashTableAccessed[nWarps];

    IndexType localWarpId = threadIdx.x / warpSize;
    IndexType globalWarpId = ( blockIdx.x * blockDim.x + threadIdx.x ) / warpSize;
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;
//IndexType numWarpsLocal  = blockDim.x / warpSize;
    IndexType numWarpsGlobal = ( blockDim.x * gridDim.x ) / warpSize;
    IndexType aRowIt = globalWarpId;

//    // TODO: Just if we launch kernels for packs of lines
    aRowIt += rowOffset;

    IndexType colB;

// Set global hash table to accessed mode, so that it gets erased in the first iteration!
    sGlobalHashTableAccessed[localWarpId] = true;

// Reset hashError
    *hashError = false;

    IndexType* hashTableIndexes = ( (IndexType*) hashTables );

    for ( ; __any( aRowIt < numRows ); aRowIt += numWarpsGlobal )
    {
        if ( aRowIt < numRows )
        {
            if ( diagonalProperty && aRowIt >= numColumns )
            {
                diagonalProperty = false;
            }

            IndexType aColIt = aIA[aRowIt] + laneId;
            IndexType aColEnd = aIA[aRowIt + 1];

            IndexType hashTableOffset = globalWarpId * numElementsHashTable;

// STEP 1: Set Initial data in hash tables
// TODO: recheck loop condition
            for ( IndexType i = 0; i < SIZE_LOCAL_HASHTABLE; i += warpSize )
            {
                if ( i + laneId < SIZE_LOCAL_HASHTABLE )
                {
                    sHashTableJa[i + laneId] = -1;
                }
            }

// Set first index of hashtable to row (diagonal property)
            if ( laneId == 0 && diagonalProperty )
            {
                sHashTableJa[0] = aRowIt;
                cIA[aRowIt]++;
            }

            for ( IndexType offset = 0; __any( aColIt < aColEnd ); aColIt += warpSize, offset += warpSize )
            {
                IndexType colA = aColIt < aColEnd ? aJA[aColIt] : -1;

                IndexType end = getNumActiveThreads( aColIt, aColEnd, aIA, aRowIt, offset );

                for ( IndexType k = 0; k < end && k < warpSize; k++ )
                {
                    if ( laneId == k )
                    {
                        sColA[localWarpId] = colA;
                    }

                    IndexType bColIt = bIA[sColA[localWarpId]] + laneId;
                    IndexType bColEnd = bIA[sColA[localWarpId] + 1];

                    for ( ; __any( bColIt < bColEnd ); bColIt += warpSize )
                    {
                        colB = bColIt < bColEnd ? bJA[bColIt] : -1;

                        if ( colB != -1 && ( !diagonalProperty || colB != aRowIt ) )
                        {
                            bool inserted = false;
                            unsigned int fx = HASH_A * colB;
                            unsigned int gx = ( fx + HASH_B ) % HASH_P;

                            for ( IndexType i = 0; i < MAX_HASH_TRIES; i++ )
                            {
                                //TODO: diagonal property
                                IndexType hash;
                                if ( diagonalProperty )
                                {
                                    hash = ( ( gx + HASH_C0 * i + HASH_C1 * (IndexType) i * i ) % ( SIZE_LOCAL_HASHTABLE - 1 ) ) + 1;
                                }
                                else
                                {
                                    hash = ( gx + HASH_C0 * i + HASH_C1 * (IndexType) i * i ) % SIZE_LOCAL_HASHTABLE;
                                }

                                IndexType val = atomicCAS( &sHashTableJa[hash], -1, colB );

                                if ( val == -1 )
                                {
                                    atomicAdd( &cIA[aRowIt], 1 );
                                    inserted = true;
                                    break;
                                }
                                if ( val == colB )
                                {
                                    inserted = true;
                                    break;
                                }
                            }

                            if ( !inserted )
                            {
                                for ( IndexType i = 0; i < MAX_HASH_TRIES; i++ )
                                {
                                    IndexType hash;
                                    if ( diagonalProperty )
                                    {
                                        hash = ( ( gx + HASH_C0 * i + HASH_C1 * (IndexType) i * i ) % ( numElementsHashTable - 1 ) )
                                               + 1;
                                    }
                                    else
                                    {
                                        hash = ( gx + HASH_C0 * i + HASH_C1 * (IndexType) i * i ) % ( numElementsHashTable );
                                    }
                                    IndexType val = atomicCAS( &hashTableIndexes[hashTableOffset + hash], -1, colB );

                                    if ( val == -1 )
                                    {
                                        inserted = true;
                                        atomicAdd( &cIA[aRowIt], 1 );
                                        sGlobalHashTableAccessed[localWarpId] = true;
                                        break;
                                    }
                                    if ( val == colB )
                                    {
                                        inserted = true;
                                        break;
                                    }
                                }
                            }
                            if ( !inserted )
                            {
                                // Wert konnte nicht in hashtable eingefügt werden
                                *hashError = true;
                            }
                        }
                    }
                }
            }

// if global hashTable was used, clean it!
            if ( sGlobalHashTableAccessed[localWarpId] )
            {
                for ( IndexType i = 0; i < numElementsHashTable; i += warpSize )
                {
                    if ( i + laneId < numElementsHashTable )
                    {
                        hashTableIndexes[hashTableOffset + i + laneId] = -1;
                    }
                }
                sGlobalHashTableAccessed[localWarpId] = false;
            }
        }
    }
}

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
    LAMA_REGION( "CUDA.CSR.matrixMultiplySizes" )

    LAMA_LOG_INFO(
        logger,
        "matrixMutliplySizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

    LAMA_CHECK_CUDA_ACCESS
// Reset cIa

    size_t free;
    size_t total;
    cuMemGetInfo( &free, &total );

//    std::cout << "free memory: " << free / 1024 / 1024 << "mb, total memory: " << total / 1024 / 1024 << "mb" << std::endl;

    thrust::device_ptr<IndexType> cIaPtr( cIa );
    thrust::fill( cIaPtr, cIaPtr + numRows, 0 );

    const unsigned int initialHashTableSize = 1024;

    unsigned int hashTableAllocatedBytes = NUM_BLOCKS * initialHashTableSize * sizeof( IndexType );

// Allocate hashTable and hashError Flag
    ContextPtr loc = ContextFactory::getContext( Context::CUDA );
// TODO: be carefull with NUM_BLOCKS here!
    IndexType* hashTable = ( IndexType* ) loc->allocate( hashTableAllocatedBytes );
    bool* hashError = (bool*) loc->allocate( sizeof(bool) );

// Reset hashTable
    thrust::device_ptr<IndexType> hashTablesPtr( hashTable );
    thrust::fill( hashTablesPtr, hashTablesPtr + NUM_BLOCKS * initialHashTableSize, -1 );

    bool hashErrorHost;
    unsigned int hashTableSize = initialHashTableSize;
    IndexType rowsPerLaunch = NUM_BLOCKS;
    IndexType iterations = 90;
    IndexType rows = std::ceil( std::ceil( numRows / (double) rowsPerLaunch ) / iterations );

    for ( IndexType i = 0; i < iterations; i++ )
    {
        IndexType startRow = rows * i * rowsPerLaunch;
        IndexType endRow = std::min( rows * ( i + 1 ) * rowsPerLaunch, numRows );
        if ( startRow >= numRows )
        {
            break;
        }

//        cudaPrintfInit();
        matrixMultiplySizesKernel<NUM_WARPS><<<NUM_BLOCKS, NUM_THREADS>>>( aIa, aJa, bIa, bJa, cIa, endRow, numColumns,
                startRow, hashTable, hashTableSize,
                hashError, diagonalProperty );
//        cudaPrintfDisplay(stdout, true);
//        cudaPrintfEnd();

        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
        ;

        cudaMemcpy( &hashErrorHost, hashError, sizeof(bool), cudaMemcpyDeviceToHost );

        if ( hashErrorHost )
        {
// repeat iteration:
            i--;

// Free old hashTable
            loc->free( (void*) hashTable, hashTableAllocatedBytes );

// Resize hashTable
            hashTableSize           *= 2;
            hashTableAllocatedBytes *= 2;
// TODO: be carefull with NUM_BLOCKS here!
            hashTable = (IndexType*) loc->allocate( hashTableAllocatedBytes );

// Reset new hashTable
            thrust::device_ptr < IndexType > hashTablesPtr( hashTable );
            thrust::fill( hashTablesPtr, hashTablesPtr + NUM_BLOCKS * hashTableSize, -1 );

// We need to clean up cIA again (for the crashed rows!)
            thrust::fill( cIaPtr + startRow, cIaPtr + endRow, 0 );
        }
    }

// Free hashTable and hashError
    loc->free( (void*) hashTable, hashTableAllocatedBytes );
    loc->free( (void*) hashError, sizeof(bool) );

//    thrust::device_ptr<IndexType> iaPtr ( cIa );
//    thrust::host_vector<IndexType> values ( iaPtr, iaPtr + numRows + 1 );
//    for(IndexType i = 0; i < numRows + 1; i++ )
//    {
//        std::cout << "cIA[" << i << "] = " << values[i] << std::endl;
//    }

// Convert sizes array to offset array
    thrust::exclusive_scan( cIaPtr, cIaPtr + numRows + 1, cIaPtr );

// Copy numValues from cIa to Host
// TODO: use cuMem cpy
    thrust::device_ptr<IndexType> iaPtr( cIa );
    thrust::host_vector<IndexType> numValues( iaPtr + numRows, iaPtr + numRows + 1 );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR

    return numValues[0];
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixAdd                                                              */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,int nWarps>
__global__
void matrixAddKernel(
    IndexType *cJA,
    ValueType *cValues,
    const IndexType *cIA,
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const ValueType alpha,
    const IndexType *aIA,
    const IndexType *aJA,
    const ValueType *aValues,
    const ValueType beta,
    const IndexType *bIA,
    const IndexType *bJA,
    const ValueType *bValues )
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
                IndexType colB = ( bColIt + bColOffset ) < bColEnd ? bJA[bColIt + bColOffset] : -1;
                ValueType valB = ( bColIt + bColOffset ) < bColEnd ? bValues[bColIt + bColOffset] : 0.0;
                if ( colB != -1 )
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
                IndexType colA = aColIt < aColEnd ? aJA[aColIt] : -1;
                ValueType valA = aColIt < aColEnd ? aValues[aColIt] : 0.0;
                IndexType end = getNumActiveThreads( aColIt, aColEnd, aIA, rowIt, aColItOffset );

                for ( IndexType k = 0; k < end && k < warpSize; k++ )
                {
                    if ( laneId == k )
                    {
                        sColA[localWarpId] = colA;
                        sValA[localWarpId] = valA;
                        sFoundJa[localWarpId] = -1;
                    }

                    for ( IndexType bColItOffset = 0; ( sFoundJa[localWarpId] == -1 ) && __any( ( bColIt + bColItOffset ) < bColEnd );
                            bColItOffset += warpSize )
                    {
                        IndexType colB = ( bColIt + bColItOffset ) < bColEnd ? bJA[bColIt + bColItOffset] : -1;
                        if ( sColA[localWarpId] == colB )
                        {
                            sFoundJa[localWarpId] = laneId + bColItOffset;
                        }
                    }
                    if ( laneId == 0 )
                    {
                        if ( sFoundJa[localWarpId] == -1 )
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
    LAMA_REGION( "CUDA.CSR.matrixAdd" )

    LAMA_LOG_INFO( logger, "matrixAdd for " << numRows << "x" << numColumns << " matrix" )

    LAMA_CHECK_CUDA_ACCESS

    matrixAddKernel<ValueType, NUM_WARPS><<<NUM_BLOCKS, NUM_THREADS>>>( cJA, cValues, cIA, numRows, numColumns,
            diagonalProperty, alpha, aIA, aJA, aValues, beta, bIA, bJA, bValues );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR
}
/* ------------------------------------------------------------------------------------------------------------------ */
/*                                             matrixMultiply                                                         */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,int nWarps>
__global__
void matrixMultiplyKernel(
    const IndexType *aIA,
    const IndexType *aJA,
    const ValueType *aValues,
    const IndexType *bIA,
    const IndexType *bJA,
    const ValueType *bValues,
    const IndexType *cIA,
    const ValueType alpha,
    IndexType *cJA,
    ValueType *cValues,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType rowOffset,
    void* hashTables,
    IndexType numElementsHashTable,
    bool* hashError,
    bool diagonalProperty )
{
//TODO: We need the numColumns of the matrix in the kernel to use diagonalProperty correct in cases of not square matrices

// TODO: Check for correct size!
    __shared__ volatile IndexType sColA[nWarps];
    __shared__ volatile ValueType sValA[nWarps];
    __shared__ volatile bool sGlobalHashTableAccessed[nWarps];

// TODO: rename hash tables to names that correlate to shared and global memory!
    __shared__ IndexType sHashTableJa[nWarps * SIZE_LOCAL_HASHTABLE];
    __shared__ ValueType sHashTableValues[nWarps * SIZE_LOCAL_HASHTABLE];

#if CUDA_ARCH < 20
    __shared__ volatile IndexType sBallot[nWarps];
#endif

    IndexType localWarpId = threadIdx.x / warpSize;
    IndexType globalWarpId = ( blockIdx.x * blockDim.x + threadIdx.x ) / warpSize;
    IndexType laneId = ( blockIdx.x * blockDim.x + threadIdx.x ) % warpSize;
//IndexType numWarpsLocal  = blockDim.x / warpSize;
    IndexType numWarpsGlobal = ( blockDim.x * gridDim.x ) / warpSize;
    IndexType aRowIt = globalWarpId;

    aRowIt += rowOffset;

// Reset hashError
    *hashError = false;

    IndexType* hashTableIndexes = ( (IndexType*) hashTables );
    ValueType* hashTableValues = (ValueType*) &hashTableIndexes[numWarpsGlobal * numElementsHashTable];

// Loop over all rows
    for ( ; __any( aRowIt < numRows ); aRowIt += numWarpsGlobal )
    {
// Check if this warp is in valid row
        if ( aRowIt < numRows )
        {
            if ( diagonalProperty && aRowIt >= numColumns )
            {
                diagonalProperty = false;
            }

            IndexType aColIt = aIA[aRowIt] + laneId;
            IndexType aColEnd = aIA[aRowIt + 1];

            IndexType hashTableOffset = globalWarpId * numElementsHashTable;

// STEP 1: Set Initial data in hash tables
// TODO: recheck loop condition
            for ( IndexType i = 0; i < SIZE_LOCAL_HASHTABLE; i += warpSize )
            {
                if ( i + laneId < SIZE_LOCAL_HASHTABLE )
                {
                    sHashTableJa[i + laneId] = -1;
                }
            }

// STEP 2: Calculate A x B in hash tables

// Set first index of hashtable to row (diagonal property)
            if ( diagonalProperty && laneId == 0 )
            {
                sHashTableJa[0] = aRowIt;
                sHashTableValues[0] = 0.0;
            }

            for ( IndexType offset = 0; __any( aColIt < aColEnd ); aColIt += warpSize, offset += warpSize )
            {
                IndexType colA = aColIt < aColEnd ? aJA[aColIt] : -1;
                ValueType valA = aColIt < aColEnd ? aValues[aColIt] : 0.0;

                IndexType end = getNumActiveThreads( aColIt, aColEnd, aIA, aRowIt, offset );

                for ( IndexType k = 0; k < end && k < warpSize; k++ )
                {
                    if ( laneId == k )
                    {
                        sColA[localWarpId] = colA;
                        sValA[localWarpId] = valA;
                    }

                    IndexType bColIt = bIA[sColA[localWarpId]] + laneId;
                    IndexType bColEnd = bIA[sColA[localWarpId] + 1];

                    for ( ; __any( bColIt < bColEnd ); bColIt += warpSize )
                    {
                        IndexType colB = bColIt < bColEnd ? bJA[bColIt] : -1;
                        ValueType valB = bColIt < bColEnd ? bValues[bColIt] : 0.0;

                        if ( colB != -1 )
                        {
                            // extra check for diagonal property
                            if ( diagonalProperty && colB == aRowIt )
                            {
                                sHashTableValues[0] += valB * sValA[localWarpId];
                            }
                            else
                            {
                                bool inserted = false;
                                unsigned int fx = HASH_A * colB;
                                unsigned int gx = ( fx + HASH_B ) % HASH_P;

                                // Hashing in shared memory
                                for ( IndexType i = 0; i < MAX_HASH_TRIES; i++ )
                                {
                                    //TODO: diagonal property
                                    IndexType hash;
                                    if ( diagonalProperty )
                                    {
                                        hash = ( ( gx + HASH_C0 * i + HASH_C1 * (IndexType) i * i ) % ( SIZE_LOCAL_HASHTABLE - 1 ) )
                                               + 1;
                                    }
                                    else
                                    {
                                        hash = ( gx + HASH_C0 * i + HASH_C1 * (IndexType) i * i ) % SIZE_LOCAL_HASHTABLE;
                                    }

                                    IndexType val = atomicCAS( &sHashTableJa[hash], -1, colB );

                                    if ( val == -1 )
                                    {
//                                        atomicAdd(&operations, 1);
                                        sHashTableValues[hash] = valB * sValA[localWarpId];
                                        inserted = true;
                                        break;
                                    }
                                    if ( val == colB )
                                    {
//                                        atomicAdd(&operations, 2);
                                        sHashTableValues[hash] += valB * sValA[localWarpId];
                                        inserted = true;
                                        break;
                                    }
                                }

                                // Hashing in global memory
                                if ( !inserted )
                                {
                                    for ( IndexType i = 0; i < MAX_HASH_TRIES; i++ )
                                    {
                                        //TODO: diagonal property
                                        IndexType hash;
                                        if ( diagonalProperty )
                                        {
                                            hash = ( ( gx + HASH_C0 * i + HASH_C1 * (IndexType) i * i ) % ( numElementsHashTable - 1 ) )
                                                   + 1;
                                        }
                                        else
                                        {
                                            hash = ( gx + HASH_C0 * i + HASH_C1 * (IndexType) i * i ) % numElementsHashTable;
                                        }

                                        IndexType val = atomicCAS( &hashTableIndexes[hashTableOffset + hash], -1, colB );

                                        if ( val == -1 )
                                        {
                                            //                                        atomicAdd(&operations, 1);
                                            hashTableValues[hashTableOffset + hash] = valB * sValA[localWarpId];
                                            inserted = true;
                                            sGlobalHashTableAccessed[localWarpId] = true;
                                            break;
                                        }
                                        if ( val == colB )
                                        {
                                            //                                        atomicAdd(&operations, 2);
                                            hashTableValues[hashTableOffset + hash] += valB * sValA[localWarpId];
                                            inserted = true;
                                            break;
                                        }
                                    }
                                }

                                if ( !inserted )
                                {
                                    // Wert konnte nicht in hashtable eingefügt werden
                                    *hashError = true;
                                }
                            }
                        }
                    }
                }
            }

// STEP 3: Copy hash Tables to cJa and cValues
// TODO: rename sColA => destinationOffset!

            sColA[localWarpId] = 0;
            IndexType rowOffset = cIA[aRowIt];
            for ( IndexType offset = 0; offset < SIZE_LOCAL_HASHTABLE; offset += warpSize )
            {
                if ( offset + laneId < SIZE_LOCAL_HASHTABLE )
                {
                    IndexType hashCol = sHashTableJa[offset + laneId];
                    ValueType hashVal = sHashTableValues[offset + laneId];

#if CUDA_ARCH >= 20
                    // TODO: be carefull here, ballot is warpsize Bit's long!
                    IndexType ballot = __ballot ( hashCol != -1 );
#else
                    if ( laneId == 0 )
                    {
                        sBallot[localWarpId] = 0;
                    }

                    if ( hashCol != -1 )
                    {
                        atomicOr( (int*) &sBallot[localWarpId], (int) ( 1 << laneId ) );
                    }
                    IndexType ballot = sBallot[localWarpId];
#endif

                    IndexType localOffset = __popc( ballot << ( warpSize - laneId ) );

                    if ( hashCol != -1 )
                    {

                        cJA[rowOffset + sColA[localWarpId] + localOffset] = hashCol;
                        cValues[rowOffset + sColA[localWarpId] + localOffset] = hashVal * alpha;
                    }

                    if ( laneId == 0 )
                    {
                        sColA[localWarpId] += __popc( ballot );
                    }
                }
            }

// copy global memory
            if ( sGlobalHashTableAccessed[localWarpId] )
            {
                for ( IndexType offset = 0; offset < numElementsHashTable; offset += warpSize )
                {
                    if ( offset + laneId < numElementsHashTable )
                    {
                        IndexType hashCol = hashTableIndexes[hashTableOffset + offset + laneId];
                        ValueType hashVal = hashTableValues[hashTableOffset + offset + laneId];

                        // Clean hashTable!
                        hashTableIndexes[hashTableOffset + offset + laneId] = -1;

#if CUDA_ARCH >= 20
                        // TODO: be carefull here, ballot is warpsize Bit's long!
                        IndexType ballot = __ballot ( hashCol != -1 );
#else
                        if ( laneId == 0 )
                        {
                            sBallot[localWarpId] = 0;
                        }

                        if ( hashCol != -1 )
                        {
                            atomicOr( (int*) &sBallot[localWarpId], (int) ( 1 << laneId ) );
                        }
                        IndexType ballot = sBallot[localWarpId];
#endif

                        IndexType localOffset = __popc( ballot << ( warpSize - laneId ) );

                        if ( hashCol != -1 )
                        {

                            cJA[rowOffset + sColA[localWarpId] + localOffset] = hashCol;
                            cValues[rowOffset + sColA[localWarpId] + localOffset] = hashVal * alpha;
                        }

                        if ( laneId == 0 )
                        {
                            sColA[localWarpId] += __popc( ballot );
                        }
                    }
                }
            }

//            // clean up global hashTable if used
//            if( sGlobalHashTableAccessed[localWarpId] )
//            {
//                for( IndexType i = 0 ; i < numElementsHashTable; i += warpSize )
//                {
//                    if( i + laneId < numElementsHashTable )
//                    {
//                        hashTableIndexes[hashTableOffset + i + laneId] = -1;
//                    }
//                }
//                sGlobalHashTableAccessed[localWarpId] = false;
//            }

        } // end check if this thread is in valid row
    } // end loop over rows
//    cuPrintf("operations: %i\n", operations);
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
    LAMA_REGION( "CUDA.CSR.matrixMultiply" )

    LAMA_LOG_INFO( logger, "matrixMultiply for " << numRows << "x" << numColumns << " matrix" )

    LAMA_CHECK_CUDA_ACCESS

    const unsigned int initialHashTableSize = 1024;

// Allocate hashTable and hashError Flag
    ContextPtr loc = ContextFactory::getContext( Context::CUDA );
// TODO: be carefull with NUM_BLOCKS here!

    unsigned int hashTableAllocatedBytes = NUM_BLOCKS * initialHashTableSize * sizeof( IndexType ) +
                                           NUM_BLOCKS * initialHashTableSize * sizeof( ValueType );

    IndexType* hashTable = (IndexType*) loc->allocate( hashTableAllocatedBytes );

    bool* hashError = (bool*) loc->allocate( sizeof(bool) );

// Reset hashTable
    thrust::device_ptr<IndexType> hashTablesPtr( hashTable );
    thrust::fill( hashTablesPtr, hashTablesPtr + initialHashTableSize, -1 );

    bool hashErrorHost;
    unsigned int hashTableSize = initialHashTableSize;
    IndexType rowsPerLaunch = NUM_BLOCKS;
    IndexType iterations = 90;
    IndexType rows = std::ceil( std::ceil( numRows / (double) rowsPerLaunch ) / iterations );
    for ( IndexType i = 0; i < iterations; i++ )
    {
        IndexType startRow = rows * i * rowsPerLaunch;
        IndexType endRow = std::min( rows * ( i + 1 ) * rowsPerLaunch, numRows );
        if ( startRow >= numRows )
        {
            break;
        }

//        cudaPrintfInit();
        matrixMultiplyKernel<ValueType, NUM_WARPS><<<NUM_BLOCKS, NUM_THREADS>>>( aIa, aJa, aValues, bIa, bJa, bValues, cIa,
                alpha, cJa, cValues, endRow, numColumns, startRow, hashTable, hashTableSize,
                hashError, diagonalProperty );
//        cudaPrintfDisplay(stdout, true);
//        cudaPrintfEnd();

        cudaStreamSynchronize( 0 );
        cudaMemcpy( &hashErrorHost, hashError, sizeof(bool), cudaMemcpyDeviceToHost );

        if ( hashErrorHost )
        {
// repeat iteration:
            i--;

// Free old hashTable
            loc->free( (void*) hashTable, hashTableAllocatedBytes );

// Resize hashTable
            hashTableSize *= 2;
            hashTableAllocatedBytes *= 2;

// TODO: be carefull with NUM_BLOCKS here!
            hashTable = (IndexType*) loc->allocate( hashTableAllocatedBytes );

// Reset new hashTable
            thrust::device_ptr < IndexType > hashTablesPtr( hashTable );
            thrust::fill( hashTablesPtr, hashTablesPtr + NUM_BLOCKS * hashTableSize, -1 );
        }
    }
    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR

// Free hashTable and hashError
    loc->free( (void*) hashTable, hashTableAllocatedBytes );
    loc->free( (void*) hashError, sizeof(bool) );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDACSRUtils::setInterface( CSRUtilsInterface& CSRUtils )
{
    LAMA_INTERFACE_REGISTER( CSRUtils, sizes2offsets )
    LAMA_INTERFACE_REGISTER( CSRUtils, offsets2sizes )
    LAMA_INTERFACE_REGISTER( CSRUtils, hasDiagonalProperty )
    LAMA_INTERFACE_REGISTER( CSRUtils, matrixAddSizes )
    LAMA_INTERFACE_REGISTER( CSRUtils, matrixMultiplySizes )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, convertCSR2CSC, float )
    LAMA_INTERFACE_REGISTER_T( CSRUtils, convertCSR2CSC, double )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( CSRUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, sparseGEMV, float )
    LAMA_INTERFACE_REGISTER_T( CSRUtils, sparseGEMV, double )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobi, float )
    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobi, double )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobiHalo, float )
    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobiHalo, double )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobiHaloWithDiag, float )
    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobiHaloWithDiag, double )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, matrixAdd, float )
    LAMA_INTERFACE_REGISTER_T( CSRUtils, matrixAdd, double )

    LAMA_INTERFACE_REGISTER_T( CSRUtils, matrixMultiply, float )
    LAMA_INTERFACE_REGISTER_T( CSRUtils, matrixMultiply, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDACSRUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.CSRUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDACSRUtils::initialized = registerInterface();


} // namespace lama
