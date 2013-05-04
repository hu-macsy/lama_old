/**
 * @file CUDAELLUtils.cu
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
 * @brief Implementation of ELL utilities with CUDA
 * @author Thomas Brandes
 * @date 04.07.2012
 * $Id$
 */

// lama cuda
#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAStreamSyncToken.hpp>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAELLUtils.hpp>
#include <lama/cuda/CUDAUtils.hpp>
#include <lama/cuda/CUDATexture.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/macros/unused.hpp>
#include <lama/tracing.hpp>

// cuda
#include <cuda.h>
#include <cuda_runtime.h>

// thrust
#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/gather.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/scatter.h>
#include <thrust/tuple.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDAELLUtils::logger, "CUDA.ELLUtils" )

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  thrust functors                                                   */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T>
struct greaterThan
{
    const T x;
    greaterThan( T _x )
        : x( _x )
    {
    }

    __host__ __device__
    T operator()( T y )
    {
        return y > x;
    }
};

template<typename T>
struct notEqual
{
    const T x;
    notEqual( T _x )
        : x( _x )
    {
    }

    __host__ __device__
    T operator()( const IndexType& value, const IndexType& index )
    {
        if ( value > x )
        {
            return index;
        }
        else
        {
            return T( -1 );
        }
    }
};

template<typename T>
struct greaterThanEqual
{
    const T x;
    greaterThanEqual( T _x )
        : x( _x )
    {
    }
    __host__ __device__
    T operator()( T y )
    {
        return y >= x;
    }
};

template<typename T>
struct identity
{
    __host__ __device__
    double operator()( thrust::tuple<T,T> x )
    {
        return thrust::get < 0 > ( x ) == thrust::get < 1 > ( x );
    }
};

template<typename ValueType,typename OtherValueType>
struct multiply
{
    __host__ __device__
    ValueType operator()( ValueType value, OtherValueType otherValue )
    {
        return value * static_cast<ValueType>( otherValue );
    }
};

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  countNonEmptyRowsBySizes                                          */
/* ------------------------------------------------------------------------------------------------------------------ */

IndexType CUDAELLUtils::countNonEmptyRowsBySizes( const IndexType sizes[], const IndexType numRows )
{
    LAMA_LOG_INFO( logger, "countNonEmptyRowsBySizes #sizes = " << sizes << " #numRows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<IndexType> sizes_ptr( const_cast<IndexType*>( sizes ) );
    IndexType counter = thrust::transform_reduce( sizes_ptr, sizes_ptr + numRows, greaterThan<IndexType>( 0 ), 0,
                        thrust::plus<IndexType>() );
    return counter;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setNonEmptyRowsBySizes                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

void CUDAELLUtils::setNonEmptyRowsBySizes(
    IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const IndexType sizes[],
    const IndexType numRows )
{
    LAMA_LOG_INFO( logger,
                   "setNonEmptyRowsBySizes" << " #rowIndexes = " << rowIndexes << ", #numNonEmptyRows = " << numNonEmptyRows << ", #sizes = " << sizes << ", #numRows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    // Create device ptr and help variables
    thrust::device_ptr<IndexType> rowIndexes_ptr( const_cast<IndexType*>( rowIndexes ) );
    thrust::device_ptr<IndexType> sizes_ptr( const_cast<IndexType*>( sizes ) );
    thrust::counting_iterator<IndexType> sequence( 0 );
    thrust::device_vector<IndexType> tmp( numRows );

    // transform array
    thrust::transform( sizes_ptr, sizes_ptr + numRows, sequence, tmp.begin(), notEqual<IndexType>( 0 ) );
    thrust::copy_if( tmp.begin(), tmp.end(), rowIndexes_ptr, greaterThanEqual<IndexType>( 0 ) );
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  hasDiagonalProperty                                               */
/* ------------------------------------------------------------------------------------------------------------------ */

bool CUDAELLUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType ellJA[] )
{
    LAMA_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<IndexType> ellJA_ptr( const_cast<IndexType*>( ellJA ) );
    thrust::counting_iterator<IndexType> sequence( 0 );

    if ( numDiagonals > 0 )
    {
        bool diagonalProperty = thrust::transform_reduce(
                                    thrust::make_zip_iterator( thrust::make_tuple( ellJA_ptr, sequence ) ),
                                    thrust::make_zip_iterator(
                                        thrust::make_tuple( ellJA_ptr + numDiagonals, sequence + numDiagonals ) ),
                                    identity<IndexType>(), true, thrust::logical_and<bool>() );

        return diagonalProperty;
    }
    else
    {
        return false;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  check                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

__global__
void checkKernel(
    const IndexType mNumRows,
    const IndexType mNumValuesPerRow,
    const IndexType mNumColumns,
    const IndexType *ia,
    const IndexType *ja,
    bool *result )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( i < mNumRows )
    {
        // check ia integrity
        result[i] = ( ia[i] <= mNumValuesPerRow );

        // check ja integrity
        for ( IndexType jj = 0; jj < ia[i]; jj++ )
        {
            IndexType j = ja[jj * mNumRows + i];
            bool jaIntegrity = ( j >= 0 && j < mNumColumns );
            result[i] = result[i] && jaIntegrity;
        }
    }
}

void CUDAELLUtils::check(
    const IndexType mNumRows,
    const IndexType mNumValuesPerRow,
    const IndexType mNumColumns,
    const IndexType *ia,
    const IndexType *ja,
    const char* msg )
{
    LAMA_LOG_INFO( logger,
                   "check # mNumRows = " << mNumRows << ", mNumValuesPerRow = " << mNumValuesPerRow << ", mNumColumns = " << mNumColumns )

    LAMA_CHECK_CUDA_ACCESS

    if ( mNumRows > 0 )
    {
        thrust::device_ptr<bool> resultPtr = thrust::device_malloc<bool>( mNumRows );
        thrust::fill( resultPtr, resultPtr + mNumRows, false );

        bool *resultRawPtr = thrust::raw_pointer_cast( resultPtr );

        const int block_size = 256;
        dim3 dimBlock( block_size, 1, 1 );
        dim3 dimGrid = makeGrid( mNumRows, dimBlock.x );

        checkKernel<<<dimGrid, dimBlock>>>( mNumRows, mNumValuesPerRow, mNumColumns, ia, ja, resultRawPtr );

        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR

        bool integrity = thrust::reduce( resultPtr, resultPtr + mNumRows, true, thrust::logical_and<bool>() );

        LAMA_ASSERT_ERROR( integrity, msg << ": ia to large, or ja out of range" )
    }
    else
    {
        LAMA_ASSERT_ERROR( mNumValuesPerRow == 0,
                           msg << ": mNumValuesPerRow should be 0, but is: " << mNumValuesPerRow )
        LAMA_ASSERT_ERROR( mNumColumns == 0, msg << ": mNumColumns should be 0, but is: " << mNumColumns )
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getRow                                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
__global__
void getRowKernel(
    OtherValueType *row,
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType rowNumColumns,
    const IndexType *ja,
    const ValueType *values )
{
    const int jj = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( jj < rowNumColumns )
    {
        IndexType pos = jj * numRows + i;
        row[ja[pos]] = static_cast<OtherValueType>( values[pos] );
    }
}

template<typename ValueType,typename OtherValueType>
void CUDAELLUtils::getRow(
    OtherValueType *row,
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType *ia,
    const IndexType *ja,
    const ValueType *values )
{
    LAMA_LOG_TRACE( logger, "get row #i = " << i )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<OtherValueType> rowPtr( const_cast<OtherValueType*>( row ) );
    thrust::fill( rowPtr, rowPtr + numColumns, 0.0 );

    thrust::device_ptr<IndexType> iaPtr( const_cast<IndexType*>( ia ) );
    thrust::host_vector<IndexType> rowNumColumns( iaPtr + i, iaPtr + i + 1 );

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( rowNumColumns[0], dimBlock.x );

    //TODO: find better CUDA / Thrust implementation
    getRowKernel<<<dimGrid, dimBlock>>>( row, i, numRows, numColumns, rowNumColumns[0], ja, values );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getValue                                                          */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void getValueKernel(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType rowNumColumns,
    const IndexType *ja,
    const ValueType *values,
    ValueType *result )
{
    const int jj = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( jj < rowNumColumns )
    {
        IndexType pos = jj * numRows + i;
        if ( ja[pos] == j )
        {
            result[jj] = values[pos];
        }
    }
}

template<typename ValueType,typename OtherValueType>
OtherValueType CUDAELLUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType *ia,
    const IndexType *ja,
    const ValueType *values )
{
    LAMA_CHECK_CUDA_ACCESS

    LAMA_LOG_TRACE( logger, "get value i = " << i << ", j = " << j << " numRows = " << numRows )

    thrust::device_ptr<IndexType> iaPtr( const_cast<IndexType*>( ia ) );
    thrust::host_vector<IndexType> rowNumColumnsVec( iaPtr + i, iaPtr + i + 1 );

    IndexType rowNumColumns = rowNumColumnsVec[0];

    if ( rowNumColumns > 0 )
    {
        thrust::device_ptr<ValueType> resultPtr = thrust::device_malloc < ValueType > ( rowNumColumns );
        thrust::fill( resultPtr, resultPtr + rowNumColumns, 0.0 );

        ValueType *resultRawPtr = thrust::raw_pointer_cast( resultPtr );

        const int block_size = 256;
        dim3 dimBlock( block_size, 1, 1 );
        dim3 dimGrid = makeGrid( rowNumColumns, dimBlock.x );

        getValueKernel<<<dimGrid, dimBlock>>>( i, j, numRows, rowNumColumns, ja, values, resultRawPtr );

        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR

        return thrust::reduce( resultPtr, resultPtr + rowNumColumns );
    }
    return 0.0;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  scaleValue                                                        */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void CUDAELLUtils::scaleValue(
    const IndexType numRows,
    const IndexType ia[],
    ValueType mValues[],
    const OtherValueType values[] )
{

    LAMA_LOG_INFO( logger,
                   "scaleValue, #numRows = " << numRows << ", ia = " << ia << ", mValues = " << mValues << ", values = " << values )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<IndexType> ia_ptr( const_cast<IndexType*>( ia ) );
    thrust::device_ptr<ValueType> mValues_ptr( const_cast<ValueType*>( mValues ) );
    thrust::device_ptr<OtherValueType> values_ptr( const_cast<OtherValueType*>( values ) );

    IndexType maxCols = CUDAUtils::maxval( ia, numRows );

    //TODO: maybe find better implementation
    for ( IndexType i = 0; i < maxCols; i++ )
    {
        thrust::transform( mValues_ptr + i * numRows, mValues_ptr + i * numRows + numRows, values_ptr,
                           mValues_ptr + i * numRows, multiply<ValueType,OtherValueType>() );
    }

}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
__global__
void ell2csrKernel(
    IndexType *csrJa,
    ValueType *csrValues,
    const IndexType * const csrIa,
    const IndexType numRows,
    const IndexType * const ellIa,
    const IndexType * const ellJa,
    const OtherValueType * const ellValues )
{

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType rowSize = ellIa[i];
        IndexType offset = csrIa[i];

        for ( IndexType jj = 0; jj < rowSize; ++jj )
        {
            IndexType pos = jj * numRows + i;
            csrJa[offset + jj] = ellJa[pos];
            csrValues[offset + jj] = static_cast<OtherValueType>( ellValues[pos] );
        }
    }

}

template<typename ELLValueType,typename CSRValueType>
void CUDAELLUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ELLValueType ellValues[] )
{
    LAMA_REGION( "CUDA.ELL->CSR_values" )

    LAMA_LOG_INFO( logger,
                   "get CSRValues<" << typeid( ELLValueType ).name() << ", " << typeid( CSRValueType ).name() << ">" << ", #rows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    //TODO: find better CUDA / Thrust implementation
    ell2csrKernel<<<dimGrid, dimBlock>>>( csrJA, csrValues, csrIA, numRows,
                                          ellSizes, ellJA, ellValues);

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T1,typename T2>
__global__
void csr2ellKernel(
    int* ell_ja,
    T1* ell_values,
    const int* const ell_ia,
    int n,
    int ellNumValuesPerRow,
    const int* const csr_ia,
    const int* const csr_ja,
    const T2* const csr_values )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        int ellOffset = i;
        int lastJ = 0;

        for ( int jj = csr_ia[i]; jj < csr_ia[i + 1]; ++jj )
        {
            lastJ = csr_ja[jj];
            ell_ja[ellOffset] = lastJ;
            ell_values[ellOffset] = csr_values[jj];
            ellOffset += n;
        }

        // fill in useful values until length of line

        for ( int jj = ell_ia[i]; jj < ellNumValuesPerRow; ++jj )
        {
            ell_ja[ellOffset] = lastJ;
            ell_values[ellOffset] = 0.0;
            ellOffset += n;
        }
    }
}

template<typename ELLValueType,typename CSRValueType>
void CUDAELLUtils::setCSRValues(
    IndexType ellJA[],
    ELLValueType ellValues[],
    const IndexType ellSizes[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[] )
{
    LAMA_REGION( "CUDA.ELL<-CSR_values" )

    LAMA_LOG_INFO( logger,
                   "set CSRValues<" << typeid( ELLValueType ).name() << ", " << typeid( CSRValueType ).name() << ">" << ", #rows = " << numRows << ", #values/row = " << numValuesPerRow )

    LAMA_LOG_DEBUG( logger,
                    "ellJA = " << ellJA << ", ellValues = " << ellValues << ", ellSizes = " << ellSizes << ", csrIA = " << csrIA << ", csrJA = " << csrJA << ", csrValues = " << csrValues )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    csr2ellKernel<<<dimGrid, dimBlock>>>( ellJA, ellValues, ellSizes, numRows, numValuesPerRow,
                                          csrIA, csrJA, csrValues);
    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR

    // throw exception if ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  SPMV                                                              */
/* ------------------------------------------------------------------------------------------------------------------ */

texture<float,1> texELLSXref;

texture<int2,1> texELLDXref;

template<typename T,bool useTexture>
__inline__     __device__ T fetch_ELLx( const T* const x, const int i )
{
    return x[i];
}

template<>
__inline__ __device__
float fetch_ELLx<float,true>( const float* const, const int i )
{
    return tex1Dfetch( texELLSXref, i );
}

template<>
__inline__ __device__
double fetch_ELLx<double,true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texELLDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<typename T,bool useTexture>
__global__
void ell_agemvpbv_kernel(
    int n,
    T alpha,
    int ellNumValuesPerRow,
    const T* ellValues,
    const int* ellJA,
    const T* const x_d,
    const T beta,
    const T* const z_d,
    T* y_d )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        T summand = 0.0;
        if ( beta != 0.0 )
        {
            summand = beta * z_d[i];
        }

        T value = 0.0;
        int pos = i;
        for ( int kk = 0; kk < ellNumValuesPerRow; ++kk )
        {
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += ellValues[pos] * fetch_ELLx<T,useTexture>( x_d, ellJA[pos] );
            pos += n;
        }
        y_d[i] = alpha * value + summand;
    }

}

template<typename ValueType>
void CUDAELLUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numNonZerosPerRow,
    const IndexType ellIA[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.ELL.normalGEMV" )

    LAMA_LOG_INFO( logger, "normalGEMV<" << typeid(ValueType).name() << ">" << ", #rows = " << numRows )

    LAMA_LOG_INFO( logger,
                   "alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", result = " << result )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    // TODO read nnc ?? TB
    //    if(nnc>134217728)//512MB * 1024 KB/MB * 1024 B/KB / (4 B/float) = 134217728
    //    {
    //        lama_setLastError(LAMA_STATUS_INPUTVECTOR_EXCEEDS_TEXTURESIZE);
    //        return;
    //    }

    // x_d and y_d need to be not aliased
    // TODO: assert( x_d != y_d );

    const int block_size = ( numRows > 8191 ? 256 : 128 );
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    const bool useTexture = CUDATexture::useTexture();

    if ( syncToken )
    {
        //useTexture = false;
    }

    if ( useTexture )
    {

        if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLSXref, x ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }
        else if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLDXref, x ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_agemvpbv_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_agemvpbv_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }

    if ( useTexture )
    {
        ell_agemvpbv_kernel<ValueType, true> <<<dimGrid, dimBlock, 0, stream>>>
        ( numRows, alpha, numNonZerosPerRow, ellValues, ellJA, x, beta, y, result );
    }
    else
    {
        ell_agemvpbv_kernel<ValueType, false> <<<dimGrid, dimBlock, 0, stream>>>
        ( numRows, alpha, numNonZerosPerRow, ellValues, ellJA, x, beta, y, result );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_SELLAGEMVPBV_CUDAKERNEL_FAILED" )

    if ( !syncToken )
    {
        LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), "LAMA_STATUS_SELLAGEMVPBV_CUDAKERNEL_FAILED" )
    }

    if ( useTexture )
    {
        if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texELLSXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" )
        }
        else if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texELLDXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" )
        }
    }
}

template<typename T,bool useTexture>
__global__
void ell_agemvpbsv_kernel(
    int n,
    T alpha,
    int nnr,
    const int* const ia_d,
    const T* ellValues,
    const int* ellJA,
    const int* const rowIndexes,
    const int nzr,
    const T* const x_d,
    T* y_d )
{
    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( id < nzr )
    {
        const int i = rowIndexes[id];

        ellValues += i;
        ellJA += i;
        T value = 0.0;
        const int noneZeros = ellValues[i];
        for ( int kk = 0; kk < noneZeros; ++kk )
        {
            const T aValue = *ellValues;
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += aValue * fetch_ELLx<T,useTexture>( x_d, *ellJA );
            ellValues += n;
            ellJA += n;
        }
        y_d[i] += alpha * value;
    }
}

template<typename ValueType>
void CUDAELLUtils::sparseGEMV(
    ValueType result[],
    const IndexType numRows,
    const IndexType numNonZerosPerRows,
    const ValueType alpha,
    const ValueType x[],
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType ellIA[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.ELL.sparseGEMV" )

    LAMA_LOG_INFO( logger, "sparseGEMV<" << typeid(ValueType).name() << ">" << ", #non-zero rows = " << numNonZeroRows )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

//TODO read nnc
//    if(nnc>134217728)//512MB * 1024 KB/MB * 1024 B/KB / (4 B/float) = 134217728
//    {
//        lama_setLastError(LAMA_STATUS_INPUTVECTOR_EXCEEDS_TEXTURESIZE);
//        return;
//    }
//x_d and y_d need to be not aliased
//TODO: assert(x_d != y_d);
    const int block_size = ( numNonZeroRows > 8191 ? 256 : 128 );
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numNonZeroRows, dimBlock.x );

    const bool useTexture = CUDATexture::useTexture();

    if ( useTexture )
    {
        if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLSXref, x ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }
        else if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLDXref, x ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_agemvpbv_kernel<ValueType, true>, cudaFuncCachePreferL1),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_agemvpbv_kernel<ValueType, false>, cudaFuncCachePreferL1),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }

//    if ( *transa == 'N'|| *transa == 'n')
//    {
    if ( useTexture )
    {
        ell_agemvpbsv_kernel<ValueType, true>
        <<<dimGrid, dimBlock, 0, stream>>>( numRows, alpha, numNonZerosPerRows, ellIA, ellValues,
                                            ellJA, rowIndexes, numNonZeroRows, x, result );
    }
    else
    {
        ell_agemvpbsv_kernel<ValueType, false>
        <<<dimGrid, dimBlock, 0, stream>>>( numRows, alpha, numNonZerosPerRows, ellIA, ellValues,
                                            ellJA, rowIndexes, numNonZeroRows, x, result );
    }
//    }
//    else if ( *transa == 'T' || *transa == 't' ||
//              *transa == 'C' || *transa == 'c')
//    {
//        //TODO: Implement this.
//        lama_setLastError( LAMA_STATUS_NOT_IMPLEMENTED );
//        return;
//    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_SELLAGEMVPBV_CUDAKERNEL_FAILED" )
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), "LAMA_STATUS_SELLAGEMVPBV_CUDAKERNEL_FAILED" )

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaUnbindTexture( texELLSXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" )
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  Jacobi                                                           */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T,bool useTexture>
__global__
void ell_jacobi_kernel(
    const int ellNumValuesPerRow,
    const int* ellJA,
    const T* ellValues,
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
        ellValues += i;
        ellJA += i;

        const T diag = *ellValues;
        ellValues += numRows;
        ellJA += numRows;

        for ( int kk = 1; kk < ellNumValuesPerRow; ++kk )
        {
            const T aValue = *ellValues;
            temp -= aValue * fetch_ELLx<T,useTexture>( oldSolution, *ellJA );
            ellValues += numRows;
            ellJA += numRows;
        }
        if ( omega == 0.5 )
        {
            solution[i] = omega * ( fetch_ELLx<T,useTexture>( oldSolution, i ) + temp / diag );
        }
        else if ( omega == 1.0 )
        {
            solution[i] = temp / diag;
        }
        else
        {
            solution[i] = omega * ( temp / diag ) + ( 1.0 - omega ) * fetch_ELLx<T,useTexture>( oldSolution, i );
        }
    }
}

template<typename ValueType>
void CUDAELLUtils::jacobi(
    ValueType solution[],
    const IndexType numRows,
    const IndexType ellNumValuesPerRow,
    const IndexType* UNUSED( ellSizes ),
    const IndexType ellJA[],
    const ValueType ellValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.ELL.jacobi" )

    LAMA_LOG_INFO( logger, "jacobi, #rows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    const int block_size = ( numRows > 8191 ? 256 : 128 );
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    const bool useTexture = CUDATexture::useTexture();

    if ( useTexture )
    {

        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLDXref, oldSolution ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLSXref, oldSolution ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_jacobi_kernel<ValueType, true>, cudaFuncCachePreferL1),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_jacobi_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }

    if ( useTexture )
    {
        ell_jacobi_kernel<ValueType, true> <<<dimGrid, dimBlock, 0, stream>>>( ellNumValuesPerRow, ellJA, ellValues,
                numRows, rhs, solution, oldSolution, omega );
    }

    else
    {
        ell_jacobi_kernel<ValueType, false> <<<dimGrid, dimBlock, 0, stream>>>( ellNumValuesPerRow, ellJA, ellValues,
                numRows, rhs, solution, oldSolution, omega );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_DCSRJACOBI_CUDAKERNEL_FAILED" )

    if ( useTexture )
    {
        LAMA_CUDA_RT_CALL( cudaUnbindTexture(texELLSXref), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" )
    }

    if ( !syncToken )
    {
        cudaStreamSynchronize( stream );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  Jacobi halo                                                       */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T,bool useTexture>
__global__
void ell_jacobi_halo_kernel(
    T* const solution,
    const T* const diagonal,
    const int* const ellSizes,
    const int* const ellJA,
    const T* const ellValues,
    const int* const rowIndexes,
    const int numNonEmptyRows,
    const int numRows,
    const T* const oldSolution,
    const T omega )
{
    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numNonEmptyRows )
    {
        int i = id;

        if ( rowIndexes )
        {
            i = rowIndexes[id];
        }

        T temp = 0.0;

        int pos = i;
        const int rowEnd = ellSizes[i];
        for ( int jj = 0; jj < rowEnd; ++jj )
        {
            temp += ellValues[pos] * fetch_ELLx<T,useTexture>( oldSolution, ellJA[pos] );
            pos += numRows;
        }

        const T diag = diagonal[i];
        solution[i] -= temp * omega / diag;
    }
}

template<typename ValueType>
void CUDAELLUtils::jacobiHalo(
    ValueType solution[],
    const IndexType numRows,
    const ValueType diagonal[],
    const IndexType UNUSED( ellNumValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const ValueType oldSolution[],
    const ValueType omega,
    SyncToken* syncToken )
{
    LAMA_REGION( "CUDA.ELL.jacobiHalo" )
 
    LAMA_LOG_INFO( logger, "jacobiHalo, #non-empty rows = " << numNonEmptyRows )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = ( numNonEmptyRows > 8191 ? 256 : 128 );
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numNonEmptyRows, dimBlock.x );

    const bool useTexture = CUDATexture::useTexture();

    if ( useTexture )
    {

        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLDXref, oldSolution), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texELLSXref, oldSolution), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }

        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_jacobi_halo_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }
    else
    {
        LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_jacobi_halo_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
    }

    if ( useTexture )
    {
        ell_jacobi_halo_kernel<ValueType, true> <<<dimGrid, dimBlock>>>(
            solution, diagonal, ellSizes, ellJA, ellValues,
            rowIndexes, numNonEmptyRows, numRows, oldSolution, omega );
    }
    else
    {
        ell_jacobi_halo_kernel<ValueType, false> <<<dimGrid, dimBlock>>>(
            solution, diagonal, ellSizes, ellJA, ellValues,
            rowIndexes, numNonEmptyRows, numRows, oldSolution, omega );
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_ELLJACOBIHALO_CUDAKERNEL_FAILED" )
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), "LAMA_STATUS_ELLJACOBIHALO_CUDAKERNEL_FAILED" )

    if ( useTexture )
    {
        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texELLDXref), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" )
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture(texELLSXref), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" )
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                Template instantiations via registration routine                                    */
/* ------------------------------------------------------------------------------------------------------------------ */

void CUDAELLUtils::setInterface( ELLUtilsInterface& ELLUtils )
{
    LAMA_INTERFACE_REGISTER( ELLUtils, countNonEmptyRowsBySizes )
    LAMA_INTERFACE_REGISTER( ELLUtils, setNonEmptyRowsBySizes )
    LAMA_INTERFACE_REGISTER( ELLUtils, hasDiagonalProperty )
    LAMA_INTERFACE_REGISTER( ELLUtils, check )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getRow, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getRow, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getRow, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getRow, double, double )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getValue, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getValue, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getValue, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getValue, double, double )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, scaleValue, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, scaleValue, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, scaleValue, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, scaleValue, double, double )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, setCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, setCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, setCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, setCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, sparseGEMV, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, sparseGEMV, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobi, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobi, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobiHalo, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobiHalo, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDAELLUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.ELLUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDAELLUtils::initialized = registerInterface();

} // namespace lama
