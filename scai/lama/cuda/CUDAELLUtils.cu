/**
 * @file CUDAELLUtils.cu
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
 * @brief Implementation of ELL utilities with CUDA
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/cuda/CUDAELLUtils.hpp>

// local library
#include <scai/lama/cuda/utils.cu.h>
#include <scai/lama/cuda/CUDAUtils.hpp>
#include <scai/lama/cuda/CUDASettings.hpp>

#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>

// internal scai library
#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Constants.hpp>

// CUDA
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

// boost
#include <boost/preprocessor.hpp>

using namespace scai::hmemo;

namespace scai
{

using common::getScalarType;

namespace lama
{

SCAI_LOG_DEF_LOGGER( CUDAELLUtils::logger, "CUDA.ELLUtils" )

/* ------------------------------------------------------------------------------------------------------------------ */

#include <scai/lama/cuda/CUDATexVector.hpp>

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
    double operator()( thrust::tuple<T, T> x )
    {
        return thrust::get < 0 > ( x ) == thrust::get < 1 > ( x );
    }
};

template<typename ValueType, typename OtherValueType>
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
    SCAI_LOG_INFO( logger, "countNonEmptyRowsBySizes #sizes = " << sizes << " #numRows = " << numRows )

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<IndexType> sizes_ptr( const_cast<IndexType*>( sizes ) );
    IndexType counter = thrust::transform_reduce( sizes_ptr, sizes_ptr + numRows,
                        greaterThan<IndexType>( 0 ), 0, thrust::plus<IndexType>() );
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
    SCAI_LOG_INFO( logger,
                   "setNonEmptyRowsBySizes" << " #rowIndexes = " << rowIndexes << ", #numNonEmptyRows = " << numNonEmptyRows << ", #sizes = " << sizes << ", #numRows = " << numRows )

    SCAI_CHECK_CUDA_ACCESS

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
    SCAI_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )

    SCAI_CHECK_CUDA_ACCESS

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
    const IndexType* ia,
    const IndexType* ja,
    bool* result )
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
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType numColumns,
    const IndexType* ia,
    const IndexType* ja,
    const char* msg )
{
    SCAI_LOG_INFO( logger,
                   "check # numRows = " << numRows << ", numValuesPerRow = " << numValuesPerRow << ", numColumns = " << numColumns )

    SCAI_CHECK_CUDA_ACCESS

    if ( numRows > 0 )
    {
        thrust::device_ptr<bool> resultPtr = thrust::device_malloc<bool>( numRows );
        thrust::fill( resultPtr, resultPtr + numRows, false );

        bool* resultRawPtr = thrust::raw_pointer_cast( resultPtr );

        const int blockSize = CUDASettings::getBlockSize( numRows );
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        checkKernel <<< dimGrid, dimBlock>>>( numRows, numValuesPerRow, numColumns, ia, ja, resultRawPtr );

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "fill result with false failed" )

        bool integrity = thrust::reduce( resultPtr, resultPtr + numRows, true, thrust::logical_and<bool>() );

        SCAI_ASSERT_ERROR( integrity, msg << ": ia to large, or ja out of range" )
    }
    else
    {
        SCAI_ASSERT_EQUAL_ERROR( 0, numValuesPerRow )
        SCAI_ASSERT_EQUAL_ERROR( 0, numColumns )
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getRow                                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, typename OtherValueType>
__global__
void getRowKernel(
    OtherValueType* row,
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType rowNumColumns,
    const IndexType* ja,
    const ValueType* values )
{
    const int jj = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( jj < rowNumColumns )
    {
        IndexType pos = jj * numRows + i;
        row[ja[pos]] = static_cast<OtherValueType>( values[pos] );
    }
}

template<typename ValueType, typename OtherValueType>
void CUDAELLUtils::getRow(
    OtherValueType* row,
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType UNUSED( numValuesPerRow ),
    const IndexType* ia,
    const IndexType* ja,
    const ValueType* values )
{
    SCAI_LOG_TRACE( logger, "get row #i = " << i )

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<OtherValueType> rowPtr( const_cast<OtherValueType*>( row ) );
    thrust::fill( rowPtr, rowPtr + numColumns, static_cast<ValueType>( 0.0 ) );

    thrust::device_ptr<IndexType> iaPtr( const_cast<IndexType*>( ia ) );
    thrust::host_vector<IndexType> rowNumColumns( iaPtr + i, iaPtr + i + 1 );

    const int blockSize = CUDASettings::getBlockSize( numRows );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( rowNumColumns[0], dimBlock.x );

    //TODO: find better CUDA / Thrust implementation
    getRowKernel <<< dimGrid, dimBlock>>>( row, i, numRows, numColumns, rowNumColumns[0], ja, values );

    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
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
    const IndexType* ja,
    const ValueType* values,
    ValueType* result )
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

template<typename ValueType>
ValueType CUDAELLUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType UNUSED( numValuesPerRow ),
    const IndexType* ia,
    const IndexType* ja,
    const ValueType* values )
{
    SCAI_CHECK_CUDA_ACCESS

    SCAI_LOG_TRACE( logger, "get value i = " << i << ", j = " << j << " numRows = " << numRows )

    thrust::device_ptr<IndexType> iaPtr( const_cast<IndexType*>( ia ) );
    thrust::host_vector<IndexType> rowNumColumnsVec( iaPtr + i, iaPtr + i + 1 );

    IndexType rowNumColumns = rowNumColumnsVec[0];

    if ( rowNumColumns > 0 )
    {
        thrust::device_ptr<ValueType> resultPtr = thrust::device_malloc < ValueType > ( rowNumColumns );
        thrust::fill( resultPtr, resultPtr + rowNumColumns, static_cast<ValueType>( 0.0 ) );

        ValueType* resultRawPtr = thrust::raw_pointer_cast( resultPtr );

        const int blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( rowNumColumns, dimBlock.x );

        getValueKernel <<< dimGrid, dimBlock>>>( i, j, numRows, rowNumColumns, ja, values, resultRawPtr );

        cudaStreamSynchronize( 0 );
        SCAI_CHECK_CUDA_ERROR

        return thrust::reduce( resultPtr, resultPtr + rowNumColumns );
    }

    return static_cast<ValueType>( 0.0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  scaleValue                                                        */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, typename OtherValueType>
void CUDAELLUtils::scaleValue(
    const IndexType numRows,
    const IndexType UNUSED( numValuesPerRow ),
    const IndexType ia[],
    ValueType ellValues[],
    const OtherValueType values[] )
{

    SCAI_LOG_INFO( logger,
                   "scaleValue, #numRows = " << numRows << ", ia = " << ia << ", ellValues = " << ellValues << ", values = " << values )

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<IndexType> ia_ptr( const_cast<IndexType*>( ia ) );
    thrust::device_ptr<ValueType> ellValues_ptr( const_cast<ValueType*>( ellValues ) );
    thrust::device_ptr<OtherValueType> values_ptr( const_cast<OtherValueType*>( values ) );

    IndexType maxCols = CUDAUtils::maxval( ia, numRows );

    //TODO: maybe find better implementation
    for ( IndexType i = 0; i < maxCols; i++ )
    {
        thrust::transform( ellValues_ptr + i * numRows, ellValues_ptr + i * numRows + numRows, values_ptr,
                           ellValues_ptr + i * numRows, multiply<ValueType, OtherValueType>() );
    }

}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, typename OtherValueType>
__global__
void ell2csrKernel(
    IndexType* csrJa,
    ValueType* csrValues,
    const IndexType* const csrIa,
    const IndexType numRows,
    const IndexType* const ellIa,
    const IndexType* const ellJa,
    const OtherValueType* const ellValues )
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

template<typename ELLValueType, typename CSRValueType>
void CUDAELLUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType UNUSED( numValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ELLValueType ellValues[] )
{
    SCAI_REGION( "CUDA.ELL->CSR_values" )

    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << getScalarType<ELLValueType>() << ", " << getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    //TODO: find better CUDA / Thrust implementation
    ell2csrKernel <<< dimGrid, dimBlock>>>( csrJA, csrValues, csrIA, numRows,
                                            ellSizes, ellJA, ellValues );

    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T1, typename T2>
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

template<typename ELLValueType, typename CSRValueType>
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
    SCAI_REGION( "CUDA.ELL<-CSR_values" )

    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << getScalarType<ELLValueType>() << ", " << getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values/row = " << numValuesPerRow )

    SCAI_LOG_DEBUG( logger,
                    "ellJA = " << ellJA << ", ellValues = " << ellValues << ", ellSizes = " << ellSizes << ", csrIA = " << csrIA << ", csrJA = " << csrJA << ", csrValues = " << csrValues )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    csr2ellKernel <<< dimGrid, dimBlock>>>( ellJA, ellValues, ellSizes, numRows, numValuesPerRow,
                                            csrIA, csrJA, csrValues );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "csr2ellKernel" );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void fillEllKernel(
    int* ell_ja,
    ValueType* ell_values,
    const int* const ell_ia,
    int n,
    int ellNumValuesPerRow )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        int lastJ = 0;

        int pos = ell_ia[i];

        int ellOffset = i + pos * n;

        if ( pos > 0 && pos < ellNumValuesPerRow )
        {
            lastJ = ell_ja[ pos - n ];
        }

        // fill in useful values until length of line

        for ( int jj = pos; jj < ellNumValuesPerRow; ++jj )
        {
            ell_ja[ellOffset] = lastJ;
            ell_values[ellOffset] = static_cast<ValueType>( 0 );
            ellOffset += n;
        }
    }
}

template<typename ValueType>
void CUDAELLUtils::fillELLValues(
    IndexType ellJA[],
    ValueType ellValues[],
    const IndexType ellSizes[],
    const IndexType numRows,
    const IndexType numValuesPerRow )
{
    SCAI_LOG_INFO( logger, "fill ELLValues<" << getScalarType<ValueType>() )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    fillEllKernel <<< dimGrid, dimBlock>>>( ellJA, ellValues, ellSizes, numRows, numValuesPerRow );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "fillEllKernel" );
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*    Kernel for  SMV + SV                                                                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gemv_kernel(
    T* result,
    const T* const x_d,
    const T* const y_d,
    T alpha,
    const T beta,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    const int* ellIA )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T summand = beta * y_d[i];

        T value = 0.0;
        int pos = i;

        for ( int kk = 0; kk < ellIA[i]; ++kk )
        {
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result[i] = alpha * value + summand;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gemv_kernel_alpha_one_beta_one(
    T* result,
    const T* const x_d,
    const T* const y_d,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    const int* ellIA )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T summand = 0.0;
        summand = y_d[i];

        T value = 0.0;
        int pos = i;

        for ( int kk = 0; kk < ellIA[i]; ++kk )
        {
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result[i] = value + summand;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gemv_kernel_alpha_one_beta_zero(
    T* result,
    const T* const x_d,
    const T* const y_d,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    const int* ellIA )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T value = 0.0;
        int pos = i;

        for ( int kk = 0; kk < ellIA[i]; ++kk )
        {
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result[i] = value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void assign_kernel(
    T* result,
    const T* const y_d,
    int numRows )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        result[i] = y_d[i];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gemv_kernel_alpha_one(
    T* result,
    const T* const x_d,
    const T* const y_d,
    const T beta,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    const int* ellIA )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T summand = beta * y_d[i];

        T value = 0.0;
        int pos = i;

        for ( int kk = 0; kk < ellIA[i]; ++kk )
        {
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result[i] = value + summand;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gemv_kernel_beta_one(
    T* result,
    const T* const x_d,
    const T* const y_d,
    T alpha,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    const int* ellIA )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T summand = y_d[i];

        T value = 0.0;
        int pos = i;

        for ( int kk = 0; kk < ellIA[i]; ++kk )
        {
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result[i] = alpha * value + summand;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void scale_kernel(
    T* result,
    const T* const y_d,
    const T beta,
    int numRows )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        result[i] = beta * y_d[i];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gemv_kernel_beta_zero(
    T* result,
    const T* const x_d,
    const T* const y_d,
    T alpha,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    const int* ellIA )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T value = 0.0;
        int pos = i;

        for ( int kk = 0; kk < ellIA[i]; ++kk )
        {
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result[i] = alpha * value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

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
    SCAI_REGION( "CUDA.ELL.normalGEMV" )

    SCAI_LOG_INFO( logger, "normalGEMV<" << getScalarType<ValueType>() << ">" <<
                   " result[ " << numRows << "] = " << alpha << " * A(ell) * x + " << beta << " * y " )

    SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
        SCAI_LOG_INFO( logger, "asyncronous execution on stream " << stream );
    }

    const int blockSize = CUDASettings::getBlockSize();

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    SCAI_LOG_INFO( logger, "Start normal_gemv_kernel<" << getScalarType<ValueType>()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ">>>" )

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_alpha_one_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_alpha_one_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::constants::ZERO && beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( assign_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            assign_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> ( result, y, numRows );
        }
        else if ( alpha == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, beta, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( scale_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            scale_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> ( result, y, beta, numRows );
        }
        else if ( beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, ellIA );
        }
        else if ( beta == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, ellIA );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, beta, ellValues, ellJA, numRows, ellIA );
        }
    }
    else
    {
        if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_alpha_one_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_alpha_one_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::constants::ZERO && beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( assign_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            assign_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> ( result, y, numRows );
        }
        else if ( alpha == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, beta, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( scale_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            scale_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> ( result, y, beta, numRows );
        }
        else if ( beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, ellIA );
        }
        else if ( beta == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, ellIA );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, beta, ellValues, ellJA, numRows, ellIA );
        }
    }

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "sync for normal_gemv_kernel failed" )

        if ( useTexture )
        {
            vectorUnbindTexture( x );
        }
    }
    else
    {
        // synchronization at SyncToken, delay unbind

        if ( useTexture )
        {
            void ( *unbind ) ( const ValueType* ) = &vectorUnbindTexture;

            syncToken->pushRoutine( common::bind( unbind, x ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*    Kernel for  SVM + SV                                                                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gevm_kernel(
    T* result,
    const T* x_d,
    const T* y_d,
    const T alpha,
    const T beta,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    int numColumns,
    const int* ellIA )
{
    // result = alpha * x_d * A + beta * y_d

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        T summand = beta * y_d[i];
        T value = 0.0;

        for ( int j = 0; j < numRows; ++j )
        {
            int pos = j;

            for ( int kk = 0; kk < ellIA[j]; ++kk )
            {
                if ( ellJA[pos] == i )
                {
                    value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, j );
                }

                pos += numRows;
            }
        }

        result[i] = alpha * value + summand;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gevm_kernel_alpha_one_beta_one(
    T* result,
    const T* x_d,
    const T* y_d,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    int numColumns,
    const int* ellIA )
{
    // result = x_d * A + y_d

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        T summand = y_d[i];
        T value = 0.0;

        for ( int j = 0; j < numRows; ++j )
        {
            int pos = j;

            for ( int kk = 0; kk < ellIA[j]; ++kk )
            {
                if ( ellJA[pos] == i )
                {
                    value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, j );
                }

                pos += numRows;
            }
        }

        result[i] = value + summand;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gevm_kernel_alpha_one_beta_zero(
    T* result,
    const T* x_d,
    const T* y_d,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    int numColumns,
    const int* ellIA )
{
    // result = x_d * A

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        T value = 0.0;

        for ( int j = 0; j < numRows; ++j )
        {
            int pos = j;

            for ( int kk = 0; kk < ellIA[j]; ++kk )
            {
                if ( ellJA[pos] == i )
                {
                    value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, j );
                }

                pos += numRows;
            }
        }

        result[i] = value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gevm_kernel_alpha_zero_beta_one(
    T* result,
    const T* x_d,
    const T* y_d,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    int numColumns,
    const int* ellIA )
{
    // result = y_d

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        result[i] = y_d[i];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gevm_kernel_alpha_one(
    T* result,
    const T* x_d,
    const T* y_d,
    const T beta,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    int numColumns,
    const int* ellIA )
{
    // result = x_d * A + beta * y_d

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        T summand = beta * y_d[i];
        T value = 0.0;

        for ( int j = 0; j < numRows; ++j )
        {
            int pos = j;

            for ( int kk = 0; kk < ellIA[j]; ++kk )
            {
                if ( ellJA[pos] == i )
                {
                    value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, j );
                }

                pos += numRows;
            }
        }

        result[i] = value + summand;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gevm_kernel_alpha_zero(
    T* result,
    const T* x_d,
    const T* y_d,
    const T beta,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    int numColumns,
    const int* ellIA )
{
    // result = beta * y_d

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        result[i] = beta * y_d[i];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gevm_kernel_beta_one(
    T* result,
    const T* x_d,
    const T* y_d,
    const T alpha,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    int numColumns,
    const int* ellIA )
{
    // result = alpha * x_d * A + y_d

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        T value = 0.0;

        for ( int j = 0; j < numRows; ++j )
        {
            int pos = j;

            for ( int kk = 0; kk < ellIA[j]; ++kk )
            {
                if ( ellJA[pos] == i )
                {
                    value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, j );
                }

                pos += numRows;
            }
        }

        result[i] = alpha * value + y_d[i];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void normal_gevm_kernel_beta_zero(
    T* result,
    const T* x_d,
    const T* y_d,
    const T alpha,
    const T* ellValues,
    const int* ellJA,
    int numRows,
    int numColumns,
    const int* ellIA )
{
    // result = alpha * x_d * A

    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numColumns )
    {
        T value = 0.0;

        for ( int j = 0; j < numRows; ++j )
        {
            int pos = j;

            for ( int kk = 0; kk < ellIA[j]; ++kk )
            {
                if ( ellJA[pos] == i )
                {
                    value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, j );
                }

                pos += numRows;
            }
        }

        result[i] = alpha * value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDAELLUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    SCAI_LOG_INFO( logger, "normalGEVM<" << getScalarType<ValueType>() << ">" <<
                   " result[ " << numColumns << "] = " << alpha << " * A(ell) * x + " << beta << " * y " )

    SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given

    const int blockSize = CUDASettings::getBlockSize();

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numColumns, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    SCAI_LOG_INFO( logger, "Start normal_gevm_kernel<" << getScalarType<ValueType>()
                   << ", useTexture = " << useTexture << ">" );

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_one_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_one_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( alpha == scai::common::constants::ZERO && beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_zero_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_zero_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( alpha == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, beta, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( alpha == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, beta, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( beta == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, beta, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
    }
    else
    {
        if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_one_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( alpha == scai::common::constants::ONE && beta == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_one_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( alpha == scai::common::constants::ZERO && beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_zero_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_zero_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( alpha == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, beta, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( alpha == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_alpha_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_alpha_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, beta, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( beta == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else if ( beta == scai::common::constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gevm_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

            normal_gevm_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, beta, ellValues, ellJA, numRows, numColumns, ellSizes );
        }
    }

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "normalGEVM, stream = " << stream )
        SCAI_LOG_DEBUG( logger, "normalGEVM<" << getScalarType<ValueType>() << "> synchronized" )
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

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void sparse_gemv_kernel(
    T* const result_d,
    const T* const x_d,
    const T alpha,
    const T* const ellValues,
    const int* const ellIA,
    const int* const ellJA,
    const int* const rowIndexes,
    const int numNonZeroRows,
    int numRows,
    int numValuesPerRow )
{
    // each thread is assigned to one non-zero row

    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numNonZeroRows )
    {
        const int i = rowIndexes[id];

        int pos = i;

        T value = 0.0;

        const int nonZeros = ellIA[i];

        for ( int kk = 0; kk < nonZeros; ++kk )
        {
            const T aValue = ellValues[pos];

            // compute capability >= 2.0: no benefits to mask with value != 0.0

            value += aValue * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result_d[i] += alpha * value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void sparse_gemv_kernel_alpha_one(
    T* const result_d,
    const T* const x_d,
    const T* const ellValues,
    const int* const ellIA,
    const int* const ellJA,
    const int* const rowIndexes,
    const int numNonZeroRows,
    int numRows,
    int numValuesPerRow )
{
    // each thread is assigned to one non-zero row

    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numNonZeroRows )
    {
        const int i = rowIndexes[id];

        int pos = i;

        T value = 0.0;

        const int nonZeros = ellIA[i];

        for ( int kk = 0; kk < nonZeros; ++kk )
        {
            const T aValue = ellValues[pos];

            // compute capability >= 2.0: no benefits to mask with value != 0.0

            value += aValue * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result_d[i] += value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDAELLUtils::sparseGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numRows,
    const IndexType numNonZerosPerRow,
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    SCAI_REGION( "CUDA.ELL.sparseGEMV" )

    SCAI_LOG_INFO( logger, "sparseGEMV<" << getScalarType<ValueType>() << ">" << ", #non-zero rows = " << numNonZeroRows )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    const int blockSize = CUDASettings::getBlockSize( numNonZeroRows );

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numNonZeroRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    if ( useTexture )
    {
        vectorBindTexture( x );
    }

    SCAI_LOG_INFO( logger, "Start ell_sparse_gemv_kernel<" << getScalarType<ValueType>()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ">>>" );

    if ( useTexture )
    {
        if ( alpha == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( sparse_gemv_kernel_alpha_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "cudaFuncSetCacheConfig failed" )

            sparse_gemv_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>>(
                result, x, ellValues, ellSizes, ellJA, rowIndexes, numNonZeroRows, numRows, numNonZerosPerRow );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( sparse_gemv_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "cudaFuncSetCacheConfig failed" )

            sparse_gemv_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>>(
                result, x, alpha, ellValues, ellSizes, ellJA, rowIndexes, numNonZeroRows, numRows, numNonZerosPerRow );
        }
    }
    else
    {
        if ( alpha == scai::common::constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( sparse_gemv_kernel_alpha_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "cudaFuncSetCacheConfig failed" )

            sparse_gemv_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>>(
                result, x, ellValues, ellSizes, ellJA, rowIndexes, numNonZeroRows, numRows, numNonZerosPerRow );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( sparse_gemv_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "cudaFuncSetCacheConfig failed" )

            sparse_gemv_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>>(
                result, x, alpha, ellValues, ellSizes, ellJA, rowIndexes, numNonZeroRows, numRows, numNonZerosPerRow );
        }
    }

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "sparse GEMV kernel failed" )

        if ( useTexture )
        {
            vectorUnbindTexture( x );
        }
    }
    else
    {
        // synchronization at SyncToken, delay unbind

        if ( useTexture )
        {
            void ( *unbind ) ( const ValueType* ) = &vectorUnbindTexture;

            syncToken->pushRoutine( common::bind( unbind, x ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void sparse_gevm_kernel(
    T* const result_d,
    const T* const x_d,
    const T alpha,
    const T* const ellValues,
    const int* const ellSizes,
    const int* const ellJA,
    const int* const rowIndexes,
    const int numNonZeroRows,
    int numRows,
    int numColumns )
{
    // each thread is assigned to one non-zero row

    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numColumns )
    {
        T value = 0.0;

        for ( int i = 0; i < numNonZeroRows; ++i )
        {
            int pos = id;

            const int nonZeros = ellSizes[pos];

            for ( int kk = 0; kk < nonZeros; ++kk )
            {
                if ( ellJA[pos] == id )
                {
                    const T aValue = ellValues[pos];

                    // compute capability >= 2.0: no benefits to mask with value != 0.0

                    value += aValue * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
                }

                pos += numRows;
            }
        }

        result_d[id] += alpha * value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void sparse_gevm_kernel_alpha_one(
    T* const result_d,
    const T* const x_d,
    const T* const ellValues,
    const int* const ellSizes,
    const int* const ellJA,
    const int* const rowIndexes,
    const int numNonZeroRows,
    int numRows,
    int numColumns )
{
    // each thread is assigned to one non-zero row

    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numColumns )
    {
        T value = 0.0;

        for ( int i = 0; i < numNonZeroRows; ++i )
        {
            int pos = id;

            const int nonZeros = ellSizes[pos];

            for ( int kk = 0; kk < nonZeros; ++kk )
            {
                if ( ellJA[pos] == id )
                {
                    const T aValue = ellValues[pos];

                    // compute capability >= 2.0: no benefits to mask with value != 0.0

                    value += aValue * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
                }

                pos += numRows;
            }
        }

        result_d[id] += value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDAELLUtils::sparseGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numNonZerosPerRow,
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    SCAI_LOG_INFO( logger,
                   "sparseGEVM<" << getScalarType<ValueType>() << ">" << ", #non-zero rows = " << numNonZeroRows )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    const int blockSize = CUDASettings::getBlockSize( numNonZeroRows );

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numNonZeroRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    if ( useTexture )
    {
        if ( alpha == scai::common::constants::ONE )
        {
            sparse_gevm_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, ellValues, ellSizes, ellJA, rowIndexes, numNonZeroRows, numRows, numColumns );
        }
        else
        {
            sparse_gevm_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, alpha, ellValues, ellSizes, ellJA, rowIndexes, numNonZeroRows, numRows, numColumns );
        }
    }
    else
    {
        if ( alpha == scai::common::constants::ONE )
        {
            sparse_gevm_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, ellValues, ellSizes, ellJA, rowIndexes, numNonZeroRows, numRows, numColumns );
        }
        else
        {
            sparse_gevm_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream >>>
            ( result, x, alpha, ellValues, ellSizes, ellJA, rowIndexes, numNonZeroRows, numRows, numColumns );
        }
    }

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "sparseGEVM, stream = " << stream )
        SCAI_LOG_INFO( logger, "sparseGEVM<" << getScalarType<ValueType>() << "> synchronized" )
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  Jacobi                                                           */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void ell_jacobi_kernel(
    const int* ellIA,
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

        for ( int kk = 1; kk < ellIA[i]; ++kk )
        {
            const T aValue = *ellValues;
            temp -= aValue * fetchVectorX<T, useTexture>( oldSolution, *ellJA );
            ellValues += numRows;
            ellJA += numRows;
        }

        if ( omega == 0.5 )
        {
            solution[i] = omega * ( fetchVectorX<T, useTexture>( oldSolution, i ) + temp / diag );
        }
        else if ( omega == 1.0 )
        {
            solution[i] = temp / diag;
        }
        else
        {
            solution[i] = omega * ( temp / diag ) + ( 1.0 - omega ) * fetchVectorX<T, useTexture>( oldSolution, i );
        }
    }
}

template<typename ValueType>
void CUDAELLUtils::jacobi(
    ValueType solution[],
    const IndexType numRows,
    const IndexType UNUSED( ellNumValuesPerRow ),
    const IndexType* ellSizes,
    const IndexType ellJA[],
    const ValueType ellValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    SyncToken* syncToken )
{
    SCAI_REGION( "CUDA.ELL.jacobi" )

    SCAI_LOG_INFO( logger, "jacobi, #rows = " << numRows )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    const bool useTexture = CUDASettings::useTexture();

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    const int blockSize = CUDASettings::getBlockSize( numRows );

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    SCAI_LOG_INFO( logger, "Start ell_jacobi_kernel<" << getScalarType<ValueType>()
                   << "> <<< block size = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ">>>" );

    if ( useTexture )
    {
        vectorBindTexture( oldSolution );

        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_jacobi_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        ell_jacobi_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>>( ellSizes, ellJA, ellValues,
                numRows, rhs, solution, oldSolution, omega );
    }

    else
    {
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_jacobi_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        ell_jacobi_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>>( ellSizes, ellJA, ellValues,
                numRows, rhs, solution, oldSolution, omega );
    }

    SCAI_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_DCSRJACOBI_CUDAKERNEL_FAILED" )

    if ( !syncToken )
    {
        // synchronize now and unbind texture if used

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "ELL: jacobiKernel FAILED" )

        if ( useTexture )
        {
            vectorUnbindTexture( oldSolution );
        }
    }
    else
    {
        if ( useTexture )
        {
            void ( *unbind ) ( const ValueType* ) = &vectorUnbindTexture;

            syncToken->pushRoutine( common::bind( unbind, oldSolution ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  Jacobi halo                                                       */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, bool useTexture>
__global__
void ell_jacobi_halo_kernel(
    ValueType* const solution,
    const ValueType* const diagonal,
    const int* const ellSizes,
    const int* const ellJA,
    const ValueType* const ellvalues,
    const int* const rowIndexes,
    const int numnonemptyrows,
    const int numrows,
    const ValueType* const oldsolution,
    const ValueType omega )
{
    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numnonemptyrows )
    {
        int i = id;

        if ( rowIndexes )
        {
            i = rowIndexes[id];
        }

        ValueType temp = 0.0;

        int pos = i;
        const int rowend = ellSizes[i];

        for ( int jj = 0; jj < rowend; ++jj )
        {
            temp += ellvalues[pos] * fetchVectorX<ValueType, useTexture>( oldsolution, ellJA[pos] );
            pos += numrows;
        }

        const ValueType diag = diagonal[i];

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
    SCAI_REGION( "CUDA.ELL.jacobiHalo" )

    SCAI_LOG_INFO( logger, "jacobiHalo, #non-empty rows = " << numNonEmptyRows )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize( numNonEmptyRows );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numNonEmptyRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    if ( useTexture )
    {
        vectorBindTexture( oldSolution );

        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_jacobi_halo_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        ell_jacobi_halo_kernel<ValueType, true> <<< dimGrid, dimBlock>>>(
            solution, diagonal, ellSizes, ellJA, ellValues,
            rowIndexes, numNonEmptyRows, numRows, oldSolution, omega );
    }
    else
    {
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( ell_jacobi_halo_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )

        ell_jacobi_halo_kernel<ValueType, false> <<< dimGrid, dimBlock>>>(
            solution, diagonal, ellSizes, ellJA, ellValues,
            rowIndexes, numNonEmptyRows, numRows, oldSolution, omega );
    }

    SCAI_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_ELLJACOBIHALO_CUDAKERNEL_FAILED" )
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "LAMA_STATUS_ELLJACOBIHALO_CUDAKERNEL_FAILED" )

    if ( useTexture )
    {
        vectorUnbindTexture( oldSolution );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                Template instantiations via registration routine                                    */
/* ------------------------------------------------------------------------------------------------------------------ */

void CUDAELLUtils::registerKernels()
{
    SCAI_LOG_INFO( logger, "register ELL routines for CUDA in KernelRegistry" )

    using kregistry::KernelRegistry;

    // ctx will contain the context for which registration is done, here Host

    common::ContextType ctx = common::context::CUDA;

    KernelRegistry::set<ELLUtilsInterface::countNonEmptyRowsBySizes>( countNonEmptyRowsBySizes, ctx );
    KernelRegistry::set<ELLUtilsInterface::setNonEmptyRowsBySizes>( setNonEmptyRowsBySizes, ctx );
    KernelRegistry::set<ELLUtilsInterface::hasDiagonalProperty>( hasDiagonalProperty, ctx );
    KernelRegistry::set<ELLUtilsInterface::check>( check, ctx );

#define LAMA_ELL_UTILS2_REGISTER(z, J, TYPE )                                             \
    KernelRegistry::set<ELLUtilsInterface::getRow<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getRow, ctx );       \
    KernelRegistry::set<ELLUtilsInterface::scaleValue<TYPE, ARITHMETIC_HOST_TYPE_##J> >( scaleValue, ctx );       \
    KernelRegistry::set<ELLUtilsInterface::setCSRValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setCSRValues, ctx );       \
    KernelRegistry::set<ELLUtilsInterface::getCSRValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getCSRValues, ctx );       \
     
#define LAMA_ELL_UTILS_REGISTER(z, I, _)                                                  \
    KernelRegistry::set<ELLUtilsInterface::normalGEMV<ARITHMETIC_HOST_TYPE_##I> >( normalGEMV, ctx );       \
    KernelRegistry::set<ELLUtilsInterface::normalGEVM<ARITHMETIC_HOST_TYPE_##I> >( normalGEVM, ctx );       \
    KernelRegistry::set<ELLUtilsInterface::sparseGEMV<ARITHMETIC_HOST_TYPE_##I> >( sparseGEMV, ctx );       \
    KernelRegistry::set<ELLUtilsInterface::sparseGEVM<ARITHMETIC_HOST_TYPE_##I> >( sparseGEVM, ctx );       \
    KernelRegistry::set<ELLUtilsInterface::jacobi<ARITHMETIC_HOST_TYPE_##I> >( jacobi, ctx );               \
    KernelRegistry::set<ELLUtilsInterface::jacobiHalo<ARITHMETIC_HOST_TYPE_##I> >( jacobiHalo, ctx );       \
    KernelRegistry::set<ELLUtilsInterface::getValue<ARITHMETIC_HOST_TYPE_##I> >( getValue, ctx );           \
    KernelRegistry::set<ELLUtilsInterface::fillELLValues<ARITHMETIC_HOST_TYPE_##I> >( fillELLValues, ctx );  \
                                                                                                             \
    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT,                                            \
                     LAMA_ELL_UTILS2_REGISTER,                                            \
                     ARITHMETIC_CUDA_TYPE_##I )                                           \
     
    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT, LAMA_ELL_UTILS_REGISTER, _ )

#undef LAMA_ELL_UTILS_REGISTER
#undef LAMA_ELL_UTILS2_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDAELLUtils::registerInterface()
{
    registerKernels();
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDAELLUtils::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
