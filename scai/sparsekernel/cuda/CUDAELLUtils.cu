/**
 * @file CUDAELLUtils.cu
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
 * @brief Implementation of ELL utilities with CUDA
 * @author Thomas Brandes
 * @date 04.07.2012
 */

// hpp
#include <scai/sparsekernel/cuda/CUDAELLUtils.hpp>

// local library
#include <scai/sparsekernel/ELLKernelTrait.hpp>

// internal scai library
#include <scai/utilskernel/cuda/CUDAUtils.hpp>
#include <scai/utilskernel/cuda/CUDAReduceUtils.hpp>
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAUtils.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>

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

#include <functional>

namespace scai
{

using common::TypeTraits;
using common::CUDASettings;
using tasking::CUDAStreamSyncToken;
using utilskernel::CUDAUtils;
using utilskernel::CUDAReduceUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( CUDAELLUtils::logger, "CUDA.ELLUtils" )

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
struct changeIndexWithZeroSize
{
    __host__ __device__
    T operator()( const IndexType& value, const IndexType& index )
    {
        if ( value > 0 )
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
struct isOkay
{
    __host__ __device__
    bool operator()( T y )
    {
        return y != T( - 1 );
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
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < mNumRows )
    {
        // check ia integrity
        result[i] = ( ia[i] <= mNumValuesPerRow );

        // check ja integrity
        for ( IndexType jj = 0; jj < ia[i]; jj++ )
        {
            IndexType j = ja[jj * mNumRows + i];
            bool jaIntegrity = common::Utils::validIndex( j, mNumColumns );
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
        SCAI_ASSERT_EQ_ERROR( 0, numValuesPerRow, "as numRows == 0" )
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getRow                                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void getRowKernel(
    ValueType* row,
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType rowNumColumns,
    const IndexType* ja,
    const ValueType* values )
{
    const IndexType jj = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( jj < rowNumColumns )
    {
        IndexType pos = jj * numRows + i;
        row[ja[pos]] = values[pos];
    }
}

template<typename ValueType>
void CUDAELLUtils::getRow(
    ValueType* row,
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType SCAI_UNUSED( numValuesPerRow ),
    const IndexType* ia,
    const IndexType* ja,
    const ValueType* values )
{
    SCAI_LOG_TRACE( logger, "get row #i = " << i )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> rowPtr( const_cast<ValueType*>( row ) );
    thrust::fill( rowPtr, rowPtr + numColumns, static_cast<ValueType>( 0.0 ) );
    thrust::device_ptr<IndexType> iaPtr( const_cast<IndexType*>( ia ) );
    thrust::host_vector<IndexType> rowNumColumns( iaPtr + i, iaPtr + i + 1 );
    const int blockSize = CUDASettings::getBlockSize( numRows );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( rowNumColumns[0], dimBlock.x );
    //TODO: find better CUDA / Thrust implementation
    getRowKernel <<< dimGrid, dimBlock>>>( row, i, numRows, numColumns, rowNumColumns[0], ja, values );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "getRowKernel failed" ) ;
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
    const IndexType jj = threadId( gridDim, blockIdx, blockDim, threadIdx );

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
    const IndexType SCAI_UNUSED( numValuesPerRow ),
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
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "getValueKernel failed" );
        SCAI_CHECK_CUDA_ERROR
        return thrust::reduce( resultPtr, resultPtr + rowNumColumns );
    }

    return static_cast<ValueType>( 0.0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  scaleValue                                                        */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDAELLUtils::scaleRows(
    ValueType ellValues[],
    const IndexType numRows,
    const IndexType SCAI_UNUSED( numValuesPerRow ),
    const IndexType ia[],
    const ValueType values[] )
{
    SCAI_LOG_INFO( logger,
                   "scaleRows, #numRows = " << numRows << ", ia = " << ia << ", ellValues = " << ellValues << ", values = " << values )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<IndexType> ia_ptr( const_cast<IndexType*>( ia ) );
    thrust::device_ptr<ValueType> ellValues_ptr( const_cast<ValueType*>( ellValues ) );
    thrust::device_ptr<ValueType> values_ptr( const_cast<ValueType*>( values ) );

    IndexType maxCols = CUDAReduceUtils::reduce( ia, numRows, IndexType( 0 ), common::BinaryOp::MAX );

    // TODO: maybe find better implementation

    for ( IndexType i = 0; i < maxCols; i++ )
    {
        thrust::transform( ellValues_ptr + i * numRows, ellValues_ptr + i * numRows + numRows, values_ptr,
                           ellValues_ptr + i * numRows, multiply<ValueType, ValueType>() );
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
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

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
    const IndexType SCAI_UNUSED( numValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ELLValueType ellValues[] )
{
    SCAI_REGION( "CUDA.ELL->CSR_values" )
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<ELLValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" << ", #rows = " << numRows )

    if ( numRows == 0 )
    {
        return;   // do not call the kernel as launch configuration params will be invalid
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    //TODO: find better CUDA / Thrust implementation
    ell2csrKernel <<< dimGrid, dimBlock>>>( csrJA, csrValues, csrIA, numRows,
                                            ellSizes, ellJA, ellValues );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "ell2csrKernel failed" );
    SCAI_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T1, typename T2>
__global__
void csr2ellKernel(
    IndexType* ell_ja,
    T1* ell_values,
    const IndexType* const ell_ia,
    IndexType n,
    IndexType ellNumValuesPerRow,
    const IndexType* const csr_ia,
    const IndexType* const csr_ja,
    const T2* const csr_values )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        IndexType ellOffset = i;
        IndexType lastJ = 0;

        for ( IndexType jj = csr_ia[i]; jj < csr_ia[i + 1]; ++jj )
        {
            lastJ = csr_ja[jj];
            ell_ja[ellOffset] = lastJ;
            ell_values[ellOffset] = csr_values[jj];
            ellOffset += n;
        }

        // fill in useful values until length of line

        for ( IndexType jj = ell_ia[i]; jj < ellNumValuesPerRow; ++jj )
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
    SCAI_REGION( "CUDA.ELL.setCSR" )
    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << TypeTraits<ELLValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">"
                   << ", #rows = " << numRows << ", #values/row = " << numValuesPerRow )

    SCAI_LOG_DEBUG( logger,
                    "ellJA = " << ellJA << ", ellValues = " << ellValues << ", ellSizes = " << ellSizes
                    << ", csrIA = " << csrIA << ", csrJA = " << csrJA << ", csrValues = " << csrValues )

    if ( numRows == 0 )
    {
        return;   // do not call the kernel as launch configuration params will be invalid
    }

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
    IndexType* ell_ja,
    ValueType* ell_values,
    const IndexType* const ell_ia,
    IndexType n,
    IndexType ellNumValuesPerRow )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        IndexType lastJ = 0;
        IndexType nRow = ell_ia[i];
        IndexType ellOffset = i + nRow * n;

        if ( nRow > 0 && nRow < ellNumValuesPerRow )
        {
            lastJ = ell_ja[ ellOffset - n ];
        }

        // fill in useful values until length of line

        for ( IndexType jj = nRow; jj < ellNumValuesPerRow; ++jj )
        {
            // ell_ja[ellOffset] = lastJ;
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
    SCAI_LOG_INFO( logger, "fill ELLValues<" << TypeTraits<ValueType>::id() )

    if ( numRows == 0 )
    {
        return;   // do not call the kernel as launch configuration params will be invalid
    }

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
    const IndexType* ellJA,
    IndexType numRows,
    const IndexType* ellIA )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T summand = beta * y_d[i];
        T value = 0.0;
        IndexType pos = i;

        for ( IndexType kk = 0; kk < ellIA[i]; ++kk )
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
    const IndexType* ellJA,
    IndexType numRows,
    const IndexType* ellIA )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T summand = 0.0;
        summand = y_d[i];
        T value = 0.0;
        IndexType pos = i;

        for ( IndexType kk = 0; kk < ellIA[i]; ++kk )
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
    const IndexType* ellJA,
    IndexType numRows,
    const IndexType* ellIA )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T value = 0.0;
        IndexType pos = i;

        for ( IndexType kk = 0; kk < ellIA[i]; ++kk )
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
    IndexType numRows )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

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
    const IndexType* ellJA,
    IndexType numRows,
    const IndexType* ellIA )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T summand = beta * y_d[i];
        T value = 0.0;
        IndexType pos = i;

        for ( IndexType kk = 0; kk < ellIA[i]; ++kk )
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
    const IndexType* ellJA,
    IndexType numRows,
    const IndexType* ellIA )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T summand = y_d[i];
        T value = 0.0;
        IndexType pos = i;

        for ( IndexType kk = 0; kk < ellIA[i]; ++kk )
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
    IndexType numRows )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

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
    const IndexType* ellJA,
    IndexType numRows,
    const IndexType* ellIA )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T value = 0.0;
        IndexType pos = i;

        for ( IndexType kk = 0; kk < ellIA[i]; ++kk )
        {
            //if (aValue != 0.0) //compute capability >= 2.0  => disadvantage
            value += ellValues[pos] * fetchVectorX<T, useTexture>( x_d, ellJA[pos] );
            pos += numRows;
        }

        result[i] = alpha * value;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*    Kernel for  result += alpha * transpose( matrix_ell ) * vector_x                                                */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void normal_gevm_kernel(
    ValueType result[],
    const ValueType x[],
    const ValueType alpha,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    IndexType numRows )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        ValueType xi  = x[i];
        IndexType pos = i;

        for ( IndexType kk = 0; kk < ellSizes[i]; ++kk )
        {
            IndexType j = ellJA[pos];
            ValueType v = alpha * ellValues[pos] * xi;
            common::CUDAUtils::atomicAdd( &result[j], v );
            pos += numRows;
        }
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
    const IndexType numColumns,
    const IndexType numNonZerosPerRow,
    const IndexType ellIA[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const common::MatrixOp op )
{
    SCAI_REGION( common::isTranspose( op ) ? "CUDA.ELL.gemv_t" : "CUDA_ELL.gemv_n" )

    const IndexType nTarget = common::isTranspose( op ) ? numColumns : numRows;

    SCAI_LOG_INFO( logger, "normalGEMV<" << TypeTraits<ValueType>::id() << ">" <<
                   " result[ " << nTarget << "] = " << alpha << " * A(ell) * x + " << beta << " * y " )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
        SCAI_LOG_INFO( logger, "asyncronous execution on stream " << stream );
    }

    const int blockSize = CUDASettings::getBlockSize();

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    SCAI_LOG_INFO( logger, "Start normal_gemv_kernel<" << TypeTraits<ValueType>::id()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ">>>" )

    // set result = beta * y, not needed if beta == 1 and y == result

    if ( common::isTranspose( op ) )
    {
        useTexture = false;    // does not help here 

        CUDAUtils::binaryOpScalar( result, y, beta, numColumns, common::BinaryOp::MULT, false );

        SCAI_LOG_DEBUG( logger, "Launch normal_gevm_kernel<" << TypeTraits<ValueType>::id() << ">" );
    
        normal_gevm_kernel<ValueType> <<< dimGrid, dimBlock, 0, stream>>> (
            result, x, alpha, ellIA, ellJA, ellValues, numRows );
    }
    else if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == scai::common::Constants::ONE && beta == scai::common::Constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::Constants::ONE && beta == scai::common::Constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_zero<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_zero<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::Constants::ZERO && beta == scai::common::Constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( assign_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            assign_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> ( result, y, numRows );
        }
        else if ( alpha == scai::common::Constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, beta, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::Constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( scale_kernel<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            scale_kernel<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> ( result, y, beta, numRows );
        }
        else if ( beta == scai::common::Constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_one<ValueType, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_beta_one<ValueType, true> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, ellIA );
        }
        else if ( beta == scai::common::Constants::ZERO )
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
        if ( alpha == scai::common::Constants::ONE && beta == scai::common::Constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::Constants::ONE && beta == scai::common::Constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one_beta_zero<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one_beta_zero<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::Constants::ZERO && beta == scai::common::Constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( assign_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            assign_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> ( result, y, numRows );
        }
        else if ( alpha == scai::common::Constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_alpha_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, beta, ellValues, ellJA, numRows, ellIA );
        }
        else if ( alpha == scai::common::Constants::ZERO )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( scale_kernel<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            scale_kernel<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> ( result, y, beta, numRows );
        }
        else if ( beta == scai::common::Constants::ONE )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( normal_gemv_kernel_beta_one<ValueType, false>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" )
            normal_gemv_kernel_beta_one<ValueType, false> <<< dimGrid, dimBlock, 0, stream>>> (
                result, x, y, alpha, ellValues, ellJA, numRows, ellIA );
        }
        else if ( beta == scai::common::Constants::ZERO )
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
            syncToken->pushRoutine( std::bind( unbind, x ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void sparse_gevm_kernel(
    ValueType result[],
    const ValueType x[],
    const ValueType alpha,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const IndexType numRows,
    const IndexType rowIndexes[],
    IndexType numNonZeroRows )
{
    const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( ii < numNonZeroRows )
    {
        IndexType i   = rowIndexes[ii];
        ValueType xi  = x[i];
        IndexType pos = i;

        for ( IndexType kk = 0; kk < ellSizes[i]; ++kk )
        {
            IndexType j = ellJA[pos];
            ValueType v = alpha * ellValues[pos] * xi;
            common::CUDAUtils::atomicAdd( &result[j], v );
            pos += numRows;
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
    const IndexType* const ellIA,
    const IndexType* const ellJA,
    const IndexType* const rowIndexes,
    const IndexType numNonZeroRows,
    IndexType numRows,
    IndexType numValuesPerRow )
{
    // each thread is assigned to one non-zero row
    const IndexType id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numNonZeroRows )
    {
        const IndexType i = rowIndexes[id];
        IndexType pos = i;
        T value = 0.0;
        const IndexType nonZeros = ellIA[i];

        for ( IndexType kk = 0; kk < nonZeros; ++kk )
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
    const IndexType* const ellIA,
    const IndexType* const ellJA,
    const IndexType* const rowIndexes,
    const IndexType numNonZeroRows,
    IndexType numRows,
    IndexType numValuesPerRow )
{
    // each thread is assigned to one non-zero row
    const IndexType id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numNonZeroRows )
    {
        const IndexType i = rowIndexes[id];
        IndexType pos = i;
        T value = 0.0;
        const IndexType nonZeros = ellIA[i];

        for ( IndexType kk = 0; kk < nonZeros; ++kk )
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
    const common::MatrixOp op )
{
    SCAI_REGION( "CUDA.ELL.sparseGEMV" )
    SCAI_LOG_INFO( logger, "sparseGEMV<" << TypeTraits<ValueType>::id() << ">" << ", #non-zero rows = " << numNonZeroRows )
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    const IndexType blockSize = CUDASettings::getBlockSize( numNonZeroRows );

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numNonZeroRows, dimBlock.x );

    bool useTexture = CUDASettings::useTexture();

    SCAI_LOG_INFO( logger, "Start ell_sparse_gemv_kernel<" << TypeTraits<ValueType>::id()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ">>>" );

    if ( common::isTranspose( op ) )
    {
        useTexture = false;

        sparse_gevm_kernel<ValueType> <<< dimGrid, dimBlock, 0, stream >>>
        ( result, x, alpha, ellSizes, ellJA, ellValues, numRows, rowIndexes, numNonZeroRows );
    }
    else if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == scai::common::Constants::ONE )
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
        if ( alpha == scai::common::Constants::ONE )
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
            syncToken->pushRoutine( std::bind( unbind, x ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  Jacobi                                                           */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T, bool useTexture>
__global__
void ell_jacobi_kernel(
    const IndexType* ellIA,
    const IndexType* ellJA,
    const T* ellValues,
    const IndexType numRows,
    const T* const rhs,
    T* const solution,
    const T* const oldSolution,
    const T omega )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        T temp = rhs[i];
        ellValues += i;
        ellJA += i;
        const T diag = *ellValues;
        ellValues += numRows;
        ellJA += numRows;

        for ( IndexType kk = 1; kk < ellIA[i]; ++kk )
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
    const IndexType SCAI_UNUSED( ellNumValuesPerRow ),
    const IndexType* ellSizes,
    const IndexType ellJA[],
    const ValueType ellValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega )
{
    SCAI_REGION( "CUDA.ELL.jacobi" )
    SCAI_LOG_INFO( logger, "jacobi, #rows = " << numRows )
    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;
    const bool useTexture = CUDASettings::useTexture();
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    const IndexType blockSize = CUDASettings::getBlockSize( numRows );

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    SCAI_LOG_INFO( logger, "Start ell_jacobi_kernel<" << TypeTraits<ValueType>::id()
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
            syncToken->pushRoutine( std::bind( unbind, oldSolution ) );
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
    const IndexType* const ellSizes,
    const IndexType* const ellJA,
    const ValueType* const ellvalues,
    const IndexType* const rowIndexes,
    const IndexType numnonemptyrows,
    const IndexType numrows,
    const ValueType* const oldsolution,
    const ValueType omega )
{
    const IndexType id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( id < numnonemptyrows )
    {
        IndexType i = id;

        if ( rowIndexes )
        {
            i = rowIndexes[id];
        }

        ValueType temp = 0.0;
        IndexType pos = i;
        const IndexType rowend = ellSizes[i];

        for ( IndexType jj = 0; jj < rowend; ++jj )
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
    const IndexType SCAI_UNUSED( ellNumValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const ValueType oldSolution[],
    const ValueType omega )
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

template<typename ValueType>
__global__
void ell_compressIA_kernel(
    const IndexType* IA,
    const IndexType* JA,
    const ValueType* ellValues,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const RealType<ValueType> eps,
    bool keepDiagonal,
    IndexType* newIA )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType length = IA[i];

        for ( IndexType jj = 0; jj < IA[i]; jj++ )
        {
            IndexType pos = jj * numRows + i;

            if ( keepDiagonal && JA[pos] == i )
            {
                continue;
            }

            if ( common::Math::abs( ellValues[pos] ) > eps )
            {
                continue;
            }

            length--;
        }

        newIA[i] = length;
    }
}

template<typename ValueType>
void CUDAELLUtils::compressIA(
    IndexType newIA[],
    const IndexType IA[],
    const IndexType JA[],
    const ValueType ellValues[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const RealType<ValueType> eps,
    const bool keepDiagonal )
{
    SCAI_LOG_INFO( logger, "compressIA with eps = " << eps )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    ell_compressIA_kernel <<< dimGrid, dimBlock>>>( IA, JA, ellValues, numRows, numValuesPerRow, eps, keepDiagonal, newIA );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "compress" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void ell_compressValues_kernel(
    const IndexType* IA,
    const IndexType* JA,
    const ValueType* values,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const RealType<ValueType> eps,
    bool keepDiagonal,
    const IndexType newNumValuesPerRow,
    IndexType* newJA,
    ValueType* newValues )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType gap = 0;

        for ( IndexType jj = 0; jj < IA[i]; jj++ )
        {
            IndexType pos = jj * numRows + i;

            // delete it if zero and not diagonal entry

            if ( common::Math::abs( values[pos] ) > eps || ( keepDiagonal && JA[pos] == i ) )
            {   
                // move entry gap positions back in this row
                
                IndexType newpos = ( jj - gap ) * numRows + i; 
                newValues[newpos] = values[pos];
                newJA[newpos] = JA[pos];
            }
            else
            {   
                gap++;
            }
        }

        // fill up to top

        for (  IndexType jj = IA[i] - gap; jj < newNumValuesPerRow; jj++ )
        {
            IndexType newpos = jj * numRows + i; 
            newValues[newpos] = 0;
            newJA[newpos] = 0;
        }
    }
}

template<typename ValueType>
void CUDAELLUtils::compressValues(
    IndexType newJA[],
    ValueType newValues[],
    const IndexType newNumValuesPerRow,
    const IndexType IA[],
    const IndexType JA[],
    const ValueType values[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const RealType<ValueType> eps,
    const bool keepDiagonal )
{
    SCAI_LOG_INFO( logger, "compressValues ( #rows = " << numRows
                   << ", values per row (old/new) = " << numValuesPerRow << " / " << newNumValuesPerRow
                   << ") with eps = " << eps )

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    ell_compressValues_kernel <<< dimGrid, dimBlock>>>( IA, JA, values, numRows, numValuesPerRow, eps, keepDiagonal,
            newNumValuesPerRow, newJA, newValues );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "compress" )
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                Template instantiations via registration routine                                    */
/* ------------------------------------------------------------------------------------------------------------------ */

void CUDAELLUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_INFO( logger, "register ELLUtils CUDA-routines for CUDA at kernel registry [" << flag << "]" )
    KernelRegistry::set<ELLKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::check>( check, ctx, flag );
    // KernelRegistry::set<ELLKernelTrait::getValuePos >( getValuePos, ctx, flag );
}

template<typename ValueType>
void CUDAELLUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register ELLUtils CUDA-routines for CUDA at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<ELLKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getRow<ValueType> >( getRow, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::fillELLValues<ValueType> >( fillELLValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::compressIA<ValueType> >( compressIA, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::compressValues<ValueType> >( compressValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::scaleRows<ValueType> >( scaleRows, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void CUDAELLUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register ELLUtils CUDA-routines for CUDA at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    KernelRegistry::set<ELLKernelTrait::setCSRValues<ValueType, OtherValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDAELLUtils::CUDAELLUtils()
{
    SCAI_LOG_INFO( logger, "register ELLUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_CUDA_LIST, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAELLUtils::~CUDAELLUtils()
{
    SCAI_LOG_INFO( logger, "unregister ELLUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_CUDA_LIST, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAELLUtils CUDAELLUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
