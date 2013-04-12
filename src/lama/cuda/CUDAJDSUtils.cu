/**
 * @file CUDAJDSUtils.cpp
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
 * @brief Implementation of JDS utilities with CUDA
 * @author Bea Hornef, Thomas Brandes
 * @date 04.07.2012
 * $Id$
 */

// hpp
#include <lama/cuda/utils.cu.h>

// others
#include <lama/cuda/CUDAStreamSyncToken.hpp>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAJDSUtils.hpp>
#include <lama/cuda/CUDAUtils.hpp>

#include <lama/exception/LAMAAssert.hpp>

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// macros
#include <lama/macros/unused.hpp>

// thrust
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
//#include <thrust/gather.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/scatter.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDAJDSUtils::logger, "CUDA.JDSUtils" )

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  thrust functors                                                   */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T>
struct identity
{
    const T x;
    identity( T _x )
        : x( _x )
    {
    }
    __host__ __device__
    T operator()( thrust::tuple<T,T> y )
    {
        if ( thrust::get < 0 > ( y ) == x )
        {
            return thrust::get < 1 > ( y );
        }
        return 0;
    }
};

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

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setDiagonalWithScalar                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDAJDSUtils::setDiagonalWithScalar( const IndexType numDiagonal, ValueType values[], Scalar scalar )
{
    LAMA_LOG_INFO( logger, "setDiagonalWithScalar with numDiagonal = " << numDiagonal << " and scalar = " << scalar )

    LAMA_CHECK_CUDA_ACCESS

    ValueType value = scalar.getValue<ValueType>();

    thrust::device_ptr<ValueType> valuesPtr( const_cast<ValueType*>( values ) );
    thrust::fill( valuesPtr, valuesPtr + numDiagonal, value );
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getRow                                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
__global__
void getRowKernel(
    OtherValueType *row,
    const IndexType i,
    const IndexType *ilg,
    const IndexType *dlg,
    const IndexType *ja,
    const ValueType *values )
{
    IndexType offset = 0;

    for ( IndexType j = 0; j < ilg[i]; j++ )
    {
        row[ja[i + offset]] = static_cast<OtherValueType>( values[i + offset] );
        offset += dlg[j];
    }
}

template<typename ValueType,typename OtherValueType>
void CUDAJDSUtils::getRow(
    OtherValueType row[],
    const IndexType i,
    const IndexType numColumns,
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType ja[],
    const ValueType values[] )
{
    LAMA_LOG_INFO( logger, "getRow with i = " << i << ", numColumns = " << numColumns << " and numRows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<OtherValueType> rowPtr( const_cast<OtherValueType*>( row ) );
    thrust::device_ptr<IndexType> permPtr( const_cast<IndexType*>( perm ) );

    thrust::fill( rowPtr, rowPtr + numColumns, 0 );

    thrust::counting_iterator<IndexType> sequence( 0 );

    // correct index with permutation array
    IndexType ii = thrust::transform_reduce(
                       thrust::make_zip_iterator( thrust::make_tuple( permPtr, sequence ) ),
                       thrust::make_zip_iterator( thrust::make_tuple( permPtr + numRows, sequence + numRows ) ),
                       identity<IndexType>( i ), 0, thrust::plus<IndexType>() );

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( 1, dimBlock.x );

    //TODO: find better CUDA / Thrust implementation
    getRowKernel<<<dimGrid, dimBlock>>>( row, ii, ilg, dlg, ja, values );

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
    const IndexType* dlg,
    const IndexType* ilg,
    const IndexType* perm,
    const IndexType* ja,
    const ValueType* values,
    ValueType* result )
{
    const int tId = threadId( gridDim, blockIdx, blockDim, threadIdx );
    if ( tId == 0 )
    {
        IndexType ii;

        // check the permutation of row i
        for ( ii = 0; ii < numRows; ii++ )
        {
            if ( perm[ii] == i )
            {
                break;
            }
        }

        IndexType k = 0;
        bool found = false;

        for ( IndexType jj = 0; jj < ilg[ii]; jj++ )
        {
            if ( ja[ii + k] == j )
            {
                result[0] = values[ii + k];
                found = true;
                break;
            }

            k += dlg[jj];
        }

        if ( !found )
        {
            result[0] = 0.0;
        }
    }
}

template<typename ValueType,typename NoType>
ValueType CUDAJDSUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType* dlg,
    const IndexType* ilg,
    const IndexType* perm,
    const IndexType* ja,
    const ValueType* values )
{
    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> resultPtr = thrust::device_malloc < ValueType > ( 1 );
    ValueType *resultRawPtr = thrust::raw_pointer_cast( resultPtr );

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( 1, dimBlock.x );

    //TODO: find better CUDA / Thrust implementation
    getValueKernel<<<dimGrid, dimBlock>>>( i, j, numRows, dlg, ilg, perm, ja, values, resultRawPtr );

    thrust::host_vector<ValueType> resultHost( resultPtr, resultPtr + 1 );

    return resultHost[0];
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  scaleValue                                                        */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
__global__
void scaleValueKernel(
    const IndexType numRows,
    const IndexType *perm,
    const IndexType *ilg,
    const IndexType *dlg,
    ValueType *mValues,
    const OtherValueType *values )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType offset = i;
        OtherValueType value = values[perm[i]];

        for ( IndexType j = 0; j < ilg[i]; j++ )
        {
            mValues[offset] *= static_cast<ValueType>( value );
            offset += dlg[j];
        }
    }
}

template<typename ValueType,typename OtherValueType>
void CUDAJDSUtils::scaleValue(
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    ValueType mValues[],
    const OtherValueType values[] )
{
    LAMA_LOG_INFO( logger, "scaleValue with numRows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    scaleValueKernel<<<dimGrid, dimBlock>>>( numRows, perm, ilg, dlg, mValues, values );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  checkDiagonalProperty                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

__global__
void checkDiagonalPropertyKernel( const IndexType numRows, bool *result, const IndexType *perm, const IndexType *ja )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        result[i] = ( ja[i] == perm[i] );
    }
}

bool CUDAJDSUtils::checkDiagonalProperty(
    const IndexType numDiagonals,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType perm[],
    const IndexType ja[],
    const IndexType dlg[] )
{
    LAMA_LOG_INFO( logger,
                   "checkDiagonalProperty with numDiagonals = " << numDiagonals << ", numRows = " << numRows << " and numColumns = " << numColumns )

    LAMA_CHECK_CUDA_ACCESS

    if ( numRows > 0 )
    {
        thrust::device_ptr<IndexType> dlgPtr( const_cast<IndexType*>( dlg ) );
        thrust::host_vector<IndexType> firstDlg( dlgPtr, dlgPtr + 1 );

        if ( firstDlg[0] < std::min( numDiagonals, numColumns ) )
        {
            return false;
        }

        thrust::device_ptr<bool> resultPtr = thrust::device_malloc<bool>( numDiagonals );
        thrust::fill( resultPtr, resultPtr + numDiagonals, false );

        bool *resultRawPtr = thrust::raw_pointer_cast( resultPtr );

        const int block_size = 256;
        dim3 dimBlock( block_size, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        checkDiagonalPropertyKernel<<<dimGrid, dimBlock>>>( numRows, resultRawPtr, perm, ja );

        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
        ;

        return thrust::reduce( resultPtr, resultPtr + numDiagonals, true, thrust::logical_and<bool>() );
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
void checkDescendingKernel( const IndexType n, const IndexType *array, bool *result )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i == 0 )
    {
        if ( n > 1 )
        {
            for ( IndexType i = 1; i < n; i++ )
            {
                if ( array[i] > array[i - 1] )
                {
                    result[0] = false;
                }
            }
        }
    }
}

bool CUDAJDSUtils::check(
    const IndexType numRows,
    const IndexType numValues,
    const IndexType numColumns,
    const IndexType ja[],
    const IndexType ilg[],
    const IndexType dlg[] )
{
    LAMA_LOG_INFO( logger, "check with numValues = " << numValues << ", numColumns = " << numColumns )

    LAMA_CHECK_CUDA_ACCESS

    if ( numRows > 0 )
    {
        thrust::device_ptr<IndexType> jaPtr( const_cast<IndexType*>( ja ) );

        bool error = false;

        error = thrust::transform_reduce( jaPtr, jaPtr + numValues, greaterThan<IndexType>( numColumns ), 0,
                                          thrust::logical_or<bool>() );
        if ( error )
        {
            return false;
        }

        thrust::device_ptr<bool> resultPtr = thrust::device_malloc<bool>( 1 );
        thrust::fill( resultPtr, resultPtr + 1, true );
        bool *resultRawPtr = thrust::raw_pointer_cast( resultPtr );

        thrust::device_ptr<IndexType> ilgPtr( const_cast<IndexType*>( ilg ) );
        thrust::host_vector<IndexType> ilgHost( ilgPtr, ilgPtr + 1 );

        const int block_size = 256;
        dim3 dimBlock( block_size, 1, 1 );
        dim3 dimGrid = makeGrid( 1, dimBlock.x );

        {
            checkDescendingKernel<<<dimGrid, dimBlock>>>( numRows, ilg, resultRawPtr );

            cudaStreamSynchronize( 0 );
            LAMA_CHECK_CUDA_ERROR
            ;

            thrust::host_vector<IndexType> result( resultPtr, resultPtr + 1 );

            if ( !result[0] )
            {
                return false;
            }
        }

        {
            checkDescendingKernel<<<dimGrid, dimBlock>>>( ilgHost[0], dlg, resultRawPtr );

            cudaStreamSynchronize( 0 );
            LAMA_CHECK_CUDA_ERROR
            ;

            thrust::host_vector<IndexType> result( resultPtr, resultPtr + 1 );

            if ( !result[0] )
            {
                return false;
            }
        }

        IndexType dlgSum = CUDAUtils::sum( dlg, ilgHost[0] );
        IndexType ilgSum = CUDAUtils::sum( ilg, numRows );

        if ( dlgSum != ilgSum )
        {
            return false;
        }

    }
    return true;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  ilg2dlg                                                           */
/* ------------------------------------------------------------------------------------------------------------------ */

__global__
void ilg2dlgKernel( IndexType *dlg, const IndexType numDiagonals, const IndexType *ilg, const IndexType numRows )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numDiagonals )
    {
        for ( IndexType j = 0; j < numRows; j++ )
        {
            if ( ilg[j] > i )
            {
                dlg[i]++;
            }
        }
    }
}

IndexType CUDAJDSUtils::ilg2dlg(
    IndexType dlg[],
    const IndexType numDiagonals,
    const IndexType ilg[],
    const IndexType numRows )
{
    LAMA_LOG_INFO( logger, "ilg2dlg with numDiagonals = " << numDiagonals << ", numRows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    if ( numDiagonals == 0 )
    {
        return 0;
    }

    // create device pointers and ilg sum

    thrust::device_ptr<IndexType> dlgPtr( const_cast<IndexType*>( dlg ) );
    thrust::device_ptr<IndexType> ilgPtr( const_cast<IndexType*>( ilg ) );
    thrust::fill( dlgPtr, dlgPtr + numDiagonals, 0 );
    IndexType sumIlg = thrust::reduce( ilgPtr, ilgPtr + numRows, 0, thrust::plus<IndexType>() );

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numDiagonals, dimBlock.x );

    ilg2dlgKernel<<<dimGrid, dimBlock>>>( dlg, numDiagonals, ilg, numRows );

    LAMA_CHECK_CUDA_ERROR

    return sumIlg;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  sortRows                                                          */
/* ------------------------------------------------------------------------------------------------------------------ */

void CUDAJDSUtils::sortRows( IndexType array[], IndexType perm[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "sort " << n << " rows by sizes" )

    // Note: this solution does not work on Tesla cards (doesent it?)

    LAMA_CHECK_CUDA_ACCESS
    thrust::device_ptr<IndexType> array_d( const_cast<IndexType*>( array ) );
    thrust::device_ptr<IndexType> perm_d( const_cast<IndexType*>( perm ) );

    thrust::stable_sort_by_key( array_d, array_d + n, perm_d, thrust::greater<IndexType>() );

    LAMA_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T1,typename T2>
__global__
void csr2jdsKernel(
    int* jdsJa,
    T1* jdsValues,
    const int* const jdsDlg,
    const int* const jdsIlg,
    const int* const jdsPerm,
    const int nrows,
    const int* const csrIa,
    const int* const csrJa,
    const T2* const csrValues )
{
    const int index = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( index < nrows )
    {
        int i = jdsPerm[index];
        int offset = index;
        for ( int jdsJJ = 0, csrJJ = csrIa[i]; jdsJJ < jdsIlg[index]; jdsJJ++, csrJJ++ )
        {
            jdsJa[offset] = csrJa[csrJJ];
            jdsValues[offset] = csrValues[csrJJ];
            offset += jdsDlg[jdsJJ]; // there is next value for row
        }
    }
}

template<typename JDSValueType,typename CSRValueType>
void CUDAJDSUtils::setCSRValues(
    IndexType jdsJA[],
    JDSValueType jdsValues[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType jdsDLG[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[] )
{
    LAMA_LOG_INFO( logger, "convert CSR to JDS, #rows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    csr2jdsKernel<<<dimGrid,dimBlock>>>( jdsJA, jdsValues, jdsDLG, jdsILG, jdsPerm, numRows, csrIA, csrJA, csrValues );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setInversePerm                                                    */
/* ------------------------------------------------------------------------------------------------------------------ */

void CUDAJDSUtils::setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "compute inverse perm, n = " << n )

    LAMA_CHECK_CUDA_ACCESS

    if ( n > 0 )
    {
        thrust::device_ptr<IndexType> inversePermPtr( const_cast<IndexType*>( inversePerm ) );
        thrust::device_ptr<IndexType> permPtr( const_cast<IndexType*>( perm ) );

        thrust::counting_iterator<IndexType> sequence( 0 );

        thrust::scatter( sequence, sequence + n, permPtr, inversePermPtr );

        LAMA_CHECK_CUDA_ERROR
        ;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename JDSValueType,typename CSRValueType>
__global__
void jds2csrKernel(
    IndexType *csrJA,
    CSRValueType *csrValues,
    const IndexType *csrIA,
    const IndexType numRows,
    const IndexType *jdsInversePerm,
    const IndexType *jdsILG,
    const IndexType *jdsDLG,
    const IndexType *jdsJA,
    const JDSValueType *jdsValues )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType ii = jdsInversePerm[i]; // where to find row i in JDS storage

        const IndexType numValuesInRow = jdsILG[ii];

        IndexType jdsOffset = ii; // run through input JDS data
        IndexType offset = csrIA[i]; // run through output data

        for ( IndexType jj = 0; jj < numValuesInRow; jj++ )
        {
            csrJA[offset + jj] = jdsJA[jdsOffset];
            csrValues[offset + jj] = static_cast<CSRValueType>( jdsValues[jdsOffset] );
            jdsOffset += jdsDLG[jj];
        }

    }
}

template<typename JDSValueType,typename CSRValueType>
void CUDAJDSUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType jdsInversePerm[],
    const IndexType jdsILG[],
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const JDSValueType jdsValues[] )
{
    LAMA_LOG_INFO( logger,
                   "get CSRValues<" << typeid( JDSValueType ).name() << ", " << typeid( CSRValueType ).name() << ">" << ", #rows = " << numRows )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR

    jds2csrKernel<<<dimGrid,dimBlock>>>( csrJA, csrValues, csrIA, numRows, jdsInversePerm, jdsILG, jdsDLG, jdsJA,
                                         jdsValues );
    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  xxxxx                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

/* --------------------------------------------------------------------------- */
/*                          Jacobi                                             */
/* --------------------------------------------------------------------------- */

texture<float,1> texJDSSXref;

texture<int2,1> texJDSDXref;

texture<int,1> texJDSdlgRef;

/* --------------------------------------------------------------------------- */

template<typename T,bool useTexture>
__inline__    __device__ T fetch_JDSx( const T* const x, const int i )
{
    return x[i];
}

template<bool useTexture,bool useSharedMemory>
__inline__ __device__
int fetch_JDSdlg( const int* const dlg_d, int[], const int i )
{
    return dlg_d[i];
}

template<>
__inline__ __device__
float fetch_JDSx<float,true>( const float* const, const int i )
{
    return tex1Dfetch( texJDSSXref, i );
}

template<>
__inline__ __device__
double fetch_JDSx<double,true>( const double* const, const int i )
{
    int2 v = tex1Dfetch( texJDSDXref, i );
    return __hiloint2double( v.y, v.x );
}

template<>
__inline__ __device__
int fetch_JDSdlg<true,false>( const int* const, int[], const int i )
{
    return tex1Dfetch( texJDSdlgRef, i );
}

template<>
__inline__ __device__
int fetch_JDSdlg<true,true>( const int* const, int dlg_sm[], const int i )
{
    return dlg_sm[i];
}

template<>
__inline__ __device__
int fetch_JDSdlg<false,true>( const int* const, int dlg_sm[], const int i )
{
    return dlg_sm[i];
}

template<typename T,bool useTexture,bool useSharedMem>
__global__
void jds_jacobi_kernel(
    const T* const jdsValues,
    const int* const jdsDlg,
    const int ndlg,
    const int* const jdsIlg,
    const int* const jdsJA,
    const int* const jdsPerm,
    const int numRows,
    const T* const rhs,
    T* const solution,
    const T* const oldSolution,
    const T omega )
{
    extern __shared__ int dlg[];
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( useSharedMem )
    {
        int k = threadIdx.x;
        while ( k < ndlg )
        {
            dlg[k] = jdsDlg[k];
            k += blockDim.x;
        }
        __syncthreads();
    }

    if ( i < numRows )
    {
        const int perm = jdsPerm[i];

        T temp = rhs[perm];

        const T aDiag = jdsValues[i];

        int pos = i + fetch_JDSdlg<useTexture,useSharedMem>( jdsDlg, dlg, 0 );
        const int rowEnd = jdsIlg[i];
        for ( int jj = 1; jj < rowEnd; ++jj )
        {
            temp -= jdsValues[pos] * fetch_JDSx<T,useTexture>( oldSolution, jdsJA[pos] );
            pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDlg, dlg, jj );
        }

        if ( omega == 0.5 )
        {
            solution[perm] = omega * ( fetch_JDSx<T,useTexture>( oldSolution, perm ) + temp / aDiag );
        }
        else if ( omega == 1.0 )
        {
            solution[perm] = temp / aDiag;
        }
        else
        {
            solution[perm] = omega * ( temp / aDiag ) + ( 1.0 - omega ) * fetch_JDSx<T,useTexture>( oldSolution, perm );
        }

    }
}

template<typename ValueType>
void CUDAJDSUtils::jacobi(
    ValueType solution[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsIlg[],
    const IndexType ndlg,
    const IndexType jdsDlg[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    SyncToken* syncToken )
{
    cudaStream_t stream = 0;

    LAMA_LOG_INFO( logger,
                   "jacobi<" << typeid(ValueType).name() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    LAMA_CHECK_CUDA_ACCESS

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        LAMA_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    const int block_size = ( numRows > 8191 ? 256 : 128 );
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    const bool useTexture = false; // lama_getUseTex_cuda();
    const bool useSharedMem = false; // maybe optimize later

    LAMA_LOG_DEBUG( logger, "useTexture = " << useTexture << ", useSharedMem = " << useSharedMem )

    if ( useTexture )
    {
        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSDXref, oldSolution ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSSXref, oldSolution ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
        }

        if ( !useSharedMem )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSdlgRef, jdsDlg ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
            LAMA_CUDA_RT_CALL(
                cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, true, false>, cudaFuncCachePreferL1 ),
                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, true, true>,cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
    }
    else
    {
        if ( !useSharedMem )
        {
            LAMA_CUDA_RT_CALL(
                cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, false, false>,cudaFuncCachePreferL1 ),
                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
        else
        {
            LAMA_CUDA_RT_CALL(
                cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, false, true>, cudaFuncCachePreferL1 ),
                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
    }

    if ( useTexture )
    {
        if ( !useSharedMem )
        {
            jds_jacobi_kernel<ValueType, true, false> <<<dimGrid, dimBlock, 0, stream>>>( jdsValues, jdsDlg, ndlg, jdsIlg, jdsJA, jdsPerm,
                    numRows, rhs, solution, oldSolution, omega );
        }
        else
        {
            const int sharedMemSize = ndlg * sizeof(int);
            jds_jacobi_kernel<ValueType, true, true> <<<dimGrid, dimBlock, sharedMemSize, stream>>>( jdsValues, jdsDlg, ndlg, jdsIlg,
                    jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega );
        }
    }
    else
    {
        if ( !useSharedMem )
        {
            jds_jacobi_kernel<ValueType, false, false> <<<dimGrid, dimBlock, 0, stream>>>( jdsValues, jdsDlg, ndlg, jdsIlg, jdsJA,
                    jdsPerm, numRows, rhs, solution, oldSolution, omega );
        }
        else
        {
            const int sharedMemSize = ndlg * sizeof(int);
            jds_jacobi_kernel<ValueType, false, true> <<<dimGrid, dimBlock, sharedMemSize, stream>>>( jdsValues, jdsDlg, ndlg, jdsIlg,
                    jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega);
        }
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_SJDSJACOBI_CUDAKERNEL_FAILED" );

    if ( useTexture )
    {

        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSDXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSSXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }

        if ( !useSharedMem )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSdlgRef ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }
    }

    if ( !syncToken )
    {
        cudaStreamSynchronize( stream );
    }
}

/* --------------------------------------------------------------------------- */
/*                          Jacobi halo                                        */
/* --------------------------------------------------------------------------- */

template<typename T,bool useTexture,bool useSharedMem>
__global__
void jds_jacobi_halo_kernel(
    const T* const diagonal,
    const T* const jdsValuesHalo,
    const int* const jdsDlgHalo,
    const int ndlg_halo,
    const int* const jdsIlgHalo,
    const int* const jdsJAHalo,
    const int* const jdsPermHalo,
    T* const solutionLocal,
    const T* const oldSolutionHalo,
    const T omega )
{
    extern __shared__ int dlg[];

    const int id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( useSharedMem )
    {
        int k = threadIdx.x;
        while ( k < ndlg_halo )
        {
            dlg[k] = jdsDlgHalo[k];
            k += blockDim.x;
        }
        __syncthreads();
    }

    if ( id < fetch_JDSdlg<useTexture,useSharedMem>( jdsDlgHalo, dlg, 0 ) )
    {
        T temp = 0.0;
        int pos = id;
        const int rowEnd = jdsIlgHalo[id];
        const int perm = jdsPermHalo[id];
        for ( int jj = 0; jj < rowEnd; ++jj )
        {
            temp += jdsValuesHalo[pos] * fetch_JDSx<T,useTexture>( oldSolutionHalo, jdsJAHalo[pos] );
            pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDlgHalo, dlg, jj );
        }

        const T aDiag = diagonal[perm];
        solutionLocal[perm] -= temp * omega / aDiag;
    }
}

template<typename ValueType>
void CUDAJDSUtils::jacobiHalo(
    ValueType solutionLocal[],
    const IndexType numRows,
    const ValueType diagonal[],
    const IndexType ndlg_halo,
    const IndexType jdsPermHalo[],
    const IndexType jdsIlgHalo[],
    const IndexType jdsDlgHalo[],
    const IndexType jdsJAHalo[],
    const ValueType jdsValuesHalo[],
    const ValueType oldSolutionHalo[],
    const ValueType omega,
    SyncToken* UNUSED(syncToken) )
{

    LAMA_LOG_INFO( logger,
                   "jacobiHalo<" << typeid(ValueType).name() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = ( numRows > 8191 ? 256 : 128 ) / 2;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x ); // TODO:numRows is too much...

    const bool useTexture = false; // lama_getUseTex_cuda();
    const bool useSharedMem = false; // maybe optimize later

    LAMA_LOG_DEBUG( logger, "useTexture = " << useTexture << ", useSharedMem = " << useSharedMem )

    if ( useTexture )
    {

        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSDXref, oldSolutionHalo ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" )
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSSXref, oldSolutionHalo ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
        }
        if ( !useSharedMem )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSdlgRef, jdsDlgHalo ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, true, false>, cudaFuncCachePreferL1),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, true, true>,cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }

    }
    else
    {
        if ( !useSharedMem )
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, false, false>, cudaFuncCachePreferL1),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
        else
        {
            LAMA_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, false, true>,cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }

    }

    if ( useTexture )
    {
        if ( !useSharedMem )
        {
            jds_jacobi_halo_kernel<ValueType, true, false> <<<dimGrid,dimBlock,0>>>( diagonal, jdsValuesHalo, jdsDlgHalo,
                    ndlg_halo, jdsIlgHalo, jdsJAHalo,
                    jdsPermHalo,
                    solutionLocal, oldSolutionHalo, omega);
        }
        else
        {
            const int sharedMemSize = ndlg_halo * sizeof(int);
            jds_jacobi_halo_kernel<ValueType, true, true> <<<dimGrid,dimBlock,sharedMemSize>>>( diagonal, jdsValuesHalo, jdsDlgHalo,
                    ndlg_halo, jdsIlgHalo, jdsJAHalo,
                    jdsPermHalo,
                    solutionLocal, oldSolutionHalo, omega);
        }

    }
    else
    {
        if ( !useSharedMem )
        {
            jds_jacobi_halo_kernel<ValueType, false, false> <<<dimGrid,dimBlock>>>( diagonal, jdsValuesHalo, jdsDlgHalo,
                    ndlg_halo, jdsIlgHalo, jdsJAHalo,
                    jdsPermHalo,
                    solutionLocal, oldSolutionHalo, omega);
        }
        else
        {
            const int sharedMemSize = ndlg_halo * sizeof(int);
            jds_jacobi_halo_kernel<ValueType, false, true> <<<dimGrid,dimBlock,sharedMemSize>>>( diagonal, jdsValuesHalo, jdsDlgHalo,
                    ndlg_halo, jdsIlgHalo, jdsJAHalo,
                    jdsPermHalo,
                    solutionLocal, oldSolutionHalo, omega);
        }
    }

    LAMA_CUDA_RT_CALL( cudaGetLastError(), "LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED" );
    LAMA_CUDA_RT_CALL( cudaStreamSynchronize(0), "LAMA_STATUS_CSRJACOBIHALO_CUDAKERNEL_FAILED" );

    if ( useTexture )
    {
        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSDXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSSXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }
        if ( !useSharedMem )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSdlgRef ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType,bool useTexture,bool useSharedMem>
__global__
void jdsgemvKernel(
    IndexType n,
    const ValueType alpha,
    const ValueType* const jdsValues,
    const IndexType* const jdsDlg,
    const IndexType ndlg,
    const IndexType* const jdsIlg,
    const IndexType* jdsJA,
    const IndexType* jdsPerm,
    const ValueType* x_d,
    const ValueType beta,
    const ValueType* y_d,
    ValueType* const result_d )
{
    extern __shared__ IndexType dlg[];
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( useSharedMem )
    {
        int k = threadIdx.x;
        while ( k < ndlg )
        {
            dlg[k] = jdsDlg[k];
            k += blockDim.x;
        }
        __syncthreads();
    }

    if ( i < n )
    {
        IndexType perm = jdsPerm[i];
        ValueType summand = 0.0;
        if ( beta != 0.0 )
        {
            summand = beta * y_d[perm];
        }

        ValueType value = 0.0;
        int k = i;
        for ( int jj = 0; jj < jdsIlg[i]; ++jj )
        {
            IndexType j = jdsJA[k];
            value += jdsValues[k] * fetch_JDSx<ValueType,useTexture>( x_d, j );
            k += fetch_JDSdlg<useTexture,useSharedMem>( jdsDlg, dlg, jj );
        }
//        for ( int jj = 0; jj < ndlg; ++jj )
//        {
//            const int incr = fetch_JDSdlg<useTexture,useSharedMem>( jdsDlg, dlg, jj );
//            if ( i < incr )
//            {
//                IndexType j = jdsJA[k];
//                value += jdsValues[k] * fetch_JDSx<ValueType,useTexture>( x_d, j );
//                k += incr;
//            }
//            else
//            {
//                break;
//            }
//        }
        result_d[perm] = alpha * value + summand;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAJDSUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType ndlg,
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    class SyncToken* /* syncToken */)
{

    LAMA_LOG_INFO( logger, "normalGEMV<" << typeid(ValueType).name() << ">" << ", #rows = " << numRows )

    LAMA_LOG_INFO(
        logger, "alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", result = " << result )

    const bool useTexture = false; //lama_getUseTex_cuda();
    const bool useSharedMem = false; // maybe optimize later

    LAMA_LOG_DEBUG( logger, "useTexture = " << useTexture << ", useSharedMem = " << useSharedMem )

    const int block_size = ( numRows > 8191 ? 256 : 128 );

    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    LAMA_CHECK_CUDA_ACCESS

    if ( useTexture )
    {
        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSDXref, x ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSSXref, x ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
        }
        if ( useSharedMem )
        {
            const int sharedMemSize = ndlg * sizeof(int);
            cudaFuncSetCacheConfig( jdsgemvKernel<ValueType,true,true>, cudaFuncCachePreferL1 );
            jdsgemvKernel<ValueType, true, true><<<dimGrid,dimBlock,sharedMemSize>>>
            ( numRows, alpha, jdsValues, jdsDLG, ndlg, jdsILG, jdsJA, jdsPerm, x, beta, y, result);
        }
        else // no sharedMem
        {
            LAMA_CUDA_RT_CALL( cudaBindTexture( NULL, texJDSdlgRef, jdsDLG ), "LAMA_STATUS_CUDA_BINDTEX_FAILED" );
            cudaFuncSetCacheConfig( jdsgemvKernel<ValueType,true,false>, cudaFuncCachePreferL1 );
            jdsgemvKernel<ValueType, true, false><<<dimGrid,dimBlock>>>
            ( numRows, alpha, jdsValues, jdsDLG, ndlg, jdsILG, jdsJA, jdsPerm, x, beta, y, result);
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSdlgRef ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }
        if ( sizeof(ValueType) == sizeof(double) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSDXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }
        else if ( sizeof(ValueType) == sizeof(float) )
        {
            LAMA_CUDA_RT_CALL( cudaUnbindTexture( texJDSSXref ), "LAMA_STATUS_CUDA_UNBINDTEX_FAILED" );
        }
    }
    else // no Texture cache
    {
        if ( useSharedMem )
        {
            const int sharedMemSize = ndlg * sizeof(int);
            cudaFuncSetCacheConfig( jdsgemvKernel<ValueType,false,true>, cudaFuncCachePreferL1 );
            jdsgemvKernel<ValueType, false, true><<<dimGrid,dimBlock,sharedMemSize>>>
            ( numRows, alpha, jdsValues, jdsDLG, ndlg, jdsILG, jdsJA, jdsPerm, x, beta, y, result);
        }
        else // no sharedMem
        {
            cudaFuncSetCacheConfig( jdsgemvKernel<ValueType,false,false>, cudaFuncCachePreferL1 );
            jdsgemvKernel<ValueType, false, false><<<dimGrid,dimBlock>>>
            ( numRows, alpha, jdsValues, jdsDLG, ndlg, jdsILG, jdsJA, jdsPerm, x, beta, y, result);
        }
    }

    LAMA_CHECK_CUDA_ERROR

    cudaStreamSynchronize( 0 );
}

/* --------------------------------------------------------------------------- */

void CUDAJDSUtils::setInterface( JDSUtilsInterface& JDSUtils )
{
    LAMA_LOG_INFO( logger, "set JDS routines for CUDA in Interface" )

    LAMA_INTERFACE_REGISTER( JDSUtils, sortRows )
    LAMA_INTERFACE_REGISTER( JDSUtils, checkDiagonalProperty )
    LAMA_INTERFACE_REGISTER( JDSUtils, check )
    LAMA_INTERFACE_REGISTER( JDSUtils, ilg2dlg )
    LAMA_INTERFACE_REGISTER( JDSUtils, setInversePerm )

    LAMA_INTERFACE_REGISTER_T( JDSUtils, setDiagonalWithScalar, float )
    LAMA_INTERFACE_REGISTER_T( JDSUtils, setDiagonalWithScalar, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, double, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, double, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, double, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobi, float )
    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobi, double )

    LAMA_INTERFACE_REGISTER_T( JDSUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( JDSUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobiHalo, float )
    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobiHalo, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDAJDSUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.JDSUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDAJDSUtils::initialized = registerInterface();


} // namespace lama
