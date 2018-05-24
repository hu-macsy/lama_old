/**
 * @file sparsekernel/cuda/CUDAJDSUtils.cu
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
 * @brief Implementation of JDS utilities with CUDA
 * @author Bea Hornef, Thomas Brandes
 * @date 04.07.2012
 */

// hpp
#include <scai/sparsekernel/cuda/CUDAJDSUtils.hpp>

// local library
#include <scai/sparsekernel/JDSKernelTrait.hpp>


// internal scai library
#include <scai/utilskernel/cuda/CUDAUtils.hpp>
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAUtils.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>

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

#include <functional>

using scai::tasking::CUDAStreamSyncToken;

namespace scai
{

using common::TypeTraits;
using common::CUDASettings;
using utilskernel::CUDAUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( CUDAJDSUtils::logger, "CUDA.JDSUtils" )

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
    T operator()( thrust::tuple<T, T> y )
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
/*                                                  getRow                                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void getRowKernel(
    ValueType* row,
    const IndexType i,
    const IndexType* ilg,
    const IndexType* dlg,
    const IndexType* ja,
    const ValueType* values )
{
    IndexType offset = 0;

    for ( IndexType j = 0; j < ilg[i]; j++ )
    {
        row[ja[i + offset]] = values[i + offset];
        offset += dlg[j];
    }
}

template<typename ValueType>
void CUDAJDSUtils::getRow(
    ValueType row[],
    const IndexType i,
    const IndexType numColumns,
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType ja[],
    const ValueType values[] )
{
    SCAI_LOG_INFO( logger, "getRow with i = " << i << ", numColumns = " << numColumns << " and numRows = " << numRows )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> rowPtr( row );
    thrust::device_ptr<IndexType> permPtr( const_cast<IndexType*>( perm ) );
    thrust::fill( rowPtr, rowPtr + numColumns, ValueType( 0 ) );
    thrust::counting_iterator<IndexType> sequence( 0 );
    // correct index with permutation array
    IndexType ii = thrust::transform_reduce(
                       thrust::make_zip_iterator( thrust::make_tuple( permPtr, sequence ) ),
                       thrust::make_zip_iterator( thrust::make_tuple( permPtr + numRows, sequence + numRows ) ),
                       identity<IndexType>( i ), 0, thrust::plus<IndexType>() );
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( 1, dimBlock.x );
    //TODO: find better CUDA / Thrust implementation
    getRowKernel <<< dimGrid, dimBlock>>>( row, ii, ilg, dlg, ja, values );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS:getRowKernel FAILED" )
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  scaleValue                                                        */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void scaleRowsKernel(
    ValueType* jdsValues,
    const IndexType numRows,
    const IndexType* perm,
    const IndexType* ilg,
    const IndexType* dlg,
    const ValueType* rowValues )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType offset = i;
        ValueType rowScale = rowValues[perm[i]];

        for ( IndexType j = 0; j < ilg[i]; j++ )
        {
            jdsValues[offset] *= rowScale;
            offset += dlg[j];
        }
    }
}

template<typename ValueType>
void CUDAJDSUtils::scaleRows(
    ValueType jdsValues[],
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    const ValueType rowValues[] )
{
    SCAI_LOG_INFO( logger, "scaleValue with numRows = " << numRows )
    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    scaleRowsKernel <<< dimGrid, dimBlock>>>( jdsValues, numRows, perm, ilg, dlg, rowValues );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS:scaleRowsKernel FAILED" )
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  ilg2dlg                                                           */
/* ------------------------------------------------------------------------------------------------------------------ */

__global__
void ilg2dlgKernel( IndexType* dlg, const IndexType numDiagonals, const IndexType* ilg, const IndexType numRows )
{
    // Entries in dlg filled every time there is a change in values of consecutive elements
    //   i:     0  1  2  3  4  5
    // ilg:     5  5  3  3  3  1
    // nd1:     5  5  3  3  3  1
    // nd2:     5  3  3  3  1  0
    //             x        x  x
    //             |        |  |->    6
    //             |        |---->       5  5
    //             |------------->             2   2
    // dlg:                           6  5  5  2   2
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType nd1 = ilg[i];
        IndexType nd2 = 0;

        if ( i + 1 < numRows )
        {
            nd2 = ilg[i + 1];
        }

        for ( IndexType j = nd2; j < nd1; j++ )
        {
            dlg[j] = i + 1;
        }
    }
}

IndexType CUDAJDSUtils::ilg2dlg(
    IndexType dlg[],
    const IndexType numDiagonals,
    const IndexType ilg[],
    const IndexType numRows )
{
    SCAI_REGION( "CUDA.JDS.ilg2dlg" )
    SCAI_LOG_INFO( logger, "ilg2dlg with numDiagonals = " << numDiagonals << ", numRows = " << numRows )
    SCAI_CHECK_CUDA_ACCESS

    if ( numDiagonals == 0 )
    {
        return 0;
    }

    // wrap raw pointer ilg to build sum, const_cast required, is safe
    thrust::device_ptr<IndexType> ilgPtr( const_cast<IndexType*>( ilg ) );
    IndexType sumIlg = thrust::reduce( ilgPtr, ilgPtr + numRows, 0, thrust::plus<IndexType>() );
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    ilg2dlgKernel <<< dimGrid, dimBlock>>>( dlg, numDiagonals, ilg, numRows );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: ilg2dlgKernel FAILED" )
    return sumIlg;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename JDSValueType, typename CSRValueType, bool useSharedMem>
__global__
void csr2jdsKernel(
    IndexType* jdsJa,
    JDSValueType* jdsValues,
    const IndexType* const jdsDLG,
    const IndexType ndlg,
    const IndexType* const jdsILG,
    const IndexType* const jdsPerm,
    const IndexType nrows,
    const IndexType* const csrIa,
    const IndexType* const csrJa,
    const CSRValueType* const csrValues )
{
    extern __shared__ IndexType dlg[];
    const IndexType iJDS = threadId( gridDim, blockIdx, blockDim, threadIdx );

    // copy DLG array into shared memory for faster access

    if ( useSharedMem )
    {
        int k = threadIdx.x;

        while ( k < ndlg )
        {
            dlg[k] = jdsDLG[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    if ( iJDS < nrows )
    {
        const IndexType iCSR = jdsPerm[iJDS]; // row index for CSR data
        const IndexType csrOffset = csrIa[iCSR];
        IndexType jdsOffset = iJDS;
        const IndexType numValuesInRow = jdsILG[iJDS];

        for ( IndexType jj = 0; jj < numValuesInRow; ++jj )
        {
            jdsJa[jdsOffset] = csrJa[csrOffset + jj];
            jdsValues[jdsOffset] = static_cast<JDSValueType>( csrValues[csrOffset + jj] );

            if ( useSharedMem )
            {
                jdsOffset += dlg[jj];
            }
            else
            {
                jdsOffset += jdsDLG[jj];
            }
        }
    }
}

template<typename JDSValueType, typename CSRValueType>
void CUDAJDSUtils::setCSRValues(
    IndexType jdsJA[],
    JDSValueType jdsValues[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType ndlg,
    const IndexType jdsDLG[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[] )
{
    // convert CSR data to JDS, ja and values
    SCAI_REGION( "CUDA.JDS.setCSR" )
    SCAI_LOG_INFO( logger, "convert CSR to JDS, #rows = " << numRows )
    SCAI_CHECK_CUDA_ACCESS
    bool useSharedMem = CUDASettings::useSharedMem();
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    SCAI_LOG_INFO( logger, "Start csr2jds_kernel<" << TypeTraits<JDSValueType>::id()
                   << ", " << TypeTraits<CSRValueType>::id()
                   << ", useSharedMem = " << useSharedMem
                   << "> ( nrows = " << numRows << ", ndiag = " << ndlg << " )" );

    if ( useSharedMem )
    {
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr2jdsKernel<JDSValueType, CSRValueType, true>,
                           cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        const int sharedMemSize = ndlg * sizeof( int );
        csr2jdsKernel<JDSValueType, CSRValueType, true> <<< dimGrid, dimBlock, sharedMemSize>>>(
            jdsJA, jdsValues, jdsDLG, ndlg, jdsILG, jdsPerm, numRows, csrIA, csrJA, csrValues );
    }
    else
    {
        SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr2jdsKernel<JDSValueType, CSRValueType, false>,
                           cudaFuncCachePreferL1 ),
                           "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        csr2jdsKernel<JDSValueType, CSRValueType, false> <<< dimGrid, dimBlock, 0>>>(
            jdsJA, jdsValues, jdsDLG, ndlg, jdsILG, jdsPerm, numRows, csrIA, csrJA, csrValues );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "csr2jdsKernel failed" );
    SCAI_LOG_INFO( logger, "Ready csr2jds_kernel<" << TypeTraits<JDSValueType>::id()
                   << ", " << TypeTraits<CSRValueType>::id() << " )" )
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  getCSRValues                                                      */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename JDSValueType, typename CSRValueType>
__global__
void jds2csrKernel(
    IndexType* csrJA,
    CSRValueType* csrValues,
    const IndexType* csrIA,
    const IndexType numRows,
    const IndexType* jdsInversePerm,
    const IndexType* jdsILG,
    const IndexType* jdsDLG,
    const IndexType* jdsJA,
    const JDSValueType* jdsValues )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType ii = jdsInversePerm[i]; // where to find row i in JDS storage
        const IndexType numValuesInRow = jdsILG[ii];
        IndexType jdsOffset = ii;// run through input JDS data
        IndexType offset = csrIA[i];// run through output data

        for ( IndexType jj = 0; jj < numValuesInRow; jj++ )
        {
            csrJA[offset + jj] = jdsJA[jdsOffset];
            csrValues[offset + jj] = static_cast<CSRValueType>( jdsValues[jdsOffset] );
            jdsOffset += jdsDLG[jj];
        }
    }
}

template<typename JDSValueType, typename CSRValueType>
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
    SCAI_REGION( "CUDA.JDS.getCSR" )
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<JDSValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" << ", #rows = " << numRows )

    if ( numRows < 1 )
    {
        return;  // do not launch kernel, invalid configuration
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    jds2csrKernel <<< dimGrid, dimBlock>>>( csrJA, csrValues, csrIA, numRows, jdsInversePerm, jdsILG, jdsDLG, jdsJA,
                                            jdsValues );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS:jds2csrKernel FAILED" )
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  xxxxx                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

template<bool useTexture, bool useSharedMemory>
__inline__ __device__
IndexType fetch_JDSdlg( const IndexType* const dlg_d, IndexType[], const IndexType i )
{
    return dlg_d[i];
}

template<>
__inline__ __device__
IndexType fetch_JDSdlg<true, false>( const IndexType* const dlg_d, IndexType[], const IndexType i )
{
    return fetchVectorX<IndexType, true>( dlg_d, i );
}

template<>
__inline__ __device__
IndexType fetch_JDSdlg<true, true>( const IndexType* const, IndexType dlg_sm[], const IndexType i )
{
    return dlg_sm[i];
}

template<>
__inline__ __device__
IndexType fetch_JDSdlg<false, true>( const IndexType* const, IndexType dlg_sm[], const IndexType i )
{
    return dlg_sm[i];
}

template<typename T, bool useTexture, bool useSharedMem>
__global__
void jds_jacobi_kernel(
    const T* const jdsValues,
    const IndexType* const jdsDLG,
    const IndexType ndlg,
    const IndexType* const jdsIlg,
    const IndexType* const jdsJA,
    const IndexType* const jdsPerm,
    const IndexType numRows,
    const T* const rhs,
    T* const solution,
    const T* const oldSolution,
    const T omega )
{
    extern __shared__ IndexType dlg[];
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < ndlg )
        {
            dlg[k] = jdsDLG[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    if ( i < numRows )
    {
        const IndexType perm = jdsPerm[i];
        T temp = rhs[perm];
        const T aDiag = jdsValues[i];
        IndexType pos = i + fetch_JDSdlg<useTexture, useSharedMem>( jdsDLG, dlg, 0 );
        const IndexType rowEnd = jdsIlg[i];

        for ( IndexType jj = 1; jj < rowEnd; ++jj )
        {
            temp -= jdsValues[pos] * fetchVectorX<T, useTexture>( oldSolution, jdsJA[pos] );
            pos += fetch_JDSdlg<useTexture, useSharedMem>( jdsDLG, dlg, jj );
        }

        if ( omega == 0.5 )
        {
            solution[perm] = omega * ( fetchVectorX<T, useTexture>( oldSolution, perm ) + temp / aDiag );
        }
        else if ( omega == 1.0 )
        {
            solution[perm] = temp / aDiag;
        }
        else
        {
            solution[perm] = omega * ( temp / aDiag ) + ( 1.0 - omega ) * fetchVectorX<T, useTexture>( oldSolution, perm );
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
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega )
{
    SCAI_REGION( "CUDA.JDS.jacobi" )
    cudaStream_t stream = 0;
    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )
    SCAI_CHECK_CUDA_ACCESS
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    const bool useTexture = CUDASettings::useTexture();

    const bool useSharedMem = CUDASettings::useSharedMem();

    const int blockSize = CUDASettings::getBlockSize( numRows );

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    SCAI_LOG_DEBUG( logger, "useTexture = " << useTexture << ", useSharedMem = " << useSharedMem )

    if ( useTexture )
    {
        vectorBindTexture( oldSolution );

        if ( !useSharedMem )
        {
            vectorBindTexture( jdsDLG );
            SCAI_CUDA_RT_CALL(
                cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, true, false>, cudaFuncCachePreferL1 ),
                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, true, true>, cudaFuncCachePreferL1 ),
                               "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
    }
    else
    {
        if ( !useSharedMem )
        {
            SCAI_CUDA_RT_CALL(
                cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, false, false>, cudaFuncCachePreferL1 ),
                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
        else
        {
            SCAI_CUDA_RT_CALL(
                cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, false, true>, cudaFuncCachePreferL1 ),
                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
        }
    }

    SCAI_LOG_INFO( logger, "Start jds_jacobi_kernel<" << TypeTraits<ValueType>::id()
                   << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">" );

    if ( useTexture )
    {
        if ( !useSharedMem )
        {
            jds_jacobi_kernel<ValueType, true, false> <<< dimGrid, dimBlock, 0, stream>>>(
                jdsValues, jdsDLG, ndlg, jdsIlg, jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega );
        }
        else
        {
            const int sharedMemSize = ndlg * sizeof( int );
            jds_jacobi_kernel<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream>>>(
                jdsValues, jdsDLG, ndlg, jdsIlg, jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega );
        }
    }
    else
    {
        if ( !useSharedMem )
        {
            jds_jacobi_kernel<ValueType, false, false> <<< dimGrid, dimBlock, 0, stream>>>(
                jdsValues, jdsDLG, ndlg, jdsIlg, jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega );
        }
        else
        {
            const int sharedMemSize = ndlg * sizeof( int );
            jds_jacobi_kernel<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream>>>(
                jdsValues, jdsDLG, ndlg, jdsIlg, jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega );
        }
    }

    SCAI_CUDA_RT_CALL( cudaGetLastError(), "jds_jacobi_kernel<" << TypeTraits<ValueType>::id()
                       << ", " << useTexture << ", " << useSharedMem << "> failed" )

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "JDS:jacobi_kernel failed" )
    }

    if ( useTexture )
    {
        if ( !syncToken )
        {
            vectorUnbindTexture( oldSolution );

            if ( !useSharedMem )
            {
                vectorUnbindTexture( jdsDLG );
            }
        }
        else
        {
            // synchronize by syncToken, delay unbind texture
            void ( *unbindV ) ( const ValueType* ) = &vectorUnbindTexture;
            void ( *unbindI ) ( const IndexType* ) = &vectorUnbindTexture;
            syncToken->pushRoutine( std::bind( unbindV, oldSolution ) );

            if ( !useSharedMem )
            {
                syncToken->pushRoutine( std::bind( unbindI, jdsDLG ) );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */
/*                          Jacobi halo                                        */
/* --------------------------------------------------------------------------- */

template<typename T, bool useTexture, bool useSharedMem>
__global__
void jds_jacobi_halo_kernel(
    const T* const diagonal,
    const T* const jdsValuesHalo,
    const IndexType* const jdsDLGHalo,
    const IndexType ndlg_halo,
    const IndexType* const jdsIlgHalo,
    const IndexType* const jdsJAHalo,
    const IndexType* const jdsPermHalo,
    T* const solutionLocal,
    const T* const oldSolutionHalo,
    const T omega )
{
    extern __shared__ IndexType dlg[];
    const IndexType id = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < ndlg_halo )
        {
            dlg[k] = jdsDLGHalo[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    if ( id < fetch_JDSdlg<useTexture, useSharedMem>( jdsDLGHalo, dlg, 0 ) )
    {
        T temp = 0.0;
        IndexType pos = id;
        const IndexType rowEnd = jdsIlgHalo[id];
        const IndexType perm = jdsPermHalo[id];

        for ( IndexType jj = 0; jj < rowEnd; ++jj )
        {
            temp += jdsValuesHalo[pos] * fetchVectorX<T, useTexture>( oldSolutionHalo, jdsJAHalo[pos] );
            pos += fetch_JDSdlg<useTexture, useSharedMem>( jdsDLGHalo, dlg, jj );
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
    const IndexType jdsDLGHalo[],
    const IndexType jdsJAHalo[],
    const ValueType jdsValuesHalo[],
    const ValueType oldSolutionHalo[],
    const ValueType omega )
{
    SCAI_REGION( "CUDA.JDS.jacobiHalo" )
    SCAI_LOG_INFO( logger, "jacobiHalo<" << TypeTraits<ValueType>::id() << ">"
                   << ", #rows = " << numRows << ", omega = " << omega )
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        COMMON_THROWEXCEPTION( "jacobiHalo not supported for aynchronous execution" )
    }

    SCAI_CHECK_CUDA_ACCESS
    const bool useTexture = CUDASettings::useTexture();
    const bool useSharedMem = CUDASettings::useSharedMem();
    const int blockSize = CUDASettings::getBlockSize( numRows );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x ); // TODO:numRows is too much...
    SCAI_LOG_DEBUG( logger, "useTexture = " << useTexture << ", useSharedMem = " << useSharedMem )

    if ( useTexture )
    {
        vectorBindTexture( oldSolutionHalo );

        if ( !useSharedMem )
        {
            vectorBindTexture( jdsDLGHalo );
        }
    }

    SCAI_LOG_INFO( logger, "Start jds_jacobi_halo_kernel<" << TypeTraits<ValueType>::id()
                   << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">" );

    if ( useTexture )
    {
        if ( !useSharedMem )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, true, false>,
                               cudaFuncCachePreferL1 ),
                               "cudaFuncSetCacheConfig jds_jacobi_halo_kernel<ValueType, true, false> failed" )
            jds_jacobi_halo_kernel<ValueType, true, false> <<< dimGrid, dimBlock, 0>>>(
                diagonal, jdsValuesHalo, jdsDLGHalo, ndlg_halo, jdsIlgHalo, jdsJAHalo,
                jdsPermHalo, solutionLocal, oldSolutionHalo, omega );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, true, true>,
                               cudaFuncCachePreferL1 ),
                               "cudaFuncSetCacheConfig jds_jacobi_halo_kernel<ValueType, true, true> failed" )
            const int sharedMemSize = ndlg_halo * sizeof( int );
            jds_jacobi_halo_kernel<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize>>>(
                diagonal, jdsValuesHalo, jdsDLGHalo, ndlg_halo, jdsIlgHalo, jdsJAHalo,
                jdsPermHalo, solutionLocal, oldSolutionHalo, omega );
        }
    }
    else
    {
        if ( !useSharedMem )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, false, false>,
                               cudaFuncCachePreferL1 ),
                               "cudaFuncSetCacheConfig jds_jacobi_halo_kernel<ValueType, false, false> failed" )
            jds_jacobi_halo_kernel<ValueType, false, false> <<< dimGrid, dimBlock>>>(
                diagonal, jdsValuesHalo, jdsDLGHalo, ndlg_halo, jdsIlgHalo, jdsJAHalo,
                jdsPermHalo, solutionLocal, oldSolutionHalo, omega );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, false, true>,
                               cudaFuncCachePreferL1 ),
                               "cudaFuncSetCacheConfig jds_jacobi_halo_kernel<ValueType, false, true> failed" )
            const int sharedMemSize = ndlg_halo * sizeof( int );
            jds_jacobi_halo_kernel<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize>>>(
                diagonal, jdsValuesHalo, jdsDLGHalo, ndlg_halo, jdsIlgHalo, jdsJAHalo,
                jdsPermHalo, solutionLocal, oldSolutionHalo, omega );
        }
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "jds_jacobi_halo_kernel" );

    if ( useTexture )
    {
        vectorUnbindTexture( oldSolutionHalo );

        if ( !useSharedMem )
        {
            vectorUnbindTexture( jdsDLGHalo );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__
void sparse_gemv_kernel(
    ValueType* const result_d,
    const ValueType* x_d,
    const ValueType alpha,
    const ValueType* const jdsValues,
    const IndexType* const jdsDLG,
    const IndexType* const jdsILG,
    const IndexType* jdsJA,
    const IndexType* jdsPerm,
    IndexType numNonZeroRows,
    const IndexType ndlg )
{
    extern __shared__ IndexType dlg[];
    const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( useSharedMem )
    {
        IndexType k = threadIdx.x;

        while ( k < ndlg )
        {
            dlg[k] = jdsDLG[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    if ( ii < numNonZeroRows )
    {
        IndexType i = jdsPerm[ii]; // row in matrix
        ValueType value = 0.0;
        IndexType pos = ii;// position in jdsJA, jdsValues
        IndexType ni = jdsILG[ii];// number entries in row

        for ( IndexType jj = 0; jj < ni; ++jj )
        {
            IndexType j = jdsJA[pos];
            value += jdsValues[pos] * fetchVectorX<ValueType, useTexture>( x_d, j );
            pos += fetch_JDSdlg<useTexture, useSharedMem>( jdsDLG, dlg, jj );
        }

        result_d[i] += alpha * value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture, bool useSharedMem>
__global__
void sparse_gevm_kernel(
    ValueType result[],
    const ValueType x[],
    const ValueType alpha,
    const ValueType jdsValues[],
    const IndexType jdsDLG[],
    const IndexType jdsILG[],
    const IndexType jdsJA[],
    const IndexType jdsPerm[],
    const IndexType numRows,
    const IndexType ndlg )
{
    extern __shared__ IndexType dlg[];

    if ( useSharedMem )
    {
        int k = threadIdx.x;

        while ( k < ndlg )
        {
            dlg[k] = jdsDLG[k];
            k += blockDim.x;
        }

        __syncthreads();
    }

    const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( ii < numRows )
    {
        IndexType i  = jdsPerm[ii];  // original row in matrix
        ValueType xi = x[i];         // same value used for all updates

        IndexType pos = ii;               // position in jdsJA, jdsValues
        IndexType ni  = jdsILG[ii];        // number entries in row

        for ( IndexType jj = 0; jj < ni; ++jj )
        {
            IndexType j = jdsJA[pos];
            ValueType v = alpha * jdsValues[pos] * xi;

            common::CUDAUtils::atomicAdd( &result[j], v );

            pos += fetch_JDSdlg<useTexture, useSharedMem>( jdsDLG, dlg, jj );
        }
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
    const IndexType numColumns,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType ndlg,
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    const common::MatrixOp op )
{
    SCAI_REGION( "CUDA.JDS.normalGEMV" )

    const IndexType nResult = common::isTranspose( op ) ? numColumns : numRows;

    SCAI_LOG_INFO( logger, "normalGEMV<" << TypeTraits<ValueType>::id() << ">"
                   << " result[ " << numRows << "] = " << alpha
                   << " * A( #jds_diags = " << ndlg << " ) * x + " << beta << " * y " )

    // set result = beta * y, not needed if beta == 1 and y == result

    if ( beta == scai::common::Constants::ONE && result == y )
    {
        SCAI_LOG_DEBUG( logger, "normalGEMV is inc, no init of result needed" )
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "normalGEMV, set result = " << beta << " * y " )

        // setScale also deals with y undefined for beta == 0

        CUDAUtils::binaryOpScalar( result, y, beta, nResult, common::BinaryOp::MULT, false );
    }

    sparseGEMV( result, alpha, x, numRows, jdsPerm, jdsILG, ndlg, jdsDLG, jdsJA, jdsValues, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAJDSUtils::sparseGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType ndlg,
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    const common::MatrixOp op )
{
    SCAI_REGION( "CUDA.JDS.sparseGEMV" )

    SCAI_LOG_INFO( logger, "sparseGEMV<" << TypeTraits<ValueType>::id() << ">"
                   << ", #rows = " << numRows << ", #diags = " << ndlg )

    if ( ndlg == 0 )
    {
        return; // nothing to do
    }

    IndexType nonEmptyRows = numRows;

    SCAI_CUDA_RT_CALL( cudaMemcpy( &nonEmptyRows, &jdsDLG[0], sizeof( IndexType ), cudaMemcpyDeviceToHost ),
                       "dlg[0] for number of non-empty rows" )

    const bool useTexture = CUDASettings::useTexture();
    const bool useSharedMem = CUDASettings::useSharedMem(); // maybe optimize later
    const int blockSize = CUDASettings::getBlockSize( nonEmptyRows );

    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( nonEmptyRows, dimBlock.x );

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;// default stream if no SyncToken is available
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    SCAI_LOG_INFO( logger, "Start jdsgemvSparseKernel<" << TypeTraits<ValueType>::id()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" 
                   << ", op = " << op )

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( useSharedMem )
        {
            const int sharedMemSize = ndlg * sizeof( int );
 
            if ( common::isTranspose( op ) )
            {
                sparse_gevm_kernel<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, nonEmptyRows, ndlg );
            }
            else
            {
                sparse_gemv_kernel<ValueType, true, true> <<< dimGrid, dimBlock, sharedMemSize, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, nonEmptyRows, ndlg );
           }
        }
        else // no sharedMem
        {
            vectorBindTexture( jdsDLG );
 
            if ( common::isTranspose( op ) )
            {
                sparse_gevm_kernel<ValueType, true, false> <<< dimGrid, dimBlock, 0, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, nonEmptyRows, ndlg );
            }
            else
            {
                sparse_gemv_kernel<ValueType, true, false> <<< dimGrid, dimBlock, 0, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, nonEmptyRows, ndlg );
            }
        }

        // skip the following in case of asynchronous execution

        if ( !syncToken )
        {
            // synchronize now, then unbind texture
            SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: gemvSparseKernel FAILED" )

            if ( !useSharedMem )
            {
                vectorUnbindTexture( jdsDLG );
            }

            vectorUnbindTexture( x );
        }
        else
        {
            // synchronize by syncToken, delay unbind texture
            void ( *unbindV ) ( const ValueType* ) = &vectorUnbindTexture;
            void ( *unbindI ) ( const IndexType* ) = &vectorUnbindTexture;

            if ( !useSharedMem )
            {
                syncToken->pushRoutine( std::bind( unbindI, jdsDLG ) );
            }

            syncToken->pushRoutine( std::bind( unbindV, x ) );
        }
    }
    else // no use of Texture cache
    {
        if ( useSharedMem )
        {
            const int sharedMemSize = ndlg * sizeof( int );

            if ( common::isTranspose( op ) )
            {
                cudaFuncSetCacheConfig( sparse_gevm_kernel<ValueType, false, true>, cudaFuncCachePreferL1 );
                sparse_gevm_kernel<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, nonEmptyRows, ndlg );
            }
            else
            {
                cudaFuncSetCacheConfig( sparse_gemv_kernel<ValueType, false, true>, cudaFuncCachePreferL1 );
                sparse_gemv_kernel<ValueType, false, true> <<< dimGrid, dimBlock, sharedMemSize, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, nonEmptyRows, ndlg );
            }
        }
        else // no use of sharedMem
        {
            if ( common::isTranspose( op ) )
            {
                cudaFuncSetCacheConfig( sparse_gevm_kernel<ValueType, false, false>, cudaFuncCachePreferL1 );
                sparse_gevm_kernel<ValueType, false, false> <<< dimGrid, dimBlock, 0, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, nonEmptyRows, ndlg );
            }
            else
            {
                cudaFuncSetCacheConfig( sparse_gemv_kernel<ValueType, false, false>, cudaFuncCachePreferL1 );
                sparse_gemv_kernel<ValueType, false, false> <<< dimGrid, dimBlock, 0, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, nonEmptyRows, ndlg );
            }
        }

        if ( !syncToken )
        {
            SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: gemvSparseKernel FAILED" )
        }
    }
}

/* --------------------------------------------------------------------------- */

void CUDAJDSUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_INFO( logger, "register JDSUtils CUDA-routines for CUDA at kernel registry [" << flag << "]" )
    KernelRegistry::set<JDSKernelTrait::ilg2dlg>( ilg2dlg, ctx, flag );
}

template<typename ValueType>
void CUDAJDSUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register JDSUtils CUDA-routines for CUDA at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<JDSKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getRow<ValueType> >( getRow, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::scaleRows<ValueType> >( scaleRows, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void CUDAJDSUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register JDSUtils CUDA-routines for CUDA at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    KernelRegistry::set<JDSKernelTrait::setCSRValues<ValueType, OtherValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDAJDSUtils::CUDAJDSUtils()
{
    SCAI_LOG_INFO( logger, "register JDSUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_CUDA_LIST, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAJDSUtils::~CUDAJDSUtils()
{
    SCAI_LOG_INFO( logger, "unregister JDSUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_CUDA_LIST, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAJDSUtils CUDAJDSUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
