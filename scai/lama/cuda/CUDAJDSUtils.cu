/**
 * @file CUDAJDSUtils.cpp
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
 * @brief Implementation of JDS utilities with CUDA
 * @author Bea Hornef, Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/cuda/CUDAJDSUtils.hpp>

// local library
#include <scai/lama/cuda/utils.cu.h>
#include <scai/lama/cuda/CUDAUtils.hpp>
#include <scai/lama/cuda/CUDASettings.hpp>

#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>

// internal scai library
#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/Assert.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/Constants.hpp>

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

// boost
#include <boost/preprocessor.hpp>

using namespace scai::hmemo;
using namespace scai::tasking;

namespace scai
{

using common::getScalarType;
using common::Constants;

namespace lama
{

    SCAI_LOG_DEF_LOGGER( CUDAJDSUtils::logger, "CUDA.JDSUtils" )

    /* ------------------------------------------------------------------------------------------------------------------ */

#include <scai/lama/cuda/CUDATexVector.hpp>

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
        SCAI_LOG_INFO( logger, "getRow with i = " << i << ", numColumns = " << numColumns << " and numRows = " << numRows )

        SCAI_CHECK_CUDA_ACCESS

        thrust::device_ptr<OtherValueType> rowPtr( row );
        thrust::device_ptr<IndexType> permPtr( const_cast<IndexType*>( perm ) );

        thrust::fill( rowPtr, rowPtr + numColumns, Constants<OtherValueType>::zero );

        thrust::counting_iterator<IndexType> sequence( Constants<IndexType>::zero );

        // correct index with permutation array
        IndexType ii = thrust::transform_reduce(
                        thrust::make_zip_iterator( thrust::make_tuple( permPtr, sequence ) ),
                        thrust::make_zip_iterator( thrust::make_tuple( permPtr + numRows, sequence + numRows ) ),
                        identity<IndexType>( i ), Constants<IndexType>::zero, thrust::plus<IndexType>() );

        const int blockSize = CUDASettings::getBlockSize();

        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( 1, dimBlock.x );

        //TODO: find better CUDA / Thrust implementation
        getRowKernel<<<dimGrid, dimBlock>>>( row, ii, ilg, dlg, ja, values );

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS:getRowKernel FAILED" )
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
        SCAI_CHECK_CUDA_ACCESS

        thrust::device_ptr<ValueType> resultPtr = thrust::device_malloc < ValueType > ( 1 );
        ValueType *resultRawPtr = thrust::raw_pointer_cast( resultPtr );

        const int blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
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
        SCAI_LOG_INFO( logger, "scaleValue with numRows = " << numRows )

        SCAI_CHECK_CUDA_ACCESS

        const int blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        scaleValueKernel<<<dimGrid, dimBlock>>>( numRows, perm, ilg, dlg, mValues, values );

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS:scaleValueKernel FAILED" )
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                                  checkDiagonalProperty                                             */
    /* ------------------------------------------------------------------------------------------------------------------ */

    __global__
    void checkDiagonalPropertyKernel(
                    bool *result,
                    const IndexType numRows,
                    const IndexType numColumns,
                    const IndexType nonEmptyRows,
                    const IndexType *perm,
                    const IndexType *ja )
    {
        const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

        if ( i >= numRows )
        {
            return;
        }

        const IndexType iRow = perm[i];

        if ( iRow >= numColumns )
        {
            // row does not count for diagonal

            return;
        }

        if ( i >= nonEmptyRows )
        {
            // iRow has no entries at all, ilg[i] is 0

            result[0] = false;
        }
        else if ( ja[i] != iRow )
        {
            result[0] = false;
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
        SCAI_LOG_INFO( logger, "checkDiagonalProperty with numDiagonals = " << numDiagonals
                        << ", numRows = " << numRows << " and numColumns = " << numColumns )

        SCAI_CHECK_CUDA_ACCESS

        if ( numRows <= 0 )
        {
            return false;
        }

        if ( numDiagonals <= 0 )
        {
            return false;
        }

        // now it is sure that dlg, perm and ja are not empty

        const IndexType diagSize = std::min( numRows, numColumns );

        IndexType nonEmptyRows = 0;

        SCAI_CUDA_RT_CALL( cudaMemcpy( &nonEmptyRows, &dlg[0], sizeof( IndexType ), cudaMemcpyDeviceToHost ),
                        "get number of non-zero rows from dlg" );

        // Be careful: numDiagonals has nothing to do with size of diagonal

        if ( nonEmptyRows < diagSize )
        {
            return false;
        }

        bool* d_hasProperty; // will be ptr to device version of hasProperty

        bool hasProperty = true;

        SCAI_CUDA_RT_CALL( cudaMalloc( (void**) &d_hasProperty, sizeof( bool ) ),
                        "allocate 4 bytes on the device for the result of hasDiagonalProperty_kernel" )

        SCAI_CUDA_RT_CALL( cudaMemcpy( d_hasProperty, &hasProperty, sizeof( bool ), cudaMemcpyHostToDevice ),
                        "copy flag for diagonalProperty to device" )

        const int blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        checkDiagonalPropertyKernel<<<dimGrid, dimBlock>>>( d_hasProperty,
                        numRows, numColumns, nonEmptyRows,
                        perm, ja );

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS:checkDiagonalPropertyKernel FAILED" )

        SCAI_CUDA_RT_CALL( cudaMemcpy( &hasProperty, d_hasProperty, sizeof( bool ), cudaMemcpyDeviceToHost ),
                        "copy flag for diagonalProperty to host" )

        SCAI_CUDA_RT_CALL( cudaFree( d_hasProperty ),
                        "free result var for diagonal property" )

        return hasProperty;
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                                  ilg2dlg                                                           */
    /* ------------------------------------------------------------------------------------------------------------------ */

    __global__
    void ilg2dlgKernel( IndexType *dlg, const IndexType numDiagonals, const IndexType *ilg, const IndexType numRows )
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
        SCAI_REGION( "CUDA.JDS:dlg<-ilg" )

        SCAI_LOG_INFO( logger, "ilg2dlg with numDiagonals = " << numDiagonals << ", numRows = " << numRows )

        SCAI_CHECK_CUDA_ACCESS

        if ( numDiagonals == 0 )
        {
            return 0;
        }

        // wrap raw pointer ilg to build sum, const_cast required, is safe 

        thrust::device_ptr<IndexType> ilgPtr( const_cast<IndexType*>( ilg ) );

        IndexType sumIlg = thrust::reduce( ilgPtr, ilgPtr + numRows, Constants<IndexType>::zero, thrust::plus<IndexType>() );

        const int blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        ilg2dlgKernel<<<dimGrid, dimBlock>>>( dlg, numDiagonals, ilg, numRows );

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: ilg2dlgKernel FAILED" )

        return sumIlg;
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                                  sortRows                                                          */
    /* ------------------------------------------------------------------------------------------------------------------ */

    void CUDAJDSUtils::sortRows( IndexType array[], IndexType perm[], const IndexType n )
    {
        SCAI_REGION( "CUDA.JDS:sortRows" )

        SCAI_LOG_INFO( logger, "sort " << n << " rows by sizes" )

        SCAI_CHECK_CUDA_ACCESS

        thrust::device_ptr<IndexType> array_d( array );
        thrust::device_ptr<IndexType> perm_d( perm );

        // stable sort, descending order, so override default comparison

        thrust::stable_sort_by_key( array_d, array_d + n, perm_d, thrust::greater<IndexType>() );

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize(0), "JDS: synchronize for sortRows FAILED" )
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
        extern __shared__ int dlg[];

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

        SCAI_REGION( "CUDA.JDS<-CSR_values" )

        SCAI_LOG_INFO( logger, "convert CSR to JDS, #rows = " << numRows )

        SCAI_CHECK_CUDA_ACCESS

        bool useSharedMem = CUDASettings::useSharedMem();

        const int blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        SCAI_LOG_INFO( logger, "Start csr2jds_kernel<" << getScalarType<JDSValueType>()
                        << ", " << getScalarType<CSRValueType>()
                        << ", useSharedMem = " << useSharedMem
                        << "> ( nrows = " << numRows << ", ndiag = " << ndlg << " )" );

        if ( useSharedMem )
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr2jdsKernel<JDSValueType, CSRValueType, true>,
                                            cudaFuncCachePreferL1 ),
                            "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );

            const int sharedMemSize = ndlg * sizeof(int);

            csr2jdsKernel<JDSValueType, CSRValueType, true><<<dimGrid, dimBlock, sharedMemSize>>>(
                            jdsJA, jdsValues, jdsDLG, ndlg, jdsILG, jdsPerm, numRows, csrIA, csrJA, csrValues );
        }
        else
        {
            SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( csr2jdsKernel<JDSValueType, CSRValueType, false>,
                                            cudaFuncCachePreferL1 ),
                            "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );

            csr2jdsKernel<JDSValueType, CSRValueType, false><<<dimGrid, dimBlock, 0>>>(
                            jdsJA, jdsValues, jdsDLG, ndlg, jdsILG, jdsPerm, numRows, csrIA, csrJA, csrValues );
        }

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "csr2jdsKernel failed" );

        SCAI_LOG_INFO( logger, "Ready csr2jds_kernel<" << getScalarType<JDSValueType>()
                        << ", " << getScalarType<CSRValueType>() << " )" )
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                                  setInversePerm                                                    */
    /* ------------------------------------------------------------------------------------------------------------------ */

    void CUDAJDSUtils::setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n )
    {
        SCAI_LOG_INFO( logger, "compute inverse perm, n = " << n )

        SCAI_CHECK_CUDA_ACCESS

        if ( n > 0 )
        {
            thrust::device_ptr<IndexType> inversePermPtr( const_cast<IndexType*>( inversePerm ) );
            thrust::device_ptr<IndexType> permPtr( const_cast<IndexType*>( perm ) );

            thrust::counting_iterator<IndexType> sequence( Constants<IndexType>::zero );

            thrust::scatter( sequence, sequence + n, permPtr, inversePermPtr );

            SCAI_CHECK_CUDA_ERROR
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
        SCAI_REGION( "CUDA.JDS->CSR_values" )

        SCAI_LOG_INFO( logger,
                        "get CSRValues<" << getScalarType<JDSValueType>() << ", " << getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows )

        SCAI_CHECK_CUDA_ACCESS

        const int blockSize = CUDASettings::getBlockSize();
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        jds2csrKernel<<<dimGrid,dimBlock>>>( csrJA, csrValues, csrIA, numRows, jdsInversePerm, jdsILG, jdsDLG, jdsJA,
                        jdsValues );

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS:jds2csrKernel FAILED" )
    }

    /* ------------------------------------------------------------------------------------------------------------------ */
    /*                                                  xxxxx                                                             */
    /* ------------------------------------------------------------------------------------------------------------------ */

    template<bool useTexture, bool useSharedMemory>
    __inline__ __device__
    int fetch_JDSdlg( const int* const dlg_d, int[], const int i )
    {
        return dlg_d[i];
    }

    template<>
    __inline__ __device__
    int fetch_JDSdlg<true, false>( const int* const dlg_d, int[], const int i )
    {
        return fetchVectorX<int, true>( dlg_d, i );
    }

    template<>
    __inline__ __device__
    int fetch_JDSdlg<true, true>( const int* const, int dlg_sm[], const int i )
    {
        return dlg_sm[i];
    }

    template<>
    __inline__ __device__
    int fetch_JDSdlg<false, true>( const int* const, int dlg_sm[], const int i )
    {
        return dlg_sm[i];
    }

    template<typename T,bool useTexture,bool useSharedMem>
    __global__
    void jds_jacobi_kernel(
                    const T* const jdsValues,
                    const int* const jdsDLG,
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
                dlg[k] = jdsDLG[k];
                k += blockDim.x;
            }
            __syncthreads();
        }

        if ( i < numRows )
        {
            const int perm = jdsPerm[i];

            T temp = rhs[perm];

            const T aDiag = jdsValues[i];

            int pos = i + fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, 0 );
            const int rowEnd = jdsIlg[i];
            for ( int jj = 1; jj < rowEnd; ++jj )
            {
                temp -= jdsValues[pos] * fetchVectorX<T,useTexture>( oldSolution, jdsJA[pos] );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            if ( omega == 0.5 )
            {
                solution[perm] = omega * ( fetchVectorX<T,useTexture>( oldSolution, perm ) + temp / aDiag );
            }
            else if ( omega == 1.0 )
            {
                solution[perm] = temp / aDiag;
            }
            else
            {
                solution[perm] = omega * ( temp / aDiag ) + ( 1.0 - omega ) * fetchVectorX<T,useTexture>( oldSolution, perm );
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
                    const ValueType omega,
                    SyncToken* syncToken )
    {
        SCAI_REGION( "CUDA.JDS.jacobi" )

        cudaStream_t stream = 0;

        SCAI_LOG_INFO( logger,
                        "jacobi<" << getScalarType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

        SCAI_CHECK_CUDA_ACCESS

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
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
                SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, true, true>,cudaFuncCachePreferL1 ),
                                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
            }
        }
        else
        {
            if ( !useSharedMem )
            {
                SCAI_CUDA_RT_CALL(
                                cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, false, false>,cudaFuncCachePreferL1 ),
                                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
            }
            else
            {
                SCAI_CUDA_RT_CALL(
                                cudaFuncSetCacheConfig( jds_jacobi_kernel<ValueType, false, true>, cudaFuncCachePreferL1 ),
                                "LAMA_STATUS_CUDA_FUNCSETCACHECONFIG_FAILED" );
            }
        }

        SCAI_LOG_INFO( logger, "Start jds_jacobi_kernel<" << getScalarType<ValueType>()
                        << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">" );

        if ( useTexture )
        {
            if ( !useSharedMem )
            {
                jds_jacobi_kernel<ValueType, true, false> <<<dimGrid, dimBlock, 0, stream>>>(
                                jdsValues, jdsDLG, ndlg, jdsIlg, jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega );
            }
            else
            {
                const int sharedMemSize = ndlg * sizeof(int);
                jds_jacobi_kernel<ValueType, true, true> <<<dimGrid, dimBlock, sharedMemSize, stream>>>(
                                jdsValues, jdsDLG, ndlg, jdsIlg, jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega );
            }
        }
        else
        {
            if ( !useSharedMem )
            {
                jds_jacobi_kernel<ValueType, false, false> <<<dimGrid, dimBlock, 0, stream>>>(
                                jdsValues, jdsDLG, ndlg, jdsIlg, jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega );
            }
            else
            {
                const int sharedMemSize = ndlg * sizeof(int);
                jds_jacobi_kernel<ValueType, false, true> <<<dimGrid, dimBlock, sharedMemSize, stream>>>(
                                jdsValues, jdsDLG, ndlg, jdsIlg, jdsJA, jdsPerm, numRows, rhs, solution, oldSolution, omega);
            }
        }

        SCAI_CUDA_RT_CALL( cudaGetLastError(), "jds_jacobi_kernel<" << getScalarType<ValueType>()
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

                syncToken->pushRoutine( common::bind( unbindV, oldSolution ) );

                if ( !useSharedMem )
                {
                    syncToken->pushRoutine( common::bind( unbindI, jdsDLG ) );
                }
            }
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
                    const int* const jdsDLGHalo,
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
                dlg[k] = jdsDLGHalo[k];
                k += blockDim.x;
            }
            __syncthreads();
        }

        if ( id < fetch_JDSdlg<useTexture,useSharedMem>( jdsDLGHalo, dlg, 0 ) )
        {
            T temp = 0.0;
            int pos = id;
            const int rowEnd = jdsIlgHalo[id];
            const int perm = jdsPermHalo[id];
            for ( int jj = 0; jj < rowEnd; ++jj )
            {
                temp += jdsValuesHalo[pos] * fetchVectorX<T,useTexture>( oldSolutionHalo, jdsJAHalo[pos] );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLGHalo, dlg, jj );
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
                    const ValueType omega,
                    SyncToken* syncToken )
    {
        SCAI_REGION( "CUDA.JDS.jacobiHalo" )

        SCAI_LOG_INFO( logger, "jacobiHalo<" << getScalarType<ValueType>() << ">"
                        << ", #rows = " << numRows << ", omega = " << omega )

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

        SCAI_LOG_INFO( logger, "Start jds_jacobi_halo_kernel<" << getScalarType<ValueType>()
                        << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">" );

        if ( useTexture )
        {
            if ( !useSharedMem )
            {
                SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, true, false>,
                                                cudaFuncCachePreferL1),
                                "cudaFuncSetCacheConfig jds_jacobi_halo_kernel<ValueType, true, false> failed" )

                jds_jacobi_halo_kernel<ValueType, true, false> <<<dimGrid, dimBlock, 0>>>(
                                diagonal, jdsValuesHalo, jdsDLGHalo, ndlg_halo, jdsIlgHalo, jdsJAHalo,
                                jdsPermHalo, solutionLocal, oldSolutionHalo, omega);
            }
            else
            {
                SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, true, true>,
                                                cudaFuncCachePreferL1),
                                "cudaFuncSetCacheConfig jds_jacobi_halo_kernel<ValueType, true, true> failed" )

                const int sharedMemSize = ndlg_halo * sizeof(int);

                jds_jacobi_halo_kernel<ValueType, true, true> <<<dimGrid, dimBlock, sharedMemSize>>>(
                                diagonal, jdsValuesHalo, jdsDLGHalo, ndlg_halo, jdsIlgHalo, jdsJAHalo,
                                jdsPermHalo, solutionLocal, oldSolutionHalo, omega);
            }

        }
        else
        {
            if ( !useSharedMem )
            {
                SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, false, false>,
                                                cudaFuncCachePreferL1),
                                "cudaFuncSetCacheConfig jds_jacobi_halo_kernel<ValueType, false, false> failed" )

                jds_jacobi_halo_kernel<ValueType, false, false> <<<dimGrid,dimBlock>>>(
                                diagonal, jdsValuesHalo, jdsDLGHalo, ndlg_halo, jdsIlgHalo, jdsJAHalo,
                                jdsPermHalo, solutionLocal, oldSolutionHalo, omega);
            }
            else
            {
                SCAI_CUDA_RT_CALL( cudaFuncSetCacheConfig( jds_jacobi_halo_kernel<ValueType, false, true>,
                                                cudaFuncCachePreferL1),
                                "cudaFuncSetCacheConfig jds_jacobi_halo_kernel<ValueType, false, true> failed" )

                const int sharedMemSize = ndlg_halo * sizeof(int);

                jds_jacobi_halo_kernel<ValueType, false, true> <<<dimGrid, dimBlock, sharedMemSize>>>(
                                diagonal, jdsValuesHalo, jdsDLGHalo, ndlg_halo, jdsIlgHalo, jdsJAHalo,
                                jdsPermHalo, solutionLocal, oldSolutionHalo, omega);
            }
        }

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize(0), "jds_jacobi_halo_kernel" );

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

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gemv_kernel(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType alpha,
                    const ValueType beta,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numRows,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( ii < numRows )
        {
            IndexType i = jdsPerm[ii]; // row in matrix

            ValueType summand = beta * y_d[i];

            ValueType value = 0.0;

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] = alpha * value + summand;
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gemv_kernel_alpha_one_beta_one(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numRows,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( ii < numRows )
        {
            IndexType i = jdsPerm[ii]; // row in matrix

            ValueType summand = y_d[i];

            ValueType value = 0.0;

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] = value + summand;
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gemv_kernel_alpha_one_beta_zero(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numRows,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( ii < numRows )
        {
            IndexType i = jdsPerm[ii]; // row in matrix

            ValueType value = 0.0;

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] = value;
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void assign_kernel(
                    ValueType* const result_d,
                    const ValueType* y_d,
                    const IndexType* jdsPerm,
                    IndexType numRows)
    {
        const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

        if ( ii < numRows )
        {
            IndexType i = jdsPerm[ii]; // row in matrix
            result_d[i] = y_d[i];
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gemv_kernel_alpha_one(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType beta,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numRows,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( ii < numRows )
        {
            IndexType i = jdsPerm[ii]; // row in matrix

            ValueType summand = beta * y_d[i];

            ValueType value = 0.0;

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] = value + summand;
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void scale_kernel(
                    ValueType* const result_d,
                    const ValueType* y_d,
                    const ValueType beta,
                    const IndexType* jdsPerm,
                    IndexType numRows )
    {
        const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

        if ( ii < numRows )
        {
            IndexType i = jdsPerm[ii]; // row in matrix
            result_d[i] = beta * y_d[i];
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gemv_kernel_beta_one(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType alpha,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numRows,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( ii < numRows )
        {
            IndexType i = jdsPerm[ii]; // row in matrix

            ValueType summand = y_d[i];

            ValueType value = 0.0;

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] = alpha * value + summand;
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gemv_kernel_beta_zero(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType alpha,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numRows,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType ii = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( ii < numRows )
        {
            IndexType i = jdsPerm[ii]; // row in matrix

            ValueType value = 0.0;

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] = alpha * value;
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
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
            int k = threadIdx.x;
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

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] += alpha * value;
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void sparse_gemv_kernel_alpha_one(
                    ValueType* const result_d,
                    const ValueType* x_d,
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
            int k = threadIdx.x;
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

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] += value;
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void sparse_gevm_kernel(
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
            int k = threadIdx.x;
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

            int pos = ii;// position in jdsJA, jdsValues

            int ni = jdsILG[ii];// number entries in row

            for ( int jj = 0; jj < ni; ++jj )
            {
                IndexType j = jdsJA[pos];
                value += jdsValues[pos] * fetchVectorX<ValueType,useTexture>( x_d, j );
                pos += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
            }

            result_d[i] += alpha * value;
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
                    SyncToken* syncToken )
    {
        if ( ( beta == Constants<ValueType>::one ) && ( result == y ) )
        {
            // result = alpha * A * x + beta * y ->  result += alpha * A * x

            sparseGEMV( result, alpha, x, numRows, jdsPerm, jdsILG, ndlg, jdsDLG, jdsJA, jdsValues, syncToken );

            return;
        }

        SCAI_REGION( "CUDA.JDS.normalGEMV" )

        SCAI_LOG_INFO( logger, "normalGEMV<" << getScalarType<ValueType>() << ">"
                        << " result[ " << numRows << "] = " << alpha
                        << " * A( #jds_diags = " << ndlg << " ) * x + " << beta << " * y " )

        SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )

        const bool useTexture = CUDASettings::useTexture();
        const bool useSharedMem = CUDASettings::useSharedMem();
        const int blockSize = CUDASettings::getBlockSize( numRows );

        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        SCAI_CHECK_CUDA_ACCESS

        cudaStream_t stream = 0; // default stream if no SyncToken is available

        SCAI_LOG_INFO( logger, "Start normal_gemv_kernel<" << getScalarType<ValueType>()
                        << "> <<< blockSize = " << blockSize << ", stream = " << stream
                        << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
        }

        int sharedMemSize = 0;
        if ( useSharedMem )
        {
            sharedMemSize = ndlg * sizeof(int);
        }

        if ( useTexture )
        {
            vectorBindTexture( x );

            if ( useSharedMem )
            {
                if( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_alpha_one_beta_one<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::zero )
                {
                    normal_gemv_kernel_alpha_one_beta_zero<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero && beta == Constants<ValueType>::one )
                {
                    assign_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, jdsPerm, numRows );
                }
                else if ( alpha == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_alpha_one<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero )
                {
                    scale_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, beta, jdsPerm, numRows);
                }
                else if ( beta == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_beta_one<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( beta == Constants<ValueType>::zero )
                {
                    normal_gemv_kernel_beta_zero<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else
                {
                    normal_gemv_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
            }
            else // no sharedMem
            {
                vectorBindTexture( jdsDLG );

                if( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_alpha_one_beta_one<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::zero )
                {
                    normal_gemv_kernel_alpha_one_beta_zero<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero && beta == Constants<ValueType>::one )
                {
                    assign_kernel<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, jdsPerm, numRows );
                }
                else if ( alpha == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_alpha_one<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero )
                {
                    scale_kernel<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, beta, jdsPerm, numRows );
                }
                else if ( beta == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_beta_one<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( beta == Constants<ValueType>::zero )
                {
                    normal_gemv_kernel_beta_zero<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else
                {
                    normal_gemv_kernel<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
            }

            if ( !syncToken )
            {
                // synchronize here and unbind texture

                SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: gemvKernel FAILED" )

                vectorUnbindTexture( x );

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

                syncToken->pushRoutine( common::bind( unbindV, x ) );

                if ( !useSharedMem )
                {
                    syncToken->pushRoutine( common::bind( unbindI, jdsDLG ) );
                }
            }
        }
        else // no Texture cache
        {
            if ( useSharedMem )
            {
                if( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_alpha_one_beta_one<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::zero )
                {
                    normal_gemv_kernel_alpha_one_beta_zero<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero && beta == Constants<ValueType>::one )
                {
                    assign_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, jdsPerm, numRows );
                }
                else if ( alpha == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_alpha_one<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero )
                {
                    scale_kernel<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, beta, jdsPerm, numRows );
                }
                else if ( beta == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_beta_one<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( beta == Constants<ValueType>::zero )
                {
                    normal_gemv_kernel_beta_zero<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else
                {
                    normal_gemv_kernel<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
            }
            else // no sharedMem
            {
                if( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_alpha_one_beta_one<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::zero )
                {
                    normal_gemv_kernel_alpha_one_beta_zero<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero && beta == Constants<ValueType>::one )
                {
                    assign_kernel<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, jdsPerm, numRows );
                }
                else if ( alpha == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_alpha_one<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero )
                {
                    scale_kernel<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, beta, jdsPerm, numRows );
                }
                else if ( beta == Constants<ValueType>::one )
                {
                    normal_gemv_kernel_beta_one<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else if ( beta == Constants<ValueType>::zero )
                {
                    normal_gemv_kernel_beta_zero<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
                else
                {
                    normal_gemv_kernel<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
                }
            }
            if ( !syncToken )
            {
                SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: gemvKernel FAILED" )
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gevm_kernel(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType alpha,
                    const ValueType beta,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numColumns,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType k = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( k < numColumns )
        {
            if( ndlg > 0 )
            {
                ValueType value = 0.0;
                ValueType summand = 0.0;

                if ( beta != 0.0 )
                {
                    summand = beta * y_d[k];
                }

                int nonEmptyRows = jdsDLG[0];
                for ( int ii = 0; ii < nonEmptyRows; ++ii ) // through the rows of the first diagonal
                {
                    IndexType i = jdsPerm[ii]; // row in matrix

                    int off = ii;// position in jdsJA, jdsValues

                    int ni = jdsILG[ii];// number entries in row

                    for ( int jj = 0; jj < ni; ++jj )
                    {
                        IndexType j = jdsJA[off];
                        if ( j == k )
                        {
                            value += jdsValues[off] * fetchVectorX<ValueType,useTexture>( x_d, i );
                        }
                        off += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
                    }
                }

                result_d[k] = alpha * value + summand;
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gevm_kernel_alpha_one_beta_one(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numColumns,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType k = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( k < numColumns )
        {
            if( ndlg > 0 )
            {
                ValueType value = 0.0;
                ValueType summand = y_d[k];

                int nonEmptyRows = jdsDLG[0];
                for ( int ii = 0; ii < nonEmptyRows; ++ii ) // through the rows of the first diagonal
                {
                    IndexType i = jdsPerm[ii]; // row in matrix

                    int off = ii;// position in jdsJA, jdsValues

                    int ni = jdsILG[ii];// number entries in row

                    for ( int jj = 0; jj < ni; ++jj )
                    {
                        IndexType j = jdsJA[off];
                        if ( j == k )
                        {
                            value += jdsValues[off] * fetchVectorX<ValueType,useTexture>( x_d, i );
                        }
                        off += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
                    }
                }

                result_d[k] = value + summand;
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gevm_kernel_alpha_one_beta_zero(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numColumns,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType k = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( k < numColumns )
        {
            if( ndlg > 0 )
            {
                ValueType value = 0.0;

                int nonEmptyRows = jdsDLG[0];
                for ( int ii = 0; ii < nonEmptyRows; ++ii ) // through the rows of the first diagonal
                {
                    IndexType i = jdsPerm[ii]; // row in matrix

                    int off = ii;// position in jdsJA, jdsValues

                    int ni = jdsILG[ii];// number entries in row

                    for ( int jj = 0; jj < ni; ++jj )
                    {
                        IndexType j = jdsJA[off];
                        if ( j == k )
                        {
                            value += jdsValues[off] * fetchVectorX<ValueType,useTexture>( x_d, i );
                        }
                        off += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
                    }
                }

                result_d[k] = value;
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gevm_kernel_alpha_one(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType beta,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numColumns,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType k = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( k < numColumns )
        {
            if( ndlg > 0 )
            {
                ValueType value = 0.0;
                ValueType summand = 0.0;

                if ( beta != 0.0 )
                {
                    summand = beta * y_d[k];
                }

                int nonEmptyRows = jdsDLG[0];
                for ( int ii = 0; ii < nonEmptyRows; ++ii ) // through the rows of the first diagonal
                {
                    IndexType i = jdsPerm[ii]; // row in matrix

                    int off = ii;// position in jdsJA, jdsValues

                    int ni = jdsILG[ii];// number entries in row

                    for ( int jj = 0; jj < ni; ++jj )
                    {
                        IndexType j = jdsJA[off];
                        if ( j == k )
                        {
                            value += jdsValues[off] * fetchVectorX<ValueType,useTexture>( x_d, i );
                        }
                        off += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
                    }
                }

                result_d[k] = value + summand;
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gevm_kernel_beta_one(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType alpha,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numColumns,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType k = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( k < numColumns )
        {
            if( ndlg > 0 )
            {
                ValueType value = 0.0;
                ValueType summand = y_d[k];

                int nonEmptyRows = jdsDLG[0];
                for ( int ii = 0; ii < nonEmptyRows; ++ii ) // through the rows of the first diagonal
                {
                    IndexType i = jdsPerm[ii]; // row in matrix

                    int off = ii;// position in jdsJA, jdsValues

                    int ni = jdsILG[ii];// number entries in row

                    for ( int jj = 0; jj < ni; ++jj )
                    {
                        IndexType j = jdsJA[off];
                        if ( j == k )
                        {
                            value += jdsValues[off] * fetchVectorX<ValueType,useTexture>( x_d, i );
                        }
                        off += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
                    }
                }

                result_d[k] = alpha * value + summand;
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType,bool useTexture,bool useSharedMem>
    __global__
    void normal_gevm_kernel_beta_zero(
                    ValueType* const result_d,
                    const ValueType* x_d,
                    const ValueType* y_d,
                    const ValueType alpha,
                    const ValueType* const jdsValues,
                    const IndexType* const jdsDLG,
                    const IndexType* const jdsILG,
                    const IndexType* jdsJA,
                    const IndexType* jdsPerm,
                    IndexType numColumns,
                    const IndexType ndlg )
    {
        extern __shared__ IndexType dlg[];
        const IndexType k = threadId( gridDim, blockIdx, blockDim, threadIdx );

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

        if ( k < numColumns )
        {
            if( ndlg > 0 )
            {
                ValueType value = 0.0;

                int nonEmptyRows = jdsDLG[0];
                for ( int ii = 0; ii < nonEmptyRows; ++ii ) // through the rows of the first diagonal
                {
                    IndexType i = jdsPerm[ii]; // row in matrix

                    int off = ii;// position in jdsJA, jdsValues

                    int ni = jdsILG[ii];// number entries in row

                    for ( int jj = 0; jj < ni; ++jj )
                    {
                        IndexType j = jdsJA[off];
                        if ( j == k )
                        {
                            value += jdsValues[off] * fetchVectorX<ValueType,useTexture>( x_d, i );
                        }
                        off += fetch_JDSdlg<useTexture,useSharedMem>( jdsDLG, dlg, jj );
                    }
                }

                result_d[k] = alpha * value;
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType>
    void CUDAJDSUtils::normalGEVM(
                    ValueType result[],
                    const ValueType alpha,
                    const ValueType x[],
                    const ValueType beta,
                    const ValueType y[],
                    const IndexType numColumns,
                    const IndexType jdsPerm[],
                    const IndexType jdsILG[],
                    const IndexType ndlg,
                    const IndexType jdsDLG[],
                    const IndexType jdsJA[],
                    const ValueType jdsValues[],
                    SyncToken* syncToken )
    {
        if ( ( beta == Constants<ValueType>::one ) && ( result == y ) )
        {
            // result = alpha * A * x + beta * y ->  result += alpha * A * x

            sparseGEVM( result, alpha, x, numColumns, jdsPerm, jdsILG, ndlg, jdsDLG, jdsJA, jdsValues, syncToken );

            return;
        }

        SCAI_REGION( "CUDA.JDS.normalGEVM" )

        SCAI_LOG_INFO( logger, "normalGEVM<" << getScalarType<ValueType>() << ">"
                        << " result[ " << numColumns << "] = " << alpha
                        << " * A( #jds_diags = " << ndlg << " ) * x + " << beta << " * y " )

        SCAI_LOG_DEBUG( logger, "x = " << x << ", y = " << y << ", result = " << result )

        const bool useTexture = CUDASettings::useTexture();
        const bool useSharedMem = CUDASettings::useSharedMem();
        const int blockSize = CUDASettings::getBlockSize( numColumns );

        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( numColumns, dimBlock.x );

        SCAI_CHECK_CUDA_ACCESS

        cudaStream_t stream = 0; // default stream if no SyncToken is available

        SCAI_LOG_INFO( logger, "Start normal_gevm_kernel<" << getScalarType<ValueType>()
                        << "> <<< blockSize = " << blockSize << ", stream = " << stream
                        << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
        }

        int sharedMemSize = 0;
        if ( useSharedMem )
        {
            sharedMemSize = ndlg * sizeof(int);
        }

        if ( useTexture )
        {
            vectorBindTexture( x );

            if ( useSharedMem )
            {
                if( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_alpha_one_beta_one<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::zero )
                {
                    normal_gevm_kernel_alpha_one_beta_zero<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero && beta == Constants<ValueType>::one )
                {
                    assign_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, jdsPerm, numColumns );
                }
                else if ( alpha == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_alpha_one<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero )
                {
                    scale_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, beta, jdsPerm, numColumns );
                }
                else if ( beta == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_beta_one<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( beta == Constants<ValueType>::zero )
                {
                    normal_gevm_kernel_beta_zero<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else
                {
                    normal_gevm_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
            }
            else // no sharedMem
            {
                vectorBindTexture( jdsDLG );

                if( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_alpha_one_beta_one<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::zero )
                {
                    normal_gevm_kernel_alpha_one_beta_zero<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero && beta == Constants<ValueType>::one )
                {
                    assign_kernel<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, jdsPerm, numColumns );
                }
                else if ( alpha == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_alpha_one<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero )
                {
                    scale_kernel<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, beta, jdsPerm, numColumns );
                }
                else if ( beta == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_beta_one<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( beta == Constants<ValueType>::zero )
                {
                    normal_gevm_kernel_beta_zero<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else
                {
                    normal_gevm_kernel<ValueType, true, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
            }

            if ( !syncToken )
            {
                // synchronize here and unbind texture

                SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: gevmKernel FAILED" )

                vectorUnbindTexture( x );

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

                syncToken->pushRoutine( common::bind( unbindV, x ) );

                if ( !useSharedMem )
                {
                    syncToken->pushRoutine( common::bind( unbindI, jdsDLG ) );
                }
            }
        }
        else // no Texture cache
        {
            if ( useSharedMem )
            {
                if( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_alpha_one_beta_one<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::zero )
                {
                    normal_gevm_kernel_alpha_one_beta_zero<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero && beta == Constants<ValueType>::one )
                {
                    assign_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, jdsPerm, numColumns );
                }
                else if ( alpha == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_alpha_one<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero )
                {
                    scale_kernel<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, beta, jdsPerm, numColumns );
                }
                else if ( beta == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_beta_one<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( beta == Constants<ValueType>::zero )
                {
                    normal_gevm_kernel_beta_zero<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else
                {
                    normal_gevm_kernel<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
            }
            else // no sharedMem
            {
                if( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_alpha_one_beta_one<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::one && beta == Constants<ValueType>::zero )
                {
                    normal_gevm_kernel_alpha_one_beta_zero<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero && beta == Constants<ValueType>::one )
                {
                    assign_kernel<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, jdsPerm, numColumns );
                }
                else if ( alpha == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_alpha_one<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( alpha == Constants<ValueType>::zero )
                {
                    scale_kernel<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, y, beta, jdsPerm, numColumns );
                }
                else if ( beta == Constants<ValueType>::one )
                {
                    normal_gevm_kernel_beta_one<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else if ( beta == Constants<ValueType>::zero )
                {
                    normal_gevm_kernel_beta_zero<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
                else
                {
                    normal_gevm_kernel<ValueType, false, false><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                    ( result, x, y, alpha, beta, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numColumns, ndlg );
                }
            }
            if ( !syncToken )
            {
                SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: normal_gevm_kernel FAILED" )
            }
        }
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
                    SyncToken* syncToken )
    {
        SCAI_REGION( "CUDA.JDS.sparseGEMV" )

        SCAI_LOG_INFO( logger, "sparseGEMV<" << getScalarType<ValueType>() << ">"
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

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
        }

        SCAI_LOG_INFO( logger, "Start jdsgemvSparseKernel<" << getScalarType<ValueType>()
                        << "> <<< blockSize = " << blockSize << ", stream = " << stream
                        << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

        if ( useTexture )
        {
            vectorBindTexture( x );

            if ( useSharedMem )
            {
                const int sharedMemSize = ndlg * sizeof(int);

                sparse_gemv_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
            }
            else // no sharedMem
            {
                vectorBindTexture( jdsDLG );

                sparse_gemv_kernel<ValueType, true, false><<<dimGrid, dimBlock, 0, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
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
                    syncToken->pushRoutine( common::bind( unbindI, jdsDLG ) );
                }

                syncToken->pushRoutine( common::bind( unbindV, x ) );
            }
        }
        else // no use of Texture cache
        {
            if ( useSharedMem )
            {
                const int sharedMemSize = ndlg * sizeof( int );

                cudaFuncSetCacheConfig( sparse_gemv_kernel<ValueType, false, true>, cudaFuncCachePreferL1 );

                sparse_gemv_kernel<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
            }
            else // no use of sharedMem
            {
                cudaFuncSetCacheConfig( sparse_gemv_kernel<ValueType, false, false>, cudaFuncCachePreferL1 );

                sparse_gemv_kernel<ValueType, false, false><<<dimGrid, dimBlock, 0, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
            }

            if ( !syncToken )
            {
                SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: gemvSparseKernel FAILED" )
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    template<typename ValueType>
    void CUDAJDSUtils::sparseGEVM(
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
                    SyncToken* syncToken )
    {
        SCAI_REGION( "CUDA.JDS.sparseGEVM" )

        SCAI_LOG_INFO( logger, "sparseGEVM<" << getScalarType<ValueType>() << ">"
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

        if ( syncToken )
        {
            CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
            SCAI_ASSERT_DEBUG( cudaStreamSyncToken, "no cuda stream sync token provided" )
            stream = cudaStreamSyncToken->getCUDAStream();
        }

        SCAI_LOG_INFO( logger, "Start jdsgevMSparseKernel<" << getScalarType<ValueType>()
                        << "> <<< blockSize = " << blockSize << ", stream = " << stream
                        << ", useTexture = " << useTexture << ", useSharedMem = " << useSharedMem << ">>>" );

        if ( useTexture )
        {
            vectorBindTexture( x );

            if ( useSharedMem )
            {
                const int sharedMemSize = ndlg * sizeof(int);

                sparse_gevm_kernel<ValueType, true, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
            }
            else // no sharedMem
            {
                vectorBindTexture( jdsDLG );

                sparse_gevm_kernel<ValueType, true, false><<<dimGrid, dimBlock, 0, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
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
                    syncToken->pushRoutine( common::bind( unbindI, jdsDLG ) );
                }

                syncToken->pushRoutine( common::bind( unbindV, x ) );
            }
        }
        else // no use of Texture cache
        {
            if ( useSharedMem )
            {
                const int sharedMemSize = ndlg * sizeof( int );

                cudaFuncSetCacheConfig( sparse_gevm_kernel<ValueType, false, true>, cudaFuncCachePreferL1 );

                sparse_gevm_kernel<ValueType, false, true><<<dimGrid, dimBlock, sharedMemSize, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
            }
            else // no use of sharedMem
            {
                cudaFuncSetCacheConfig( sparse_gevm_kernel<ValueType, false, false>, cudaFuncCachePreferL1 );

                sparse_gevm_kernel<ValueType, false, false><<<dimGrid, dimBlock, 0, stream>>>
                ( result, x, alpha, jdsValues, jdsDLG, jdsILG, jdsJA, jdsPerm, numRows, ndlg );
            }

            if ( !syncToken )
            {
                SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "JDS: gemvSparseKernel FAILED" )
            }
        }
    }

    /* --------------------------------------------------------------------------- */

    void CUDAJDSUtils::setInterface( JDSUtilsInterface& JDSUtils )
    {
        SCAI_LOG_INFO( logger, "set JDS routines for CUDA in Interface" )

        LAMA_INTERFACE_REGISTER( JDSUtils, sortRows )
        LAMA_INTERFACE_REGISTER( JDSUtils, checkDiagonalProperty )
        LAMA_INTERFACE_REGISTER( JDSUtils, ilg2dlg )
        LAMA_INTERFACE_REGISTER( JDSUtils, setInversePerm )

#define LAMA_JDS_UTILS2_REGISTER(z, J, TYPE )                                             \
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, TYPE, ARITHMETIC_CUDA_TYPE_##J )        \
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, TYPE, ARITHMETIC_CUDA_TYPE_##J )      \
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, TYPE, ARITHMETIC_CUDA_TYPE_##J )    \
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, TYPE, ARITHMETIC_CUDA_TYPE_##J )  \
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, TYPE, ARITHMETIC_CUDA_TYPE_##J )  \

#define LAMA_JDS_UTILS_REGISTER(z, I, _)                                                  \
    LAMA_INTERFACE_REGISTER_T( JDSUtils, normalGEMV, ARITHMETIC_CUDA_TYPE_##I )           \
    LAMA_INTERFACE_REGISTER_T( JDSUtils, normalGEVM, ARITHMETIC_CUDA_TYPE_##I )           \
    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobi, ARITHMETIC_CUDA_TYPE_##I )               \
    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobiHalo, ARITHMETIC_CUDA_TYPE_##I )           \
                                                                                          \
    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT,                                            \
                     LAMA_JDS_UTILS2_REGISTER,                                            \
                     ARITHMETIC_CUDA_TYPE_##I )                                           \

    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT, LAMA_JDS_UTILS_REGISTER, _ )

#undef LAMA_JDS_UTILS_REGISTER
#undef LAMA_JDS_UTILS2_REGISTER

    }

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDAJDSUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( context::CUDA );
    setInterface( interface.JDSUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDAJDSUtils::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
