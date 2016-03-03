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
 */

// hpp
#include <scai/sparsekernel/cuda/CUDACOOUtils.hpp>

// local library
#include <scai/sparsekernel/COOKernelTrait.hpp>

// internal scai library
#include <scai/utilskernel/cuda/CUDAUtils.hpp>

#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/bind.hpp>

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/Container.hpp>

// thrust
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

using namespace scai::tasking;

namespace scai
{

using common::TypeTraits;
using common::CUDASettings;
using utilskernel::CUDAUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( CUDACOOUtils::logger, "CUDA.COOUtils" )

/* --------------------------------------------------------------------------- */

__device__ inline void cooAtomicAdd( double* address, double val )
{
    unsigned long long int* address_as_ull =
        ( unsigned long long int* ) address;

    unsigned long long int old = *address_as_ull, assumed;

    do
    {
        assumed = old;
        old = atomicCAS( address_as_ull, assumed,
                         __double_as_longlong( val +
                                               __longlong_as_double( assumed ) ) );
    }
    while ( assumed != old );
}

__device__ inline void cooAtomicAdd( float* address, float val )

{
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 200

    // CUDA runtime offers faster solution for capability >= 2.0

    atomicAdd( address, val );

#else

    // old slow solution

    int i_val = __float_as_int( val );

    int tmp0 = 0;

    int tmp1;

    while ( ( tmp1 = atomicCAS( ( int* ) address, tmp0, i_val ) ) != tmp0 )

    {
        tmp0 = tmp1;
        i_val = __float_as_int( val + __int_as_float( tmp1 ) );
    }

#endif
}

#ifdef SCAI_COMPLEX_SUPPORTED

__device__ inline void cooAtomicAdd( ComplexFloat* address, ComplexFloat val )
{
    cooAtomicAdd( ( float* )&address[0], val.real() );
    cooAtomicAdd( ( float* )&address[1], val.imag() );
}

__device__ inline void cooAtomicAdd( ComplexDouble* address, ComplexDouble val )
{
    cooAtomicAdd( ( double* )&address[0], val.real() );
    cooAtomicAdd( ( double* )&address[1], val.imag() );
}

#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__ void cooGemvKernel(
    ValueType* result,
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

        const ValueType resultUpdate = alpha * cooValues[k] * fetchVectorX<ValueType, useTexture>( x, j );

        // atomic add required, solution above

        cooAtomicAdd( &result[i], resultUpdate );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__ void cooGemvKernel_alpha_one(
    ValueType* result,
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

        const ValueType resultUpdate = cooValues[k] * fetchVectorX<ValueType, useTexture>( x, j );

        // atomic add required, solution above

        cooAtomicAdd( &result[i], resultUpdate );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__ void cooGevmKernel(
    ValueType* result,
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

        const ValueType resultUpdate = alpha * cooValues[k] * fetchVectorX<ValueType, useTexture>( x, i );

        // atomic add required, solution above

        cooAtomicAdd( &result[j], resultUpdate );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__ void cooGevmKernel_alpha_one(
    ValueType* result,
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

        const ValueType resultUpdate = cooValues[k] * fetchVectorX<ValueType, useTexture>( x, i );

        // atomic add required, solution above

        cooAtomicAdd( &result[j], resultUpdate );
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
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[] )
{
    SCAI_REGION( "CUDA.COO.normalGEMV" )

    SCAI_LOG_INFO( logger, "normalGEMV<" << TypeTraits<ValueType>::id() << ">, "
                   << "result[ " << numRows << "] = " << alpha
                   << " COO( #vals = " << numValues << " ) * x + " << beta << " * y" )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
        SCAI_LOG_INFO( logger, "asyncronous execution on stream " << stream );
    }

    bool useTexture = CUDASettings::useTexture();

    IndexType blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numValues, dimBlock.x );

    // set result = beta * y, not needed if beta == 1 and y == result

    if ( beta == scai::common::constants::ONE && result == y )
    {
        SCAI_LOG_DEBUG( logger, "normalGEMV is sparse, no init of result needed" )
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "normalGEMV, set result = " << beta << " * y " )
        // setScale also deals with y undefined for beta == 0
        CUDAUtils::setScale( result, beta, y, numRows );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "COO: initGemvKernel FAILED" )

    if ( numValues == 0 )
    {
        return;
    }

    blockSize = CUDASettings::getBlockSize( numValues );
    dimBlock = dim3( blockSize, 1, 1 );
    dimGrid = makeGrid( numValues, dimBlock.x );

    SCAI_LOG_INFO( logger, "Start cooGemvKernel<" << TypeTraits<ValueType>::id()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", alpha = " << alpha
                   << ", useTexture = " << useTexture << ">>>" )

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == scai::common::constants::ONE )
        {
            cooGemvKernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock>>>
            ( result, x, numValues, cooIA, cooJA, cooValues );
        }
        else
        {
            cooGemvKernel<ValueType, true> <<< dimGrid, dimBlock>>>
            ( result, alpha, x, numValues, cooIA, cooJA, cooValues );
        }
    }
    else
    {
        if ( alpha == scai::common::constants::ONE )
        {
            cooGemvKernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock>>>
            ( result, x, numValues, cooIA, cooJA, cooValues );
        }
        else
        {
            cooGemvKernel<ValueType, false> <<< dimGrid, dimBlock>>>
            ( result, alpha, x, numValues, cooIA, cooJA, cooValues );
        }
    }

    if ( !syncToken )
    {
        // synchronization now, unbind texture if it has been used

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "COO: gemvKernel FAILED" )

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

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDACOOUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[] )
{
    SCAI_REGION( "CUDA.COO.normalGEVM" )

    SCAI_LOG_INFO( logger, "normalGEVM, #rows = " << numRows << ", #vals = " << numValues )

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;

    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
        SCAI_LOG_INFO( logger, "asyncronous execution on stream " << stream );
    }

    bool useTexture = CUDASettings::useTexture();

    IndexType blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numValues, dimBlock.x );

    // set result = beta * y, not needed if beta == 1 and y == result

    if ( beta == scai::common::constants::ONE && result == y )
    {
        SCAI_LOG_DEBUG( logger, "normalGEVM is sparse, no init of result needed" )
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "normalGEMV, set result = " << beta << " * y " )
        CUDAUtils::setScale( result, beta, y, numRows );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "COO: initGevmKernel FAILED" )

    if ( numValues == 0 )
    {
        return;
    }

    blockSize = CUDASettings::getBlockSize( numValues );
    dimBlock = dim3( blockSize, 1, 1 );
    dimGrid = makeGrid( numValues, dimBlock.x );

    SCAI_LOG_INFO( logger, "Start cooGevmKernel<" << TypeTraits<ValueType>::id()
                   << "> <<< blockSize = " << blockSize << ", stream = " << stream
                   << ", useTexture = " << useTexture << ">>>" )

    if ( useTexture )
    {
        vectorBindTexture( x );

        if ( alpha == scai::common::constants::ONE )
        {
            cooGevmKernel_alpha_one<ValueType, true> <<< dimGrid, dimBlock>>>
            ( result, x, numValues, cooIA, cooJA, cooValues );
        }
        else
        {
            cooGevmKernel<ValueType, true> <<< dimGrid, dimBlock>>>
            ( result, alpha, x, numValues, cooIA, cooJA, cooValues );
        }
    }
    else
    {
        if ( alpha == scai::common::constants::ONE )
        {
            cooGevmKernel_alpha_one<ValueType, false> <<< dimGrid, dimBlock>>>
            ( result, x, numValues, cooIA, cooJA, cooValues );
        }
        else
        {
            cooGevmKernel<ValueType, false> <<< dimGrid, dimBlock>>>
            ( result, alpha, x, numValues, cooIA, cooJA, cooValues );
        }
    }

    if ( !syncToken )
    {
        // synchronization now, unbind texture if it has been used

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "COO: gevmKernel FAILED" )

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

/* --------------------------------------------------------------------------- */

__global__
static void offsets2ia_kernel( IndexType* cooIA, const IndexType* csrIA, const IndexType numRows, const IndexType numDiagonals )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType csrOffset = csrIA[i];
        IndexType cooOffset = 0; // additional offset due to diagonals

        if ( i < numDiagonals )
        {
            // diagonal elements will be the first nrows entries

            cooIA[i] = i;
            csrOffset += 1;// do not fill diagonal element again
            cooOffset = numDiagonals - i - 1;// offset in coo moves
        }

        // now fill remaining part of row i

        for ( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
        {
            cooIA[ jj + cooOffset] = i;
        }
    }
}

/* --------------------------------------------------------------------------- */

void CUDACOOUtils::offsets2ia(
    IndexType cooIA[],
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numDiagonals )
{
    SCAI_LOG_INFO( logger,
                   "build cooIA( " << numValues << " ) from csrIA( " << ( numRows + 1 )
                   << " ), #diagonals = " << numDiagonals )

    SCAI_CHECK_CUDA_ACCESS

    // make grid

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    offsets2ia_kernel <<< dimGrid, dimBlock>>>( cooIA, csrIA, numRows, numDiagonals );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "sync for offsets2ia_kernel" )
}

/* --------------------------------------------------------------------------- */

__global__
static void build_offset_kernel(
    IndexType* offsets,
    const IndexType n,
    const IndexType* ia,
    const IndexType nz )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    // Entries in offset filled every time there is a change in values of consecutive elements
    //   i:     0  1  2  3  4  5
    //  ia:     0  0  1  1  1  3
    // nd1:     0  0  1  1  1  3
    // nd2:     0  1  1  1  3  4
    //             x        x  x
    //             |        |  |->                6
    //             |        |---->          5  5
    //             |------------->       2
    // offset:                        0  2  5  5  6

    if ( i < nz )
    {
        IndexType nd1 = ia[i];
        IndexType nd2 = n;

        if ( i + 1 < nz )
        {
            nd2 = ia[i + 1];
        }

        for ( IndexType j = nd1; j < nd2; j++ )
        {
            offsets[j + 1] = i + 1;
        }

        if ( i == 0 )
        {
            for ( IndexType i = 0; i <= nd1; i++ )
            {
                offsets[i] = 0;
            }
        }
    }
}

__global__
static void add_diagonals_kernel(
    IndexType* offsets,
    const IndexType numRows,
    const IndexType numDiagonals )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    //  #diagonals = 3               |
    //  offsets( in ) :  0  3  4   7   9  10
    //  offsets( out ):  0  4  6  10  12  13
    //  i                   0  1   2   3   4

    if ( i < numRows )
    {
        if ( i < numDiagonals )
        {
            offsets[i + 1] += i + 1;
        }
        else
        {
            offsets[i + 1] += numDiagonals;
        }
    }
}

/* --------------------------------------------------------------------------- */

void CUDACOOUtils::ia2offsets(
    IndexType csrIA[],
    const IndexType numRows,
    const IndexType numDiagonals,
    const IndexType cooIA[],
    const IndexType numValues )
{
    SCAI_LOG_INFO( logger,
                   "build csrIA( " << numRows + 1 << " ) from cooIA( " << ( numValues )
                   << " ), #diagonals = " << numDiagonals )

    // Note: the array cooIA is assumed to be sorted after the diagonal elements

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0;// default stream, asynchronous execution not supported here

    const int blockSize = CUDASettings::getBlockSize();
    const dim3 dimBlock( blockSize, 1, 1 );
    const dim3 dimGrid = makeGrid( numValues, dimBlock.x );

    build_offset_kernel <<< dimGrid, dimBlock>>>( csrIA, numRows,
            cooIA + numDiagonals, numValues - numDiagonals );

    // increment offsets for the diagonal elements

    if ( numDiagonals > 0 )
    {
        const dim3 dimGrid = makeGrid( numRows, dimBlock.x );

        add_diagonals_kernel <<< dimGrid, dimBlock>>>( csrIA, numRows, numDiagonals );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "normalGEMV, stream = " << stream )
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType, typename CSRValueType>
__global__
static void csr2coo_kernel( COOValueType* cooValues, const CSRValueType* csrValues,
                            const IndexType* csrIA, const IndexType numRows, const IndexType numDiagonals )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType csrOffset = csrIA[i];
        IndexType cooOffset = 0; // additional offset due to diagonals

        if ( i < numDiagonals )
        {
            // diagonal elements will be the first nrows entries

            cooValues[i] = csrValues[csrOffset];

            csrOffset += 1;// do not fill diagonal element again
            cooOffset = numDiagonals - i - 1;// offset in coo moves
        }

        // now fill remaining part of row i

        for ( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
        {
            cooValues[ jj + cooOffset] = static_cast<COOValueType>( csrValues[ jj ] );
        }
    }
}

template<typename COOValueType, typename CSRValueType>
void CUDACOOUtils::setCSRData(
    COOValueType cooValues[],
    const CSRValueType csrValues[],
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numDiagonals )
{
    SCAI_LOG_INFO( logger,
                   "build cooValues( << " << numValues << " from csrValues + csrIA( " << ( numRows + 1 )
                   << " ), #diagonals = " << numDiagonals )

    SCAI_CHECK_CUDA_ACCESS

    // make grid

    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );

    csr2coo_kernel <<< dimGrid, dimBlock>>>( cooValues, csrValues, csrIA, numRows, numDiagonals );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "sync for csr2coo_kernel" )
}

/* --------------------------------------------------------------------------- */

void CUDACOOUtils::Registrator::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::CUDA;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register COOUtils CUDA-routines for CUDA at kernel registry [" << flag << "]" )

    KernelRegistry::set<COOKernelTrait::offsets2ia>( CUDACOOUtils::offsets2ia, CUDA, flag );
    KernelRegistry::set<COOKernelTrait::setCSRData<IndexType, IndexType> >( CUDACOOUtils::setCSRData, CUDA, flag );
}

template<typename ValueType>
void CUDACOOUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::CUDA;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register COOUtils CUDA-routines for CUDA at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<COOKernelTrait::normalGEMV<ValueType> >( CUDACOOUtils::normalGEMV, CUDA, flag );
    KernelRegistry::set<COOKernelTrait::normalGEVM<ValueType> >( CUDACOOUtils::normalGEVM, CUDA, flag );
}

template<typename ValueType, typename OtherValueType>
void CUDACOOUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::CUDA;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register COOUtils CUDA-routines for CUDA at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<COOKernelTrait::setCSRData<ValueType, OtherValueType> >( CUDACOOUtils::setCSRData, CUDA, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDACOOUtils::CUDACOOUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_CUDA_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, ARITHMETIC_CUDA_LIST, ARITHMETIC_CUDA_LIST>::call( flag );
}

CUDACOOUtils::~CUDACOOUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_CUDA> ValueTypes;
    typedef common::mepr::ContainerVO<RegistratorVO, ARITHMETIC_CUDA> MoreValueTypes;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_CUDA_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, ARITHMETIC_CUDA_LIST, ARITHMETIC_CUDA_LIST>::call( flag );
}

CUDACOOUtils CUDACOOUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
