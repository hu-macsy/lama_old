/**
 * @file sparsekernel/cuda/CUDACOOUtils.cu
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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

#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAUtils.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>

// thrust
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

#include <functional>

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
        common::CUDAUtils::atomicAdd( &result[i], resultUpdate );
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
        common::CUDAUtils::atomicAdd( &result[i], resultUpdate );
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
        common::CUDAUtils::atomicAdd( &result[j], resultUpdate );
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
        common::CUDAUtils::atomicAdd( &result[j], resultUpdate );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
static inline void launchGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],  
    common::MatrixOp op,
    cudaStream_t stream )
{
    IndexType blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numValues, dimBlock.x );

    if ( alpha == common::Constants::ONE )
    {
        if ( common::isTranspose( op ) )
        {
            cooGevmKernel_alpha_one<ValueType, useTexture> <<< dimGrid, dimBlock, 0, stream>>>
            ( result, x, numValues, cooIA, cooJA, cooValues );
        }
        else
        {
            cooGemvKernel_alpha_one<ValueType, useTexture> <<< dimGrid, dimBlock, 0, stream>>>
            ( result, x, numValues, cooIA, cooJA, cooValues );
        }
    }
    else
    {
        if ( common::isTranspose( op ) )
        {
            cooGevmKernel<ValueType, useTexture> <<< dimGrid, dimBlock, 0, stream>>>
            ( result, alpha, x, numValues, cooIA, cooJA, cooValues );
        }
        else
        {
            cooGemvKernel<ValueType, useTexture> <<< dimGrid, dimBlock, 0, stream>>>
            ( result, alpha, x, numValues, cooIA, cooJA, cooValues );
        }
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
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const common::MatrixOp op )
{
    SCAI_REGION( "CUDA.COO.normalGEMV" )

    const IndexType nResult = common::isTranspose( op ) ? numColumns : numRows;

    SCAI_LOG_INFO( logger, "normalGEMV<" << TypeTraits<ValueType>::id() << ">, "
                   << "result[ " << nResult << "] = " << alpha
                   << " COO( #vals = " << numValues << " ) * x + " << beta << " * y" )

    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
        SCAI_LOG_INFO( logger, "asyncronous execution on stream " << stream );
    }

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

    if ( numValues == 0 )
    {
        return;
    }

    bool useTexture = CUDASettings::useTexture();

    SCAI_LOG_INFO( logger, "Start cooGemvKernel<" << TypeTraits<ValueType>::id()
                   << "> <<< stream = " << stream
                   << ", alpha = " << alpha
                   << ", useTexture = " << useTexture << ">>>" )

    if ( useTexture )
    {
        vectorBindTexture( x );

        launchGEMV<ValueType, true>( result, alpha, x, numValues, cooIA, cooJA, cooValues, op, stream );
    }
    else
    {
        launchGEMV<ValueType, false>( result, alpha, x, numValues, cooIA, cooJA, cooValues, op, stream );
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
            syncToken->pushRoutine( std::bind( unbind, x ) );
        }
    }
}

/* --------------------------------------------------------------------------- */

__global__
static void offsets2ia_kernel( IndexType* cooIA, const IndexType* csrIA, const IndexType numRows )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numRows )
    {
        IndexType csrOffset = csrIA[i];

        for ( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
        {
            cooIA[ jj ] = i;
        }
    }
}

/* --------------------------------------------------------------------------- */

void CUDACOOUtils::offsets2ia(
    IndexType cooIA[],
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType numRows )
{
    SCAI_REGION( "CUDA.COO.offsets2ia" )

    SCAI_LOG_INFO( logger,
                   "build cooIA( " << numValues << " ) from csrIA( " << ( numRows + 1 ) << " )" )
    SCAI_CHECK_CUDA_ACCESS
    // make grid
    const int blockSize = CUDASettings::getBlockSize();
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( numRows, dimBlock.x );
    offsets2ia_kernel <<< dimGrid, dimBlock>>>( cooIA, csrIA, numRows );
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

/* --------------------------------------------------------------------------- */

void CUDACOOUtils::ia2offsets(
    IndexType csrIA[],
    const IndexType numRows,
    const IndexType cooIA[],
    const IndexType numValues )
{
    SCAI_REGION( "CUDA.COO.ia2offsets" )

    SCAI_LOG_INFO( logger,
                   "build csrIA( " << numRows + 1 << " ) from cooIA( " << ( numValues ) << " )" )

    // Note: the array cooIA is assumed to be sorted

    SCAI_CHECK_CUDA_ACCESS
    cudaStream_t stream = 0;// default stream, asynchronous execution not supported here

    const int blockSize = CUDASettings::getBlockSize();
    const dim3 dimBlock( blockSize, 1, 1 );
    const dim3 dimGrid = makeGrid( numValues, dimBlock.x );
    build_offset_kernel <<< dimGrid, dimBlock>>>( csrIA, numRows, cooIA, numValues );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "ia2offsets, stream = " << stream )
}

/* --------------------------------------------------------------------------- */

void CUDACOOUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_INFO( logger, "register COOUtils CUDA-routines for CUDA at kernel registry [" << flag << "]" )
    KernelRegistry::set<COOKernelTrait::offsets2ia>( CUDACOOUtils::offsets2ia, ctx, flag );
}

template<typename ValueType>
void CUDACOOUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register COOUtils CUDA-routines for CUDA at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<COOKernelTrait::normalGEMV<ValueType> >( CUDACOOUtils::normalGEMV, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDACOOUtils::CUDACOOUtils()
{
    SCAI_LOG_INFO( logger, "register COOUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDACOOUtils::~CUDACOOUtils()
{
    SCAI_LOG_INFO( logger, "unregister COOUtilsKernel CUDA version" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDACOOUtils CUDACOOUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
