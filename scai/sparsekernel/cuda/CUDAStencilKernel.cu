/**
 * @file sparsekernel/cuda/CUDAStencilKernel.cu
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
 * @brief CUDA Implementations on GPU for stencil kernels.
 * @author Thomas Brandes
 * @date 04.05.2017
 */

#include <scai/hmemo/HArray.hpp>
#include <scai/tracing.hpp>

#include <scai/common/Grid.hpp>
#include <scai/sparsekernel/cuda/CUDAStencilKernel.hpp>
#include <scai/sparsekernel/StencilKernelTrait.hpp>
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/common/Settings.hpp>

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAUtils.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/Grid.hpp>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>

#include <functional>

namespace scai
{

using common::Grid;

namespace sparsekernel
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( CUDAStencilKernel::logger, "CUDA.StencilKernel" )

/* --------------------------------------------------------------------------- */

/** Help routine to determine the left position in a dimension with a certain boUnaryOp type. */

__inline__ __device__ 
bool getBorderPosL( IndexType& pos, const IndexType offset, const IndexType size, const IndexType border )
{
    bool valid = true;

    if ( pos >= offset )
    {
        pos = pos - offset;  // is a valid pos
    }
    else if ( border == 0 )
    {
        valid = false;
    }
    else if ( border == 1 )
    {
        pos = ( pos + size ) - offset;
    }

    return valid;
}

/** Help routine to determine the right position in a dimension with a certain boUnaryOp type. */

__inline__ __device__ 
bool getBorderPosR( IndexType& pos, const IndexType offset, const IndexType size, const IndexType border )
{
    bool valid = true;

    if ( pos + offset < size )
    {
        pos += offset;  // is a valid pos
    }
    else if ( border == Grid::BORDER_ABSORBING )
    {
        valid = false;
    }
    else if ( border == Grid::BORDER_PERIODIC )
    {
        pos = ( pos + offset ) - size; 
    }
    return valid;
}

/** Help routine to determine the correct stencil position depending on border types.
 *
 *  Note: gridSizes and gridBorders are available in constant memory.
 */

__inline__ __device__ 
bool getOffsetPos( 
    IndexType pos[],
    const int stencilPositions[],
    const IndexType p,
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType gridBorders[] )
{
    bool valid = true;

    for ( IndexType iDim = 0; iDim < nDims; ++iDim )
    {
        int offs = stencilPositions[ nDims * p + iDim ];

        if ( offs < 0 )
        {
            valid = getBorderPosL( pos[iDim], static_cast<IndexType>( -offs ), gridSizes[iDim],  gridBorders[2 * iDim] );
        }
        else if ( offs > 0 )
        {
            valid = getBorderPosR( pos[iDim], static_cast<IndexType>( offs ), gridSizes[iDim], gridBorders[ 2 * iDim + 1] );
        }
        if ( !valid )
        {
            break;
        }
    }

    return valid;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv1Kernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[] )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= gridSizes[0] )
    {
        return;   // might happen if gridSizes[0] is not multiple of blockDim.x
    }

    IndexType gridPos = i * gridDistances[0];

    ValueType v = 0;

    if ( ( i >= gridStencilWidth[0] ) && ( i < gridSizes[0] - gridStencilWidth[1] ) )
    {
        // gridPoint ( i ) is inner point, we have not to check for valid stencil points

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            v += stencilVal[ p ] * x[ gridPos + stencilOffset[ p ] ];
        }
    }
    else 
    {
        // gridPoint ( i, j) is border point, we have to check each stencil neighbor individually

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            IndexType pos[] = { i };

            bool valid = getOffsetPos( pos, stencilPositions, p, 1, gridSizes, gridBorders );

            if ( !valid )
            {
                continue;
            }

            IndexType stencilLinearPos = pos[0] * gridDistances[0];

            v += stencilVal[ p ] * x[ stencilLinearPos ];
        }
    }

    if ( beta == 0 )
    {
        result[ gridPos] = alpha * v;
    }
    else
    {
        result[ gridPos] = alpha * v + beta * y[ gridPos ];
    }
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV1(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType hostGridSizes[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV1" )

    IndexType n0 = hostGridSizes[0];

    SCAI_LOG_INFO( logger,  "stencilGEMV1<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " grid" )

    dim3 threadsPerBlock( 256 );

    dim3 numBlocks( ( n0 + threadsPerBlock.x - 1 ) / threadsPerBlock.x );

    gemv1Kernel<ValueType><<< numBlocks, threadsPerBlock>>>( 
        result, alpha, x, beta, y, nPoints, stencilPositions, stencilVal, stencilOffset,
        gridSizes, gridDistances, gridBorders, gridStencilWidth );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv2Kernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[] )
{
    const IndexType j = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i = blockIdx.y * blockDim.y + threadIdx.y;

    if ( j >= gridSizes[1]  || i >= gridSizes[0] )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

    ValueType v = 0;

    if (    ( i >= gridStencilWidth[0] ) && ( i < gridSizes[0] - gridStencilWidth[1] ) 
         && ( j >= gridStencilWidth[2] ) && ( j < gridSizes[1] - gridStencilWidth[3] ) )
    {
        // gridPoint(i,j) is inner point, all stencil points can be applied

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            v += stencilVal[ p ] * x[ gridPos + stencilOffset[ p ] ];
        }
    }
    else 
    {
        // gridPoint(i,j) is border point, each stencil neighbor is checked individually

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            IndexType pos[] = { i, j };

            bool valid = getOffsetPos( pos, stencilPositions, p, 2, gridSizes, gridBorders );

            if ( !valid )
            {
                continue;
            }

            IndexType stencilLinearPos = pos[0] * gridDistances[0] + pos[1] * gridDistances[1];

            v += stencilVal[ p ] * x[ stencilLinearPos ];
        }
    }

    if ( beta == 0 )
    {
        result[ gridPos] = alpha * v;
    }
    else
    {
        result[ gridPos] = alpha * v + beta * y[ gridPos ];
    }
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV2(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType hostGridSizes[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV2" )

    IndexType n0 = hostGridSizes[0];
    IndexType n1 = hostGridSizes[1];

    SCAI_LOG_INFO( logger,  "stencilGEMV2<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << n0 << " x " << n1 )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    gemv2Kernel<ValueType><<< numBlocks, threadsPerBlock>>>( 
        result, alpha, x, beta, y, nPoints, stencilPositions, stencilVal, stencilOffset,
        gridSizes, gridDistances, gridBorders, gridStencilWidth );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv3Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv3Kernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[] )
{
    __shared__ IndexType smGridInfo[18];

    __shared__ int smStencilOffset[64];
    __shared__ ValueType smStencilVal[64];

    IndexType* smGridSize = smGridInfo;
    IndexType* smGridDistance = smGridInfo + 3;
    IndexType* smGridBorders  = smGridInfo + 6;
    IndexType* smGridStencilWidth = smGridInfo + 12;

    IndexType tid = threadIdx.x;

    if ( tid < 3 )
    {
        smGridSize[tid] = gridSizes[tid];
        smGridDistance[tid] = gridDistances[tid];
    }

    if ( tid < 6 )
    {
        smGridBorders[tid] = gridBorders[tid];
        smGridStencilWidth[tid] = gridStencilWidth[tid];
    }

    if ( tid < nPoints )
    {
        smStencilOffset[tid] = stencilOffset[tid];
        smStencilVal[tid] = stencilVal[tid];
    }

    __syncthreads();

    const IndexType k = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType j = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i = blockIdx.z * blockDim.z + threadIdx.z;

    if ( j >= smGridSize[1]  || k >= smGridSize[2] || i >= smGridSize[0] )
    {
        return;
    }

    IndexType gridPos = i * smGridDistance[0] + j * smGridDistance[1] + k * smGridDistance[2];

    ValueType v = 0;

    if (    ( i >= smGridStencilWidth[0] ) && ( i < smGridSize[0] - smGridStencilWidth[1] ) 
         && ( j >= smGridStencilWidth[2] ) && ( j < smGridSize[1] - smGridStencilWidth[3] ) 
         && ( k >= smGridStencilWidth[4] ) && ( k < smGridSize[2] - smGridStencilWidth[5] ) )
    {
        // gridPoint ( i, j, k ) is inner point, we have not to check for valid stencil points

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            v += smStencilVal[ p ] * x[ gridPos + smStencilOffset[ p ] ];
        }
    }
    else 
    {
        // gridPoint ( i, j, k ) is border point, we have to check each stencil neighbor individually

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            IndexType pos[] = { i, j, k };

            bool valid = getOffsetPos( pos, stencilPositions, p, 3, smGridSize, smGridBorders );

            if ( !valid )
            {
                continue;
            }

            IndexType stencilLinearPos =   pos[0] * smGridDistance[0] + pos[1] * smGridDistance[1] 
                                         + pos[2] * smGridDistance[2];

            v += stencilVal[ p ] * x[ stencilLinearPos ];
        }
    }

    if ( beta == 0 )
    {
        result[ gridPos] = alpha * v;
    }
    else
    {
        result[ gridPos] = alpha * v + beta * y[ gridPos ];
    }
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV3(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType hostGridSizes[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV3" )

    SCAI_CHECK_CUDA_ACCESS

    IndexType n0 = hostGridSizes[0];
    IndexType n1 = hostGridSizes[1];
    IndexType n2 = hostGridSizes[2];

    SCAI_LOG_INFO( logger,  "stencilGEMV3<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << n0 << " x " << n1 << " x " << n2 )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n2 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n1 + threadsPerBlock.y - 1 ) / threadsPerBlock.y, 
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.z );

    cudaStream_t stream = 0; // default stream if no syncToken is given

    tasking::CUDAStreamSyncToken* syncToken = tasking::CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // asynchronous execution takes other stream and will not synchronize later
        stream = syncToken->getCUDAStream();
    }

    gemv3Kernel<ValueType><<< numBlocks, threadsPerBlock, 0, stream>>>( 
        result, alpha, x, beta, y, nPoints, stencilPositions, stencilVal, stencilOffset,
        gridSizes, gridDistances, gridBorders, gridStencilWidth );

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv3Kernel failed" ) ;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv4Kernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[] )
{
    const IndexType m = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType k = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType ij = ( blockIdx.z * blockDim.z ) + threadIdx.z;
    const IndexType i  = ij / gridSizes[1];
    const IndexType j  = ij - i * gridSizes[1];

    // as grid sizes are not always multiple of the corresponding threads, check for valid grid point

    if ( m >= gridSizes[3]  || k >= gridSizes[2] || j >= gridSizes[1] || i >= gridSizes[0] )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];

    ValueType v = 0;

    if (    ( i >= gridStencilWidth[0] ) && ( i < gridSizes[0] - gridStencilWidth[1] ) 
         && ( j >= gridStencilWidth[2] ) && ( j < gridSizes[1] - gridStencilWidth[3] ) 
         && ( k >= gridStencilWidth[4] ) && ( k < gridSizes[2] - gridStencilWidth[5] ) 
         && ( m >= gridStencilWidth[6] ) && ( m < gridSizes[3] - gridStencilWidth[7] ) )
    {
        // gridPoint(i, j, k, m) is inner point, we have not to check for valid stencil points

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            v += stencilVal[ p ] * x[ gridPos + stencilOffset[ p ] ];
        }
    }
    else 
    {
        // gridPoint ( i, j, k ) is border point, we have to check each stencil neighbor individually

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            IndexType pos[] = { i, j, k, m };

            bool valid = getOffsetPos( pos, stencilPositions, p, 4, gridSizes, gridBorders );

            if ( !valid )
            {
                continue;
            }

            IndexType stencilLinearPos =   pos[0] * gridDistances[0] + pos[1] * gridDistances[1] 
                                         + pos[2] * gridDistances[2] + pos[3] * gridDistances[3];

            v += stencilVal[ p ] * x[ stencilLinearPos ];
        }
    }

    if ( beta == 0 )
    {
        result[ gridPos] = alpha * v;
    }
    else
    {
        result[ gridPos] = alpha * v + beta * y[ gridPos ];
    }
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV4(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType hostGridSizes[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV4" )

    IndexType n0 = hostGridSizes[0];
    IndexType n1 = hostGridSizes[1];
    IndexType n2 = hostGridSizes[2];
    IndexType n3 = hostGridSizes[3];

    SCAI_LOG_INFO( logger, "stencilGEMV4<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << n0 << " x " << n1 << " x " << n2 << " x " << n3 )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n3 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n2 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n1 * n0 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    gemv4Kernel<ValueType><<< numBlocks, threadsPerBlock>>>( 
        result, alpha, x, beta, y, nPoints, stencilPositions, stencilVal, stencilOffset,
        gridSizes, gridDistances, gridBorders, gridStencilWidth );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv4Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAStencilKernel::normalGEMV( 
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const ValueType beta,  
    const ValueType y[],
    const IndexType nDims, 
    const IndexType hostGridSizes[],
    const IndexType gridSizes[],
    const IndexType gridDistances[],
    const IndexType gridBorders[],
    const IndexType gridStencilWidth[],
    const IndexType nPoints,
    const int stencilPositions[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_LOG_INFO( logger, "normalGEMV" << nDims << ", #points = " << nPoints )

    switch ( nDims ) 
    {
        case 1 : stencilGEMV1( result, alpha, x, beta, y, hostGridSizes,
                               nPoints, stencilPositions, stencilVal, stencilOffset,
                               gridSizes, gridDistances, gridBorders, gridStencilWidth );
                 break;

        case 2 : stencilGEMV2( result, alpha, x, beta, y, hostGridSizes,
                               nPoints, stencilPositions, stencilVal, stencilOffset,
                               gridSizes, gridDistances, gridBorders, gridStencilWidth );
                 break;

        case 3 : stencilGEMV3( result, alpha, x, beta, y, hostGridSizes,
                               nPoints, stencilPositions, stencilVal, stencilOffset,
                               gridSizes, gridDistances, gridBorders, gridStencilWidth );
                 break;

        case 4 : stencilGEMV4( result, alpha, x, beta, y, hostGridSizes,
                               nPoints, stencilPositions, stencilVal, stencilOffset,
                               gridSizes, gridDistances, gridBorders, gridStencilWidth );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilGEMV for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAStencilKernel::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::ContextType ctx = common::ContextType::CUDA;

    SCAI_LOG_DEBUG( logger,
                    "register StencilKernel CUDA-routines for Host at kernel registry [" << flag 
                    << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<StencilKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
}

/* --------------------------------------------------------------------------- */

CUDAStencilKernel::CUDAStencilKernel()
{
    SCAI_LOG_INFO( logger, "register StencilKernel CUDA-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    // Registrator::registerKernels( flag );

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */

CUDAStencilKernel::~CUDAStencilKernel()
{
    SCAI_LOG_INFO( logger, "unregister StencilKernel CUDA-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    // Registrator::registerKernels( flag );

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */

CUDAStencilKernel CUDAStencilKernel::guard;

/* --------------------------------------------------------------------------- */

}

}
