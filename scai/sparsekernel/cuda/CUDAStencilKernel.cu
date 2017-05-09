/**
 * @file CUDAStencilKernel.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief CUDA Implementations on GPU for stencil kernels.
 * @author Thomas Brandes
 * @date May 04, 2017
 */

#include <scai/hmemo/HArray.hpp>
#include <scai/tracing.hpp>

#include <scai/common/Grid.hpp>
#include <scai/sparsekernel/cuda/CUDAStencilKernel.hpp>
#include <scai/sparsekernel/StencilKernelTrait.hpp>

#include <scai/common/Settings.hpp>

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAUtils.hpp>
#include <scai/common/cuda/launchHelper.hpp>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>

namespace scai
{
namespace sparsekernel
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( CUDAStencilKernel::logger, "CUDA.StencilKernel" )

/* --------------------------------------------------------------------------- */

__constant__ IndexType gridDistancesD[SCAI_GRID_MAX_DIMENSION];
__constant__ IndexType gridSizesD[SCAI_GRID_MAX_DIMENSION];
__constant__ IndexType gridLB[SCAI_GRID_MAX_DIMENSION];
__constant__ IndexType gridUB[SCAI_GRID_MAX_DIMENSION];


/** This routine checks whether a point pos is an inner point 
 *
 *  @param[in] pos is the position, in range 0 .. size-1
 *  @param[in] offset specifies the relative offset to pos, might be positive or negative
 *  @param[in] size is the size of the range, only needed for offset > 0
 */
#

__inline__ __device__ 
bool isInner( const IndexType pos, const int offset, const IndexType size )
{
    if ( offset == 0 )
    {
        return true;
    }
    if ( offset < 0 )
    {
        return pos >= static_cast<IndexType>( -offset );
    }

    return pos + static_cast<IndexType>( offset ) < size;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void gemv1Kernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= gridSizesD[0] )
    {
        return;   // might happen if gridSizesD[0] is not multiple of blockDim.x
    }

    IndexType gridPos = i * gridDistancesD[0];

    ValueType v = 0;

    if ( ( i >= gridLB[0] ) && ( i < gridSizesD[0] - gridUB[0] ) )
    {
        // gridPoint ( i ) is inner point, we have not to check for valid stenil points

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            v += fetchVectorX<ValueType, useTexture>( stencilVal, p )
                 * x[ gridPos + fetchVectorX<int, useTexture>( stencilOffset, p + nPoints ) ];
        }
    }
    else 
    {
        // gridPoint ( i, j) is border point, we have to check each stencil neighbor individually

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( !isInner( i, fetchVectorX<int, useTexture>( stencilOffset, p ), gridSizesD[0] ) ) 
            {
                continue;
            }
   
            v += fetchVectorX<ValueType, useTexture>( stencilVal, p )
                 * x[ gridPos + fetchVectorX<int, useTexture>( stencilOffset, p + nPoints ) ];
        }
    }

    result[ gridPos] += alpha * v;
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV1(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV2" )

    IndexType n0 = gridSizes[0];

    SCAI_LOG_INFO( logger,  "stencilGEMV1<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " grid" )

    dim3 threadsPerBlock( 256 );

    dim3 numBlocks( ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    bool useTexture = true;  // defaut is true, might be disabled explicitly for test

    common::Settings::getEnvironment( useTexture, "SCAI_CUDA_USE_TEXTURE" );

    if ( useTexture )
    {
        vectorBindTexture( stencilOffset );
        vectorBindTexture( stencilVal );

        gemv1Kernel<ValueType, true><<< numBlocks, threadsPerBlock>>>( 
            result, alpha, x, nPoints, stencilVal, stencilOffset );

        vectorUnbindTexture( stencilVal );
        vectorUnbindTexture( stencilOffset );
    }
    else
    {
        gemv1Kernel<ValueType, false><<< numBlocks, threadsPerBlock>>>( 
            result, alpha, x, nPoints, stencilVal, stencilOffset );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void gemv2Kernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    const IndexType j = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i = blockIdx.y * blockDim.y + threadIdx.y;

    if ( j >= gridSizesD[1]  || i >= gridSizesD[0] )
    {
        return;
    }

    IndexType gridPos = i * gridDistancesD[0] + j * gridDistancesD[1];

    ValueType v = 0;

    if (    ( i >= gridLB[0] ) && ( i < gridSizesD[0] - gridUB[0] ) 
         && ( j >= gridLB[1] ) && ( j < gridSizesD[1] - gridUB[1] ) )
    {
        // gridPoint(i,j) is inner point, all stencil points can be applied

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            v += fetchVectorX<ValueType, useTexture>( stencilVal, p )
                 * x[ gridPos + fetchVectorX<int, useTexture>( stencilOffset, p + 2 * nPoints ) ];
        }
    }
    else 
    {
        // gridPoint(i,j) is border point, each stencil neighbor is checked individually

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( !isInner( i, fetchVectorX<int, useTexture>( stencilOffset, 2*p ), gridSizesD[0] ) ) 
            {
                continue;
            }
            if ( !isInner( j, fetchVectorX<int, useTexture>( stencilOffset, 2*p+1 ), gridSizesD[1] ) ) 
            {
                continue;
            }
   
            v += fetchVectorX<ValueType, useTexture>( stencilVal, p )
                 * x[ gridPos + fetchVectorX<int, useTexture>( stencilOffset, p + 2 * nPoints ) ];
        }
    }

    result[ gridPos] += alpha * v;
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV2(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV2" )

    IndexType n0 = gridSizes[0];
    IndexType n1 = gridSizes[1];

    SCAI_LOG_INFO( logger,  "stencilGEMV2<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << n0 << " x " << n1 )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    bool useTexture = true;  // defaut is true, might be disabled explicitly for test

    common::Settings::getEnvironment( useTexture, "SCAI_CUDA_USE_TEXTURE" );

    if ( useTexture )
    {
        vectorBindTexture( stencilOffset );
        vectorBindTexture( stencilVal );

        gemv2Kernel<ValueType, true><<< numBlocks, threadsPerBlock>>>( 
            result, alpha, x, nPoints, stencilVal, stencilOffset );

        vectorUnbindTexture( stencilVal );
        vectorUnbindTexture( stencilOffset );
    }
    else
    {
        gemv2Kernel<ValueType, false><<< numBlocks, threadsPerBlock>>>( 
            result, alpha, x, nPoints, stencilVal, stencilOffset );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv3Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void gemv3Kernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    const IndexType k = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType j = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i = blockIdx.z;

    if ( j >= gridSizesD[1]  || k >= gridSizesD[2] )
    {
        return;
    }

    IndexType gridPos = i * gridDistancesD[0] + j * gridDistancesD[1] + k * gridDistancesD[2];

    ValueType v = 0;

    if (    ( i >= gridLB[0] ) && ( i < gridSizesD[0] - gridUB[0] ) 
         && ( j >= gridLB[1] ) && ( j < gridSizesD[1] - gridUB[1] ) 
         && ( k >= gridLB[2] ) && ( k < gridSizesD[2] - gridUB[2] ) )
    {
        // gridPoint ( i, j, k ) is inner point, we have not to check for valid stenil points

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            v += fetchVectorX<ValueType, useTexture>( stencilVal, p )
                 * x[ gridPos + fetchVectorX<int, useTexture>( stencilOffset, p + 3 * nPoints ) ];
        }
    }
    else 
    {
        // gridPoint ( i, j, k ) is border point, we have to check each stencil neighbor individually

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( !isInner( i, fetchVectorX<int, useTexture>( stencilOffset, 3*p ), gridSizesD[0] ) ) 
            {
                continue;
            }
            if ( !isInner( j, fetchVectorX<int, useTexture>( stencilOffset, 3*p+1 ), gridSizesD[1] ) ) 
            {
                continue;
            }
            if ( !isInner( k, fetchVectorX<int, useTexture>( stencilOffset, 3*p+2 ), gridSizesD[2] ) ) 
            {
                continue;
            }
   
            v += fetchVectorX<ValueType, useTexture>( stencilVal, p )
                 * x[ gridPos + fetchVectorX<int, useTexture>( stencilOffset, p + 3 * nPoints ) ];
        }
    }

    result[ gridPos] += alpha * v;
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV3(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV3" )

    IndexType n0 = gridSizes[0];
    IndexType n1 = gridSizes[1];
    IndexType n2 = gridSizes[2];

    SCAI_LOG_INFO( logger,  "stencilGEMV3<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << n0 << " x " << n1 << " x " << n2 )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n2 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n1 + threadsPerBlock.y - 1 ) / threadsPerBlock.y, 
                      n0 );

    bool useTexture = true;  // defaut is true, might be disabled explicitly for test

    common::Settings::getEnvironment( useTexture, "SCAI_CUDA_USE_TEXTURE" );

    if ( useTexture )
    {
        vectorBindTexture( stencilOffset );
        vectorBindTexture( stencilVal );

        gemv3Kernel<ValueType, true><<< numBlocks, threadsPerBlock>>>( 
            result, alpha, x, nPoints, stencilVal, stencilOffset );

        vectorUnbindTexture( stencilVal );
        vectorUnbindTexture( stencilOffset );
    }
    else
    {
        gemv3Kernel<ValueType, false><<< numBlocks, threadsPerBlock>>>( 
            result, alpha, x, nPoints, stencilVal, stencilOffset );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv3Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, bool useTexture>
__global__
void gemv4Kernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    const IndexType m = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType k = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType ij = ( blockIdx.z * blockDim.z ) + threadIdx.z;
    const IndexType i  = ij / gridSizesD[1] ;
    const IndexType j  = ij - i * gridSizesD[1];

    if ( j >= gridSizesD[1]  || k >= gridSizesD[2] )
    {
        return;
    }

    IndexType gridPos = i * gridDistancesD[0] + j * gridDistancesD[1] + k * gridDistancesD[2];

    ValueType v = 0;

    if (    ( i >= gridLB[0] ) && ( i < gridSizesD[0] - gridUB[0] ) 
         && ( j >= gridLB[1] ) && ( j < gridSizesD[1] - gridUB[1] ) 
         && ( k >= gridLB[2] ) && ( k < gridSizesD[2] - gridUB[2] ) 
         && ( m >= gridLB[3] ) && ( m < gridSizesD[3] - gridUB[3] ) )
    {
        // gridPoint ( i, j, k ) is inner point, we have not to check for valid stenil points

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            v += fetchVectorX<ValueType, useTexture>( stencilVal, p )
                 * x[ gridPos + fetchVectorX<int, useTexture>( stencilOffset, p + 4 * nPoints ) ];
        }
    }
    else 
    {
        // gridPoint ( i, j, k ) is border point, we have to check each stencil neighbor individually

        for ( IndexType p = 0; p < nPoints; ++p )
        {
            if ( !isInner( i, fetchVectorX<int, useTexture>( stencilOffset, 4*p ), gridSizesD[0] ) ) 
            {
                continue;
            }
            if ( !isInner( j, fetchVectorX<int, useTexture>( stencilOffset, 4*p+1 ), gridSizesD[1] ) ) 
            {
                continue;
            }
            if ( !isInner( k, fetchVectorX<int, useTexture>( stencilOffset, 4*p+2 ), gridSizesD[2] ) ) 
            {
                continue;
            }
            if ( !isInner( m, fetchVectorX<int, useTexture>( stencilOffset, 4*p+3 ), gridSizesD[3] ) ) 
            {
                continue;
            }
   
            v += fetchVectorX<ValueType, useTexture>( stencilVal, p )
                 * x[ gridPos + fetchVectorX<int, useTexture>( stencilOffset, p + 4 * nPoints ) ];
        }
    }

    result[ gridPos] += alpha * v;
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV4(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV3" )

    IndexType n0 = gridSizes[0];
    IndexType n1 = gridSizes[1];
    IndexType n2 = gridSizes[2];
    IndexType n3 = gridSizes[3];

    SCAI_LOG_INFO( logger,  "stencilGEMV3<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << n0 << " x " << n1 << " x " << n2 << " x " << n3 )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n3 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n2 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n0 * n1 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    bool useTexture = true;  // defaut is true, might be disabled explicitly for test

    common::Settings::getEnvironment( useTexture, "SCAI_CUDA_USE_TEXTURE" );

    if ( useTexture )
    {
        vectorBindTexture( stencilOffset );
        vectorBindTexture( stencilVal );

        gemv4Kernel<ValueType, true><<< numBlocks, threadsPerBlock>>>( 
            result, alpha, x, nPoints, stencilVal, stencilOffset );

        vectorUnbindTexture( stencilVal );
        vectorUnbindTexture( stencilOffset );
    }
    else
    {
        gemv4Kernel<ValueType, false><<< numBlocks, threadsPerBlock>>>( 
            result, alpha, x, nPoints, stencilVal, stencilOffset );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv4Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV( 
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims, 
    const IndexType gridSizes[],
    const IndexType lb[],
    const IndexType ub[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    // allocate arrays for stencil Data on GPU

    thrust::device_vector<ValueType> dStencilVal( nPoints );
    thrust::device_vector<int> dStencilOffset( nPoints + nDims * nPoints );

    // copy stencil data to GPU, stencilOffsetD = [ stencilNodes, stencilOffset ]

    thrust::copy( stencilNodes, stencilNodes + nDims * nPoints, dStencilOffset.begin() );
    thrust::copy( stencilVal, stencilVal + nPoints, dStencilVal.begin() );
    thrust::copy( stencilOffset, stencilOffset + nPoints, dStencilOffset.begin() + nDims * nPoints );

    // copy grid data to GPU, constant memory

    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( gridDistancesD, gridDistances, 3 * sizeof(IndexType), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( gridSizesD, gridSizes, 3 * sizeof(IndexType), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( gridLB, lb, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( gridUB, ub, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );


    const IndexType* dStencilOffsetPtr = dStencilOffset.data().get();
    const ValueType* dStencilValPtr = dStencilVal.data().get();

    switch ( nDims ) 
    {
        case 1 : stencilGEMV1( result, alpha, x, gridSizes,
                               nPoints, dStencilValPtr, dStencilOffsetPtr );
                 break;

        case 2 : stencilGEMV2( result, alpha, x, gridSizes,
                               nPoints, dStencilValPtr, dStencilOffsetPtr );
                 break;

        case 3 : stencilGEMV3( result, alpha, x, gridSizes,
                               nPoints, dStencilValPtr, dStencilOffsetPtr );
                 break;

        case 4 : stencilGEMV4( result, alpha, x, gridSizes,
                               nPoints, dStencilValPtr, dStencilOffsetPtr );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilGEMV for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAStencilKernel::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    common::context::ContextType ctx = common::context::CUDA;

    SCAI_LOG_DEBUG( logger,
                    "register StencilKernel CUDA-routines for Host at kernel registry [" << flag 
                    << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<StencilKernelTrait::stencilGEMV<ValueType> >( stencilGEMV, ctx, flag );
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
