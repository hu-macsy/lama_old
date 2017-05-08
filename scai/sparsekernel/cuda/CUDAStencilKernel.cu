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

#include <scai/common/cuda/CUDATexVector.hpp>
#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/cuda/CUDAError.hpp>
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

/** Inline predicate to check if a stencil point is still in the 2D grid */

__inline__ __device__
bool isInner2( const IndexType pos0, const IndexType pos1, 
               const int offset[2], const IndexType size[2] )
{
    if ( !isInner( pos0, offset[0], size[0] ) ) return false;
    if ( !isInner( pos1, offset[1], size[1] ) ) return false;
    return true;
}

/** Inline predicate to check if a stencil point is still in the 3D grid */

__inline__ __device__
bool isInner3( const IndexType pos0, const IndexType pos1, const IndexType pos2, 
                             const int offset[3], const IndexType size[3] )
{
    if ( !isInner( pos0, offset[0], size[0] ) ) return false;
    if ( !isInner( pos1, offset[1], size[1] ) ) return false;
    if ( !isInner( pos2, offset[2], size[2] ) ) return false;
    return true;
}

/** Inline predicate to check if a 4-dimensional stencil point is still in the grid */

__inline__ __device__
bool isInner4( const IndexType pos0, const IndexType pos1, 
               const IndexType pos2, const IndexType pos3, 
               const int offset[4], const IndexType size[4] )
{
    if ( !isInner( pos0, offset[0], size[0] ) ) return false;
    if ( !isInner( pos1, offset[1], size[1] ) ) return false;
    if ( !isInner( pos2, offset[2], size[2] ) ) return false;
    if ( !isInner( pos3, offset[3], size[3] ) ) return false;
    return true;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv1InnerKernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{    
    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx ) + i0;

    if ( i >= i1 )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0];

    ValueType v = 0;

    for ( IndexType p = 0; p < nPoints; ++p )
    {
        v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
    }

    result[ gridPos] += alpha * v;
}

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV1Inner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV1Inner" )

    SCAI_LOG_INFO( logger,  "stencilGEMV1Inner<" << common::TypeTraits<ValueType>::id() 
                             << " on " << gridBounds[0] << " - " << gridBounds[1] )


    const IndexType n1 = gridBounds[1] - gridBounds[0];

    const int blockSize = common::CUDASettings::getBlockSize( n1 );
        
    dim3 dimBlock( blockSize, 1, 1 );
    
    dim3 dimGrid = makeGrid( n1, dimBlock.x );
    
    // allocate arrays on device for stencil data
    
    thrust::device_vector<IndexType> dGridBounds(2);
    thrust::device_vector<IndexType> dGridDistances(1);

    // copy gridSizes, gridBounds, gridDistances, stencilNodes, stencilVal, stencilOffset to GPU */

    thrust::copy( gridBounds, gridBounds + 2, dGridBounds.begin() );
    thrust::copy( gridDistances, gridDistances + 1, dGridDistances.begin() );

    gemv1InnerKernel<<< dimGrid, dimBlock>>>( result, alpha, x,
                                              dGridBounds.data().get(), dGridDistances.data().get(),
                                              nPoints, stencilVal, stencilOffset );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv1InnerKernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv2InnerKernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[2],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{    
    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];

    const IndexType j = j0 + ( blockIdx.x * blockDim.x ) + threadIdx.x;
    const IndexType i = i0 + ( blockIdx.y * blockDim.y ) + threadIdx.y;


    if ( i >= i1 || j >= j1 )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

    ValueType v = 0;

    for ( IndexType p = 0; p < nPoints; ++p )
    {
        v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
    }

    result[ gridPos] += alpha * v;
}

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV2Inner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV2Inner" )

    SCAI_LOG_INFO( logger,  "stencilGEMV2Inner<" << common::TypeTraits<ValueType>::id() 
                             << " on " << gridBounds[0] << " - " << gridBounds[1] 
                             << " x " << gridBounds[2] << " - " << gridBounds[3] )


    // allocate arrays on device for stencil data
    
    thrust::device_vector<IndexType> dGridBounds(4);
    thrust::device_vector<IndexType> dGridDistances(2);

    // copy gridSizes, gridBounds, gridDistances, stencilNodes, stencilVal, stencilOffset to GPU */

    thrust::copy( gridBounds, gridBounds + 4, dGridBounds.begin() );
    thrust::copy( gridDistances, gridDistances + 2, dGridDistances.begin() );

    IndexType n0 = gridBounds[1] - gridBounds[0];
    IndexType n1 = gridBounds[3] - gridBounds[2];

    dim3 threadsPerBlock( 16, 16 );

    // Note: gridSizes do not have to be multiple of threads in one block dimension

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x, 
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    gemv2InnerKernel<<< numBlocks, threadsPerBlock>>>( result, alpha, x,
                                                       dGridBounds.data().get(), dGridDistances.data().get(),
                                                       nPoints, stencilVal, stencilOffset );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv2InnerKernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv3InnerKernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{    
    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];
    
    const IndexType k = k0 + ( blockIdx.x * blockDim.x ) + threadIdx.x;
    const IndexType j = j0 + ( blockIdx.y * blockDim.y ) + threadIdx.y;
    const IndexType i = i0 + blockIdx.z;

    if ( i >= i1 || j >= j1  || k >= k1 )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2];

    ValueType v = 0;

    for ( IndexType p = 0; p < nPoints; ++p )
    {
        v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
    }

    result[ gridPos] += alpha * v;
}

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV3Inner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV3Inner" )

    SCAI_LOG_INFO( logger,  "stencilGEMV3Inner<" << common::TypeTraits<ValueType>::id() 
                             << " on " << gridBounds[0] << " - " << gridBounds[1] 
                             << " x " << gridBounds[2] << " - " << gridBounds[3] 
                             << " x " << gridBounds[4] << " - " << gridBounds[5] )

    // allocate arrays on device for stencil data
    
    thrust::device_vector<IndexType> dGridBounds(6);
    thrust::device_vector<IndexType> dGridDistances(3);

    // copy gridSizes, gridBounds, gridDistances, stencilNodes, stencilVal, stencilOffset to GPU */

    thrust::copy( gridBounds, gridBounds + 6, dGridBounds.begin() );
    thrust::copy( gridDistances, gridDistances + 3, dGridDistances.begin() );

    dim3 threadsPerBlock( 16, 16 );

    // Note: gridSizes do not have to be multiple of threads in one block dimension

    IndexType n0 = gridBounds[1] - gridBounds[0];
    IndexType n1 = gridBounds[3] - gridBounds[2];
    IndexType n2 = gridBounds[5] - gridBounds[4];

    dim3 numBlocks( ( n2 + threadsPerBlock.x - 1 ) / threadsPerBlock.x, 
                    ( n1 + threadsPerBlock.y - 1 ) / threadsPerBlock.y, n0 );

    gemv3InnerKernel<<< numBlocks, threadsPerBlock>>>( result, alpha, x,
                                                       dGridBounds.data().get(), dGridDistances.data().get(),
                                                       nPoints, stencilVal, stencilOffset );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv3InnerKernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv4InnerKernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{    
    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];
    const IndexType m0 = gridBounds[6];
    const IndexType m1 = gridBounds[7];
    
    const IndexType n1 = j1 - j0;

    const IndexType m = m0 + ( blockIdx.x * blockDim.x ) + threadIdx.x;
    const IndexType k = k0 + ( blockIdx.y * blockDim.y ) + threadIdx.y;
    IndexType i = blockIdx.z / n1 ;
    const IndexType j = blockIdx.z - i * n1 + j0;
    i += i0;

    if ( i >= i1 || j >= j1  || k >= k1  || m >= m1 )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];

    ValueType v = 0;

    for ( IndexType p = 0; p < nPoints; ++p )
    {
        v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
    }

    result[gridPos] += alpha * v;
}

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV4Inner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV4Inner" )

    SCAI_LOG_INFO( logger,  "stencilGEMV4Inner<" << common::TypeTraits<ValueType>::id() 
                             << " on " << gridBounds[0] << " - " << gridBounds[1] 
                             << " x " << gridBounds[2] << " - " << gridBounds[3] 
                             << " x " << gridBounds[4] << " - " << gridBounds[5] 
                             << " x " << gridBounds[6] << " - " << gridBounds[7] )

    // allocate arrays on device for stencil data
    
    thrust::device_vector<IndexType> dGridBounds(8);
    thrust::device_vector<IndexType> dGridDistances(4);

    // copy gridSizes, gridBounds, gridDistances, stencilNodes, stencilVal, stencilOffset to GPU */

    thrust::copy( gridBounds, gridBounds + 8, dGridBounds.begin() );
    thrust::copy( gridDistances, gridDistances + 4, dGridDistances.begin() );

    dim3 threadsPerBlock( 16, 16 );

    // Note: gridSizes do not have to be multiple of threads in one block dimension

    IndexType n0 = gridBounds[1] - gridBounds[0];
    IndexType n1 = gridBounds[3] - gridBounds[2];
    IndexType n2 = gridBounds[5] - gridBounds[4];
    IndexType n3 = gridBounds[7] - gridBounds[6];

    dim3 numBlocks( ( n3 + threadsPerBlock.x - 1 ) / threadsPerBlock.x, 
                    ( n2 + threadsPerBlock.y - 1 ) / threadsPerBlock.y, n1 * n0 );

    gemv4InnerKernel<<< numBlocks, threadsPerBlock>>>( result, alpha, x,
                                                       dGridBounds.data().get(), dGridDistances.data().get(),
                                                       nPoints, stencilVal, stencilOffset );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv4InnerKernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAStencilKernel::stencilGEMVInner(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims,
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    switch ( nDims ) 
    {
          case 1 : stencilGEMV1Inner( result, alpha, x, gridBounds, gridDistances,
                                      nPoints, stencilVal, stencilOffset );
                   break;

          case 2 : stencilGEMV2Inner( result, alpha, x, gridBounds, gridDistances,
                                      nPoints, stencilVal, stencilOffset );
                   break;

          case 3 : stencilGEMV3Inner( result, alpha, x, gridBounds, gridDistances,
                                      nPoints, stencilVal, stencilOffset );
                   break;

          case 4 : stencilGEMV4Inner( result, alpha, x, gridBounds, gridDistances,
                                      nPoints, stencilVal, stencilOffset );
                   break;

        default: COMMON_THROWEXCEPTION( "stencilGEMVInner for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv1BorderKernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];

    if ( i < i0 || i >= i1 )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0];

    ValueType v = 0;

    for ( IndexType p = 0; p < nPoints; ++p )
    {
        if ( isInner( i, stencilNodes[ p ], gridSizes [0] ) )
        {
            v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
        }
    }

    result[ gridPos] += alpha * v;
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV1Border(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV1Border" )

    SCAI_LOG_INFO( logger,  "stencilGEMV1Border<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << gridBounds[0] << " - " << gridBounds[1] )

    const IndexType n1 = gridSizes[0];

    const int blockSize = common::CUDASettings::getBlockSize( n1 );

    dim3 dimBlock( blockSize, 1, 1 );

    dim3 dimGrid = makeGrid( n1, dimBlock.x );

    // allocate arrays on device for stencil data

    thrust::device_vector<IndexType> dGridSizes(1);     
    thrust::device_vector<IndexType> dGridBounds(2);     
    thrust::device_vector<IndexType> dGridDistances(1);     

    // copy gridSizes, gridBounds, gridDistances, stencilNodes, stencilVal, stencilOffset to GPU */

    thrust::copy( gridSizes, gridSizes + 1, dGridSizes.begin() );
    thrust::copy( gridBounds, gridBounds + 2, dGridBounds.begin() );
    thrust::copy( gridDistances, gridDistances + 1, dGridDistances.begin() );

    gemv1BorderKernel<<< dimGrid, dimBlock>>>( result, alpha, x, 
                                               dGridSizes.data().get(), dGridBounds.data().get(),
                                               dGridDistances.data().get(),
                                               nPoints, stencilNodes, stencilVal, stencilOffset );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv1BorderKernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv2BorderKernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];

    const IndexType j = j0 + ( blockIdx.x * blockDim.x ) + threadIdx.x;
    const IndexType i = i0 + ( blockIdx.y * blockDim.y ) + threadIdx.y;

    if ( i >= i1 || j >= j1 ) 
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1];

    ValueType v = 0;

    for ( IndexType p = 0; p < nPoints; ++p )
    {
        if ( isInner2( i, j, &stencilNodes[ 2 * p ], gridSizes ) )
        {
            v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
        }
    }

    result[ gridPos] += alpha * v;
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV2Border(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV2Border" )

    SCAI_LOG_INFO( logger,  "stencilGEMV2Border<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << gridBounds[0] << " - " << gridBounds[1] 
                             << " x " << gridBounds[2] << " - " << gridBounds[3] )

    // allocate arrays on device for stencil data

    thrust::device_vector<IndexType> dGridSizes(2);     
    thrust::device_vector<IndexType> dGridBounds(4);     
    thrust::device_vector<IndexType> dGridDistances(2);     

    // copy gridSizes, gridBounds, gridDistances, stencilNodes, stencilVal, stencilOffset to GPU */

    thrust::copy( gridSizes, gridSizes + 2, dGridSizes.begin() );
    thrust::copy( gridBounds, gridBounds + 4, dGridBounds.begin() );
    thrust::copy( gridDistances, gridDistances + 2, dGridDistances.begin() );

    IndexType* dg = dGridSizes.data().get();

    // Note: gridSizes do not have to be multiple of threads in one block dimension

    IndexType n0 = gridBounds[1] - gridBounds[0];
    IndexType n1 = gridBounds[3] - gridBounds[2];

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    gemv2BorderKernel<<< numBlocks, threadsPerBlock>>>( result, alpha, x, 
                                                        dGridSizes.data().get(), dGridBounds.data().get(),
                                                        dGridDistances.data().get(),
                                                        nPoints, stencilNodes, stencilVal, stencilOffset );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv2BorderKernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv3BorderKernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];

    const IndexType k = k0 + ( blockIdx.x * blockDim.x ) + threadIdx.x;
    const IndexType j = j0 + ( blockIdx.y * blockDim.y ) + threadIdx.y;
    const IndexType i = i0 + ( blockIdx.z * blockDim.z ) + threadIdx.z;

    if ( i >= i1 || j >= j1  || k >= k1 )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2];

    ValueType v = 0;

    for ( IndexType p = 0; p < nPoints; ++p )
    {
        if ( isInner3( i, j, k, &stencilNodes[ 3 * p ], gridSizes ) )
        {
            v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
        }
    }

    result[ gridPos] += alpha * v;
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV3Border(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV3Border" )

    SCAI_LOG_INFO( logger,  "stencilGEMV3Border<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << gridBounds[0] << " - " << gridBounds[1] 
                             << " x " << gridBounds[2] << " - " << gridBounds[3] 
                             << " x " << gridBounds[4] << " - " << gridBounds[5] )

    // allocate arrays on device for stencil data

    thrust::device_vector<IndexType> dGridSizes(3);     
    thrust::device_vector<IndexType> dGridBounds(6);     
    thrust::device_vector<IndexType> dGridDistances(3);     

    // copy gridSizes, gridBounds, gridDistances, stencilNodes, stencilVal, stencilOffset to GPU */

    thrust::copy( gridSizes, gridSizes + 3, dGridSizes.begin() );
    thrust::copy( gridBounds, gridBounds + 6, dGridBounds.begin() );
    thrust::copy( gridDistances, gridDistances + 3, dGridDistances.begin() );


    // Note: gridSizes do not have to be multiple of threads in one block dimension

    IndexType n0 = gridBounds[1] - gridBounds[0];
    IndexType n1 = gridBounds[3] - gridBounds[2];
    IndexType n2 = gridBounds[5] - gridBounds[4];

    // ToDo: usually one of the values n0, n1 or n2 is very small and we could shape
    //       threadsPerBlock correspondingly.
  
    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n2 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n1 + threadsPerBlock.y - 1 ) / threadsPerBlock.y, 
                    ( n0 + threadsPerBlock.z - 1 ) / threadsPerBlock.z  );

    gemv3BorderKernel<<< numBlocks, threadsPerBlock>>>( result, alpha, x, 
                                                        dGridSizes.data().get(), dGridBounds.data().get(),
                                                        dGridDistances.data().get(),
                                                        nPoints, stencilNodes, stencilVal, stencilOffset );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv3BorderKernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void gemv4BorderKernel(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    const IndexType i0 = gridBounds[0];
    const IndexType i1 = gridBounds[1];
    const IndexType j0 = gridBounds[2];
    const IndexType j1 = gridBounds[3];
    const IndexType k0 = gridBounds[4];
    const IndexType k1 = gridBounds[5];
    const IndexType m0 = gridBounds[6];
    const IndexType m1 = gridBounds[7];

    const IndexType n1 = j1 - j0;

    const IndexType m  = m0 + ( blockIdx.x * blockDim.x ) + threadIdx.x;
    const IndexType k  = k0 + ( blockIdx.y * blockDim.y ) + threadIdx.y;

    const IndexType ij = ( blockIdx.z * blockDim.z ) + threadIdx.z;

    IndexType i = ij / n1 ;

    const IndexType j = j0 + ij - i * n1;

    i += i0;

    if ( i >= i1 || j >= j1  || k >= k1  || m >= m1 )
    {
        return;
    }

    IndexType gridPos = i * gridDistances[0] + j * gridDistances[1] + k * gridDistances[2] + m * gridDistances[3];

    ValueType v = 0;

    for ( IndexType p = 0; p < nPoints; ++p )
    {
        if ( isInner4( i, j, k, m, &stencilNodes[ 4 * p ], gridSizes ) )
        {
            v += stencilVal[p] * x[ gridPos + stencilOffset[p] ];
        }
    }

    result[ gridPos] += alpha * v;
}   

template<typename ValueType>
void CUDAStencilKernel::stencilGEMV4Border(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_REGION( "CUDA.Stencil.GEMV4Border" )

    SCAI_LOG_INFO( logger,  "stencilGEMV4Border<" << common::TypeTraits<ValueType>::id() << "> on " 
                             << gridBounds[0] << " - " << gridBounds[1] 
                             << " x " << gridBounds[2] << " - " << gridBounds[3] 
                             << " x " << gridBounds[4] << " - " << gridBounds[5] 
                             << " x " << gridBounds[6] << " - " << gridBounds[7] )

    // allocate arrays on device for stencil data

    thrust::device_vector<IndexType> dGridSizes(4);     
    thrust::device_vector<IndexType> dGridBounds(8);     
    thrust::device_vector<IndexType> dGridDistances(4);     

    // copy gridSizes, gridBounds, gridDistances, stencilNodes, stencilVal, stencilOffset to GPU */

    thrust::copy( gridSizes, gridSizes + 4, dGridSizes.begin() );
    thrust::copy( gridBounds, gridBounds + 8, dGridBounds.begin() );
    thrust::copy( gridDistances, gridDistances + 4, dGridDistances.begin() );

    IndexType n0 = gridBounds[1] - gridBounds[0];
    IndexType n1 = gridBounds[3] - gridBounds[2];
    IndexType n2 = gridBounds[5] - gridBounds[4];
    IndexType n3 = gridBounds[7] - gridBounds[6];

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n3 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n2 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n0 * n1 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    gemv4BorderKernel<<< numBlocks, threadsPerBlock>>>( result, alpha, x, 
                                               dGridSizes.data().get(), dGridBounds.data().get(),
                                               dGridDistances.data().get(),
                                               nPoints, stencilNodes, stencilVal, stencilOffset );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "gemv4BorderKernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAStencilKernel::stencilGEMVBorder(
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType gridBounds[],
    const IndexType gridDistances[],
    const IndexType nPoints,
    const int stencilNodes[],
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    switch ( nDims ) 
    {
        case 1 : stencilGEMV1Border( result, alpha, x, gridSizes, gridBounds, gridDistances,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 2 : stencilGEMV2Border( result, alpha, x, gridSizes, gridBounds, gridDistances,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 3 : stencilGEMV3Border( result, alpha, x, gridSizes, gridBounds, gridDistances,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        case 4 : stencilGEMV4Border( result, alpha, x, gridSizes, gridBounds, gridDistances,
                                     nPoints, stencilNodes, stencilVal, stencilOffset );
                 break;

        default: COMMON_THROWEXCEPTION( "stencilGEMVBorder for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAStencilKernel::stencilGEMVCaller(
    IndexType gridBounds[],
    ValueType result[], 
    const ValueType alpha,  
    const ValueType x[],
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType lb[],
    const IndexType ub[],
    const IndexType gridDistances[],
    const IndexType currentDim,
    const IndexType nPoints,
    const int stencilNodes[], 
    const ValueType stencilVal[],
    const int stencilOffset[] )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( currentDim, nDims, "illegal current dimension" )

    if ( gridSizes[currentDim] <= lb[currentDim] + ub[currentDim] )
    {
        // no inner part in current dimension, call boundary routine, afterwards all is done

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
 
        return;
    }

    // if left boundary, handle it

    if ( lb[currentDim] > 0 )
    {
        gridBounds[2 * currentDim] = 0;
        gridBounds[2 * currentDim + 1] = lb[currentDim];

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
    }

    // set inner boundaries,

    gridBounds[2 * currentDim] = lb[currentDim];
    gridBounds[2 * currentDim + 1] = gridSizes[currentDim] - ub[currentDim];

    if ( currentDim + 1 == nDims )
    {
        // all boundaries are now inner ones, we can call the routine for the inner points

        stencilGEMVInner( result, alpha, x, nDims, gridBounds, gridDistances, nPoints, stencilVal, stencilOffset );
    }
    else
    { 
        // recursive call to set grid bounds for the next dimension 

        stencilGEMVCaller( gridBounds, result, alpha, x, nDims, gridSizes, lb, ub, gridDistances, currentDim + 1, 
                           nPoints, stencilNodes, stencilVal, stencilOffset );
    }

    // if right boundary, travere it

    if ( ub[currentDim] > 0 )
    {
        gridBounds[2 * currentDim] = gridSizes[currentDim] - ub[currentDim];
        gridBounds[2 * currentDim + 1] = gridSizes[currentDim];

        stencilGEMVBorder( result, alpha, x, nDims, gridSizes, gridBounds, gridDistances, nPoints, stencilNodes, stencilVal, stencilOffset );
    }

    // reset the boundaries

    gridBounds[2 * currentDim] = 0;
    gridBounds[2 * currentDim + 1] = gridSizes[currentDim];
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
    // prepare array with gridBounds for the recursive traverser routine

    IndexType gridBounds[ 2 * SCAI_GRID_MAX_DIMENSION ];

    for ( IndexType i = 0; i < nDims; ++i )
    {
        gridBounds[2 * i] = 0;
        gridBounds[2 * i + 1] = gridSizes[i];
    }

    IndexType currentDim = 0;

    thrust::device_vector<int> dStencilNodes( nDims * nPoints );
    thrust::device_vector<ValueType> dStencilVal( nPoints );
    thrust::device_vector<int> dStencilOffset( nPoints );

    // copy stencil data to GPU

    thrust::copy( stencilNodes, stencilNodes + nDims * nPoints, dStencilNodes.begin() );
    thrust::copy( stencilVal, stencilVal + nPoints, dStencilVal.begin() );
    thrust::copy( stencilOffset, stencilOffset + nPoints, dStencilOffset.begin() );

    stencilGEMVCaller( gridBounds, result, alpha, x, nDims, gridSizes, lb, ub, gridDistances, currentDim,
                       nPoints, dStencilNodes.data().get(), dStencilVal.data().get(), dStencilOffset.data().get() );
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
