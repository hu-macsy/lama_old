/**
 * @file CUDASection.cu
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
 * @brief Implementation of kernel operations for multidimension sections with CUDA
 * @author Thomas Brandes
 * @date 17.05.2017
 */

// hpp
#include <scai/utilskernel/cuda/CUDASection.hpp>

// local library
#include <scai/utilskernel/SectionKernelTrait.hpp>

#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>

#include <scai/common/Grid.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

namespace scai
{

using common::binary;
using common::unary;
using common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( CUDASection::logger, "CUDA.Section" )

/* --------------------------------------------------------------------------- */

__constant__ IndexType sectionSizesD[SCAI_GRID_MAX_DIMENSION];
__constant__ IndexType targetDistancesD[SCAI_GRID_MAX_DIMENSION];
__constant__ IndexType sourceDistancesD[SCAI_GRID_MAX_DIMENSION];

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assign0Kernel( 
    ValueType target[],
    const ValueType source[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i > 0 )
    {
        return;  // only one element
    }

    if ( swapOperands )
    {   
        target[i] = applyBinary( source[i], op, target[i] );
    }
    else
    {   
        target[i] = applyBinary( target[i], op, source[i] );
    }
}

template<typename ValueType>
void CUDASection::assign0(
    ValueType targetSection[],
    const ValueType sourceSection[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assign0" )

    dim3 threadsPerBlock( 32 );
    dim3 numBlocks( 1 );

    assign0Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assign0Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assign1Kernel( 
    ValueType targetSection[],
    const ValueType sourceSection[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = targetSection[i * targetDistancesD[0]];
    const ValueType& source = sourceSection[i * sourceDistancesD[0]];

    if ( swapOperands )
    {   
        target = applyBinary( source, op, target );
    }
    else
    {   
        target = applyBinary( target, op, source );
    }
}

template<typename ValueType>
void CUDASection::assign1(
    ValueType targetSection[],
    const ValueType sourceSection[],
    const IndexType sizes[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assign1" )

    IndexType n0 = sizes[0];

    SCAI_LOG_INFO( logger,  "assign<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " grid" )

    dim3 threadsPerBlock( 256 );
    dim3 numBlocks( ( n0 + threadsPerBlock.x - 1 ) / threadsPerBlock.x );

    assign1Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assign1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assign2Kernel( 
    ValueType targetSection[],
    const ValueType sourceSection[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType j = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i = blockIdx.y * blockDim.y + threadIdx.y;

    if ( j >= sectionSizesD[1]  || i >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i * targetDistancesD[0] + j * targetDistancesD[1]];
    const ValueType& source = sourceSection[i * sourceDistancesD[0] + j * sourceDistancesD[1]];

    if ( swapOperands )
    {   
        target = applyBinary( source, op, target );
    }
    else
    {   
        target = applyBinary( target, op, source );
    }
}

template<typename ValueType>
void CUDASection::assign2(
    ValueType targetSection[],
    const ValueType sourceSection[],
    const IndexType sizes[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assign2" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];

    SCAI_LOG_INFO( logger,  "assign<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " x " << n1 )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    assign2Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assign2Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assign3Kernel( 
    ValueType targetSection[],
    const ValueType sourceSection[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType k = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType j = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i = blockIdx.z;

    if ( j >= sectionSizesD[1]  || k >= sectionSizesD[2] )
    {
        return;
    }

    ValueType& target = targetSection[i * targetDistancesD[0] + j * targetDistancesD[1] + k * targetDistancesD[2]];
    const ValueType& source = sourceSection[i * sourceDistancesD[0] + j * sourceDistancesD[1] + k * sourceDistancesD[2]];

    if ( swapOperands )
    {   
        target = applyBinary( source, op, target );
    }
    else
    {   
        target = applyBinary( target, op, source );
    }
}

template<typename ValueType>
void CUDASection::assign3(
    ValueType targetSection[],
    const ValueType sourceSection[],
    const IndexType sizes[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assign3" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];
    IndexType n2 = sizes[2];

    SCAI_LOG_INFO( logger,  "assign<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " x " << n1 )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n2 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n1 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                      n0 );

    assign3Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assign3Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assign4Kernel( 
    ValueType targetSection[],
    const ValueType sourceSection[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType m = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType k = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType ij = ( blockIdx.z * blockDim.z ) + threadIdx.z;
    const IndexType i  = ij / sectionSizesD[1] ;
    const IndexType j  = ij - i * sectionSizesD[1];

    // i and j will always be legal, but m or might be out of range

    if ( m >= sectionSizesD[3]  || k >= sectionSizesD[2] || j >= sectionSizesD[1] || i >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i * targetDistancesD[0] + j * targetDistancesD[1] 
                                    + k * targetDistancesD[2] + m * targetDistancesD[3]];

    const ValueType& source = sourceSection[i * sourceDistancesD[0] + j * sourceDistancesD[1] 
                                          + k * sourceDistancesD[2] + m * sourceDistancesD[3]];

    if ( swapOperands )
    {   
        target = applyBinary( source, op, target );
    }
    else
    {   
        target = applyBinary( target, op, source );
    }
}

template<typename ValueType>
void CUDASection::assign4(
    ValueType targetSection[],
    const ValueType sourceSection[],
    const IndexType sizes[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assign3" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];
    IndexType n2 = sizes[2];
    IndexType n3 = sizes[3];

    SCAI_LOG_INFO( logger,  "assign<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " x " << n1 << " x " << n2 << " x " << n3 )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n3 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n2 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n0 * n1 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    assign4Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assign4Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASection::assign(
    ValueType targetSection[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType targetDistances[],
    const ValueType sourceSection[],
    const IndexType sourceDistances[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    if ( nDims == 0 )
    {
        // apply kernel on scalar arguments, no section sizes, differences are needed

        assign0( targetSection, sourceSection, op, swapOperands );
        return;
    }

    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( sectionSizesD, sizes, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( targetDistancesD, targetDistances, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( sourceDistancesD, sourceDistances, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );

    switch ( nDims )
    {
        case 1 : assign1( targetSection, sourceSection, sizes, op, swapOperands );
                 break;

        case 2 : assign2( targetSection, sourceSection, sizes, op, swapOperands );
                 break;

        case 3 : assign3( targetSection, sourceSection, sizes, op, swapOperands );
                 break;

        case 4 : assign4( targetSection, sourceSection, sizes, op, swapOperands );
                 break;

        default: COMMON_THROWEXCEPTION( "assign for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
__global__
void unaryOp1Kernel( 
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const common::unary::UnaryOp op )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    TargetValueType& target = targetSection[i * targetDistancesD[0]];
    const SourceValueType& source = sourceSection[i * sourceDistancesD[0]];

    target = static_cast<TargetValueType>( applyUnary( op, source ) );
}

template<typename TargetValueType, typename SourceValueType>
void CUDASection::unaryOp1(
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const IndexType sizes[],
    const common::unary::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.unaryOp1" )

    IndexType n0 = sizes[0];

    SCAI_LOG_INFO( logger,  "unaryOp1<" << common::TypeTraits<TargetValueType>::id() << 
                           ", " << common::TypeTraits<SourceValueType>::id() << "> on " << n0 << " grid" )

    dim3 threadsPerBlock( 256 );
    dim3 numBlocks( ( n0 + threadsPerBlock.x - 1 ) / threadsPerBlock.x );

    unaryOp1Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "unaryOp1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
__global__
void unaryOp2Kernel( 
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const common::unary::UnaryOp op )
{
    const IndexType j = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i = blockIdx.y * blockDim.y + threadIdx.y;
      
    if ( j >= sectionSizesD[1]  || i >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    TargetValueType& target = targetSection[i * targetDistancesD[0] + j * targetDistancesD[1]];
    const SourceValueType& source = sourceSection[i * sourceDistancesD[0] + j * sourceDistancesD[1]];

    target = static_cast<TargetValueType>( applyUnary( op, source ) );
}

template<typename TargetValueType, typename SourceValueType>
void CUDASection::unaryOp2(
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const IndexType sizes[],
    const common::unary::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.unaryOp2" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];

    SCAI_LOG_INFO( logger, "unaryOp2<" << common::TypeTraits<TargetValueType>::id() << 
                           ", " << common::TypeTraits<SourceValueType>::id() 
                           << "> on " << n0 << " x " << n1 << " grid" )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    unaryOp2Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "unaryOp2Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void CUDASection::unaryOp( 
    TargetValueType targetSection[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType targetDistances[],
    const SourceValueType sourceSection[],
    const IndexType sourceDistances[],
    const common::unary::UnaryOp op )
{
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( sectionSizesD, sizes, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( targetDistancesD, targetDistances, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( sourceDistancesD, sourceDistances, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );

    switch ( nDims )
    {
        case 1 : unaryOp1( targetSection, sourceSection, sizes, op );
                 break;

        case 2 : unaryOp2( targetSection, sourceSection, sizes, op );
                 break;

        default: COMMON_THROWEXCEPTION( "unaryOp for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */
/*    assignScalar0 : kernel + launcher                                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType, bool swap>
__global__
void assignScalar0Kernel( 
    ValueType target[],
    const ValueType source,
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i > 0 )
    {
        return;  // only one element
    }

    if ( swap )
    {   
        target[i] = applyBinary( source, op, target[i] );
    }
    else
    {   
        target[i] = applyBinary( target[i], op, source );
    }
}

template<typename ValueType>
void CUDASection::assignScalar0(
    ValueType targetSection[],
    const ValueType source,
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assignScalar0" )

    dim3 threadsPerBlock( 32 );
    dim3 numBlocks( 1 );

    if ( swapOperands )
    {
        assignScalar0Kernel<ValueType, true><<< numBlocks, threadsPerBlock>>>
            ( targetSection, source, op, swapOperands );
    }
    else
    {
        assignScalar0Kernel<ValueType, false><<< numBlocks, threadsPerBlock>>>
            ( targetSection, source, op, swapOperands );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assignScalar0Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */
/*    assignScalar1 : kernel + launcher                                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assignScalar1Kernel(
    ValueType targetSection[],
    const ValueType val,
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= sectionSizesD[0] )
    {   
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = targetSection[i * targetDistancesD[0]];

    if ( swapOperands )
    {
        target = applyBinary( val, op, target );
    }
    else
    {
        target = applyBinary( target, op, val );
    }
}

template<typename ValueType>
void CUDASection::assignScalar1(
    ValueType targetSection[],
    const ValueType val,
    const IndexType sizes[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assignScalar1" )

    IndexType n0 = sizes[0];

    SCAI_LOG_INFO( logger,  "assignScalar1<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " grid" )

    dim3 threadsPerBlock( 256 );
    dim3 numBlocks( ( n0 + threadsPerBlock.x - 1 ) / threadsPerBlock.x );

    assignScalar1Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, val, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assignScalar1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */
/*    assignScalar2 : kernel + launcher                                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assignScalar2Kernel( 
    ValueType targetSection[],
    const ValueType val,
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType j = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i = blockIdx.y * blockDim.y + threadIdx.y;

    if ( j >= sectionSizesD[1]  || i >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i * targetDistancesD[0] + j * targetDistancesD[1]];

    if ( swapOperands )
    {   
        target = applyBinary( val, op, target );
    }
    else
    {   
        target = applyBinary( target, op, val );
    }
}

template<typename ValueType>
void CUDASection::assignScalar2(
    ValueType targetSection[],
    const ValueType val,
    const IndexType sizes[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assignScalar2" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];

    SCAI_LOG_INFO( logger,  "assignScalar2<" << common::TypeTraits<ValueType>::id() 
                            << "> on " << n0 << " x " << n1 << " grid" )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    assignScalar2Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, val, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assignScalar1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */
/*    assignScalar3 : kernel + launcher                                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assignScalar3Kernel( 
    ValueType targetSection[],
    const ValueType val,
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType k = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType j = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i = blockIdx.z * blockDim.z + threadIdx.z;

    if ( k >= sectionSizesD[2] || j >= sectionSizesD[1]  || i >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i * targetDistancesD[0] + j * targetDistancesD[1] + k * targetDistancesD[2]];

    if ( swapOperands )
    {   
        target = applyBinary( val, op, target );
    }
    else
    {   
        target = applyBinary( target, op, val );
    }
}

template<typename ValueType>
void CUDASection::assignScalar3(
    ValueType targetSection[],
    const ValueType val,
    const IndexType sizes[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assignScalar3" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];
    IndexType n2 = sizes[2];

    SCAI_LOG_INFO( logger,  "assignScalar2<" << common::TypeTraits<ValueType>::id() 
                            << "> on " << n0 << " x " << n1 << " x " << n2 << " grid" )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n2 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n1 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n0 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    assignScalar3Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, val, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assignScalar1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */
/*    assignScalar4 : kernel + launcher                                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void assignScalar4Kernel( 
    ValueType targetSection[],
    const ValueType val,
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    const IndexType m = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType k = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType ij = ( blockIdx.z * blockDim.z ) + threadIdx.z;
    const IndexType i  = ij / sectionSizesD[1] ;
    const IndexType j  = ij - i * sectionSizesD[1];

    if ( m >= sectionSizesD[3] || k >= sectionSizesD[2] || j >= sectionSizesD[1]  || i >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i * targetDistancesD[0] + j * targetDistancesD[1]
                                    + k * targetDistancesD[2] + m * targetDistancesD[3]];

    if ( swapOperands )
    {   
        target = applyBinary( val, op, target );
    }
    else
    {   
        target = applyBinary( target, op, val );
    }
}

template<typename ValueType>
void CUDASection::assignScalar4(
    ValueType targetSection[],
    const ValueType val,
    const IndexType sizes[],
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    SCAI_REGION( "CUDA.Section.assignScalar4" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];
    IndexType n2 = sizes[2];
    IndexType n3 = sizes[3];

    SCAI_LOG_INFO( logger,  "assignScalar2<" << common::TypeTraits<ValueType>::id() 
                            << "> on " << n0 << " x " << n1 << " x " << n2 << " x " << n3 << " grid" )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n3 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n2 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n0 * n1 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    assignScalar3Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, val, op, swapOperands );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "assignScalar1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASection::assignScalar( 
    ValueType section[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType distances[],
    ValueType val,
    const common::binary::BinaryOp op,
    const bool swapOperands )
{
    if ( nDims == 0 )
    {
        // apply kernel on scalar arguments, no section sizes, differences are needed

        assignScalar0( section, val, op, swapOperands );
        return;
    }

    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( sectionSizesD, sizes, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( targetDistancesD, distances, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );

    switch ( nDims )
    {
        case 1 : assignScalar1( section, val, sizes, op, swapOperands );
                 break;

        case 2 : assignScalar2( section, val, sizes, op, swapOperands );
                 break;

        case 3 : assignScalar3( section, val, sizes, op, swapOperands );
                 break;

        case 4 : assignScalar4( section, val, sizes, op, swapOperands );
                 break;

        default: COMMON_THROWEXCEPTION( "assign for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void unary1Kernel( 
    ValueType section[],
    const common::unary::UnaryOp op )
{
    const IndexType i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = section[i * targetDistancesD[0]];
    target = applyUnary( op, target );
}

template<typename ValueType>
void CUDASection::unary1(
    ValueType section[],
    const IndexType sizes[],
    const common::unary::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.unary1" )

    IndexType n0 = sizes[0];

    SCAI_LOG_INFO( logger,  "unary1<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " grid" )

    dim3 threadsPerBlock( 256 );
    dim3 numBlocks( ( n0 + threadsPerBlock.x - 1 ) / threadsPerBlock.x );

    unary1Kernel<<< numBlocks, threadsPerBlock>>>( section, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "unary1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void unary2Kernel( 
    ValueType section[],
    const common::unary::UnaryOp op )
{
    const IndexType j = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i = blockIdx.y * blockDim.y + threadIdx.y;

    if ( j >= sectionSizesD[1]  || i >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = section[i * targetDistancesD[0] + j * targetDistancesD[1]];
    target = applyUnary( op, target );
}

template<typename ValueType>
void CUDASection::unary2(
    ValueType section[],
    const IndexType sizes[],
    const common::unary::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.unary2" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[0];

    SCAI_LOG_INFO( logger,  "unary2<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " x " << n1 << " grid" )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    unary2Kernel<<< numBlocks, threadsPerBlock>>>( section, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "unary2Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASection::unary( 
    ValueType section[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType distances[],
    const common::unary::UnaryOp op )
{
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( sectionSizesD, sizes, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( targetDistancesD, distances, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    switch ( nDims )
    {
        case 1 : unary1( section, sizes, op );
                 break;

        case 2 : unary2( section, sizes, op );
                 break;

        default: COMMON_THROWEXCEPTION( "unary for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */
/*    struct ArrayKernels methods                                              */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASection::ArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_DEBUG( logger, "register SectionKernel CUDA routines for GPU at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<SectionKernelTrait::assign<ValueType> >( assign, ctx, flag );
    KernelRegistry::set<SectionKernelTrait::assignScalar<ValueType> >( assignScalar, ctx, flag );
    KernelRegistry::set<SectionKernelTrait::unary<ValueType> >( unary, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void CUDASection::BinOpKernels<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;

    SCAI_LOG_DEBUG( logger, "register SectionKernel CUDA-routines for GPUs at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<SectionKernelTrait::unaryOp<ValueType, OtherValueType> >( unaryOp, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDASection::CUDASection()
{
    SCAI_LOG_INFO( logger, "register section CUDA kernels for GPUs" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    kregistry::mepr::RegistratorV<ArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<BinOpKernels, SCAI_ARRAY_TYPES_CUDA_LIST, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDASection::~CUDASection()
{
    SCAI_LOG_INFO( logger, "unregister section CUDA kernels for GPUs" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<ArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<BinOpKernels, SCAI_ARRAY_TYPES_CUDA_LIST, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

CUDASection CUDASection::guard;

} /* end namespace utilskernel */

} /* end namespace scai */
