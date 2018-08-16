/**
 * @file CUDASection.cu
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

using common::BinaryOp;
using common::UnaryOp;
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i0 = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i0 > 0 )
    {
        return;  // only one element
    }

    if ( swapOperands )
    {   
        target[i0] = applyBinary( source[i0], op, target[i0] );
    }
    else
    {   
        target[i0] = applyBinary( target[i0], op, source[i0] );
    }
}

template<typename ValueType>
void CUDASection::assign0(
    ValueType targetSection[],
    const ValueType sourceSection[],
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i0 = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i0 >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = targetSection[i0 * targetDistancesD[0]];
    const ValueType& source = sourceSection[i0 * sourceDistancesD[0]];

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
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i1 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i0 = blockIdx.y * blockDim.y + threadIdx.y;

    if ( i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1]];
    const ValueType& source = sourceSection[i0 * sourceDistancesD[0] + i1 * sourceDistancesD[1]];

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
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i2 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i1 = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i0 = blockIdx.z;

    if ( i1 >= sectionSizesD[1]  || i2 >= sectionSizesD[2] )
    {
        return;
    }

    ValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1] + i2 * targetDistancesD[2]];
    const ValueType& source = sourceSection[i0 * sourceDistancesD[0] + i1 * sourceDistancesD[1] + i2 * sourceDistancesD[2]];

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
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i3 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i2 = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i10 = ( blockIdx.z * blockDim.z ) + threadIdx.z;
    const IndexType i0  = i10 / sectionSizesD[1] ;
    const IndexType i1  = i10 - i0 * sectionSizesD[1];

    // i and i1 will always be legal, but i3 or might be out of range

    if ( i3 >= sectionSizesD[3]  || i2 >= sectionSizesD[2] || i1 >= sectionSizesD[1] || i0 >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1] 
                                    + i2 * targetDistancesD[2] + i3 * targetDistancesD[3]];

    const ValueType& source = sourceSection[i0 * sourceDistancesD[0] + i1 * sourceDistancesD[1] 
                                          + i2 * sourceDistancesD[2] + i3 * sourceDistancesD[3]];

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
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    SCAI_CHECK_CUDA_ACCESS

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
    const common::UnaryOp op )
{
    const IndexType i0 = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i0 >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    TargetValueType& target = targetSection[i0 * targetDistancesD[0]];
    const SourceValueType& source = sourceSection[i0 * sourceDistancesD[0]];

    target = static_cast<TargetValueType>( applyUnary( op, source ) );
}

template<typename TargetValueType, typename SourceValueType>
void CUDASection::unaryOp1(
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const IndexType sizes[],
    const common::UnaryOp op )
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
    const common::UnaryOp op )
{
    const IndexType i1 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i0 = blockIdx.y * blockDim.y + threadIdx.y;
      
    if ( i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    TargetValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1]];
    const SourceValueType& source = sourceSection[i0 * sourceDistancesD[0] + i1 * sourceDistancesD[1]];

    target = static_cast<TargetValueType>( applyUnary( op, source ) );
}

template<typename TargetValueType, typename SourceValueType>
void CUDASection::unaryOp2(
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const IndexType sizes[],
    const common::UnaryOp op )
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
__global__
void unaryOp3Kernel( 
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const common::UnaryOp op )
{
    const IndexType i2 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i1 = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i0 = blockIdx.z * blockDim.z + threadIdx.z;

    if ( i2 >= sectionSizesD[2] || i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;
    }

    TargetValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1] + i2 * targetDistancesD[2]];
    const SourceValueType& source = sourceSection[i0 * sourceDistancesD[0] + i1 * sourceDistancesD[1] + i2 * sourceDistancesD[2]];

    target = static_cast<TargetValueType>( applyUnary( op, source ) );
}

template<typename TargetValueType, typename SourceValueType>
void CUDASection::unaryOp3(
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const IndexType sizes[],
    const common::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.unaryOp3" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];
    IndexType n2 = sizes[2];

    SCAI_LOG_INFO( logger, "unaryOp3<" << common::TypeTraits<TargetValueType>::id() << 
                           ", " << common::TypeTraits<SourceValueType>::id() 
                           << "> on " << n0 << " x " << n1 << " x " << n2 << " grid" )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n2 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n1 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n0 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    unaryOp3Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "unaryOp2Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
__global__
void unaryOp4Kernel( 
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const common::UnaryOp op )
{
    const IndexType i3  = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i2  = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i10 = ( blockIdx.z * blockDim.z ) + threadIdx.z;
    const IndexType i0  = i10 / sectionSizesD[1] ;
    const IndexType i1  = i10 - i0 * sectionSizesD[1];

    if ( i3 >= sectionSizesD[3] || i2 >= sectionSizesD[2] || i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;
    }

    TargetValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1]
                                          + i2 * targetDistancesD[2] + i3 * targetDistancesD[3]];

    const SourceValueType& source = sourceSection[i0 * sourceDistancesD[0] + i1 * sourceDistancesD[1] 
                                                + i2 * sourceDistancesD[2] + i3 * sourceDistancesD[3]];

    target = static_cast<TargetValueType>( applyUnary( op, source ) );
}

template<typename TargetValueType, typename SourceValueType>
void CUDASection::unaryOp4(
    TargetValueType targetSection[],
    const SourceValueType sourceSection[],
    const IndexType sizes[],
    const common::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.unaryOp3" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];
    IndexType n2 = sizes[2];
    IndexType n3 = sizes[3];

    SCAI_LOG_INFO( logger, "unaryOp3<" << common::TypeTraits<TargetValueType>::id() << 
                           ", " << common::TypeTraits<SourceValueType>::id() 
                           << "> on " << n0 << " x " << n1 << " x " << n2 << " x " << n3 << " grid" )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n3 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n2 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n0 * n1 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    unaryOp4Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, sourceSection, op );

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
    const common::UnaryOp op )
{
    SCAI_CHECK_CUDA_ACCESS

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

        case 3 : unaryOp3( targetSection, sourceSection, sizes, op );
                 break;

        case 4 : unaryOp4( targetSection, sourceSection, sizes, op );
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i0 = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i0 > 0 )
    {
        return;  // only one element
    }

    if ( swap )
    {   
        target[i0] = applyBinary( source, op, target[i0] );
    }
    else
    {   
        target[i0] = applyBinary( target[i0], op, source );
    }
}

template<typename ValueType>
void CUDASection::assignScalar0(
    ValueType targetSection[],
    const ValueType source,
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i0 = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i0 >= sectionSizesD[0] )
    {   
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = targetSection[i0 * targetDistancesD[0]];

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
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i1 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i0 = blockIdx.y * blockDim.y + threadIdx.y;

    if ( i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1]];

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
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i2 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i1 = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i0 = blockIdx.z * blockDim.z + threadIdx.z;

    if ( i2 >= sectionSizesD[2] || i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1] + i2 * targetDistancesD[2]];

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
    const BinaryOp op,
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
    const BinaryOp op,
    const bool swapOperands )
{
    const IndexType i3 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i2 = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i10 = ( blockIdx.z * blockDim.z ) + threadIdx.z;
    const IndexType i0  = i10 / sectionSizesD[1] ;
    const IndexType i1  = i10 - i0 * sectionSizesD[1];

    if ( i3 >= sectionSizesD[3] || i2 >= sectionSizesD[2] || i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;
    }

    ValueType& target = targetSection[i0 * targetDistancesD[0] + i1 * targetDistancesD[1]
                                    + i2 * targetDistancesD[2] + i3 * targetDistancesD[3]];

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
    const BinaryOp op,
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

    assignScalar4Kernel<<< numBlocks, threadsPerBlock>>>( targetSection, val, op, swapOperands );

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
    const BinaryOp op,
    const bool swapOperands )
{
    SCAI_CHECK_CUDA_ACCESS

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
void UnaryOp1Kernel( 
    ValueType section[],
    const common::UnaryOp op )
{
    const IndexType i0 = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i0 >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = section[i0 * targetDistancesD[0]];
    target = applyUnary( op, target );
}

template<typename ValueType>
void CUDASection::UnaryOp1(
    ValueType section[],
    const IndexType sizes[],
    const common::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.UnaryOp1" )

    IndexType n0 = sizes[0];

    SCAI_LOG_INFO( logger,  "UnaryOp1<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " grid" )

    dim3 threadsPerBlock( 256 );
    dim3 numBlocks( ( n0 + threadsPerBlock.x - 1 ) / threadsPerBlock.x );

    UnaryOp1Kernel<<< numBlocks, threadsPerBlock>>>( section, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "UnaryOp1Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void UnaryOp2Kernel( 
    ValueType section[],
    const common::UnaryOp op )
{
    const IndexType i1 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i0 = blockIdx.y * blockDim.y + threadIdx.y;

    if ( i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = section[i0 * targetDistancesD[0] + i1 * targetDistancesD[1]];
    target = applyUnary( op, target );
}

template<typename ValueType>
void CUDASection::UnaryOp2(
    ValueType section[],
    const IndexType sizes[],
    const common::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.UnaryOp2" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];

    SCAI_LOG_INFO( logger,  "UnaryOp2<" << common::TypeTraits<ValueType>::id() << "> on " << n0 << " x " << n1 << " grid" )

    dim3 threadsPerBlock( 16, 16 );

    dim3 numBlocks( ( n1 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n0 + threadsPerBlock.y - 1 ) / threadsPerBlock.y );

    UnaryOp2Kernel<<< numBlocks, threadsPerBlock>>>( section, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "UnaryOp2Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void UnaryOp3Kernel( 
    ValueType section[],
    const common::UnaryOp op )
{
    const IndexType i2 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i1 = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i0 = blockIdx.z * blockDim.z + threadIdx.z;

    if ( i2 >= sectionSizesD[2] || i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;   // might happen if sectionSizesD[0] is not multiple of blockDim.x
    }

    ValueType& target = section[i0 * targetDistancesD[0] + i1 * targetDistancesD[1] + i2 * targetDistancesD[2]];
    target = applyUnary( op, target );
}

template<typename ValueType>
void CUDASection::UnaryOp3(
    ValueType section[],
    const IndexType sizes[],
    const common::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.UnaryOp3" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];
    IndexType n2 = sizes[2];

    SCAI_LOG_INFO( logger,  "UnaryOp3<" << common::TypeTraits<ValueType>::id() << "> on " 
                                      << n0 << " x " << n1 << " x " << n2 << " grid" )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n2 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n1 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n0 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    UnaryOp3Kernel<<< numBlocks, threadsPerBlock>>>( section, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "UnaryOp3Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void UnaryOp4Kernel( 
    ValueType section[],
    const common::UnaryOp op )
{
    const IndexType i3 = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType i2 = blockIdx.y * blockDim.y + threadIdx.y;
    const IndexType i10 = ( blockIdx.z * blockDim.z ) + threadIdx.z;
    const IndexType i0  = i10 / sectionSizesD[1] ;
    const IndexType i1  = i10 - i0 * sectionSizesD[1];

    if ( i3 >= sectionSizesD[3] || i2 >= sectionSizesD[2] || i1 >= sectionSizesD[1]  || i0 >= sectionSizesD[0] )
    {
        return;   // might happen if sizes[i0] is not multiple of threads in that dim
    }

    ValueType& target = section[i0 * targetDistancesD[0] + i1 * targetDistancesD[1] 
                              + i2 * targetDistancesD[2] + i3 * targetDistancesD[3]];

    target = applyUnary( op, target );
}

template<typename ValueType>
void CUDASection::UnaryOp4(
    ValueType section[],
    const IndexType sizes[],
    const common::UnaryOp op )
{
    SCAI_REGION( "CUDA.Section.UnaryOp4" )

    IndexType n0 = sizes[0];
    IndexType n1 = sizes[1];
    IndexType n2 = sizes[2];
    IndexType n3 = sizes[3];

    SCAI_LOG_INFO( logger,  "UnaryOp4<" << common::TypeTraits<ValueType>::id() << "> on " 
                                      << n0 << " x " << n1 << " x " << n2 << " x " << n3 << " grid" )

    dim3 threadsPerBlock( 16, 4, 4 );

    dim3 numBlocks( ( n3 + threadsPerBlock.x - 1 ) / threadsPerBlock.x,
                    ( n2 + threadsPerBlock.y - 1 ) / threadsPerBlock.y,
                    ( n1 * n0 + threadsPerBlock.z - 1 ) / threadsPerBlock.z );

    UnaryOp4Kernel<<< numBlocks, threadsPerBlock>>>( section, op );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "UnaryOp4Kernel failed" ) ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASection::UnaryOp( 
    ValueType section[],
    const IndexType nDims,
    const IndexType sizes[],
    const IndexType distances[],
    const common::UnaryOp op )
{
    SCAI_CHECK_CUDA_ACCESS

    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( sectionSizesD, sizes, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    SCAI_CUDA_RT_CALL( cudaMemcpyToSymbol( targetDistancesD, distances, nDims * sizeof( IndexType ), 0, cudaMemcpyHostToDevice ),
                       "copy2Device failed" );
    switch ( nDims )
    {
        case 1 : UnaryOp1( section, sizes, op );
                 break;

        case 2 : UnaryOp2( section, sizes, op );
                 break;

        case 3 : UnaryOp3( section, sizes, op );
                 break;

        case 4 : UnaryOp4( section, sizes, op );
                 break;

        default: COMMON_THROWEXCEPTION( "UnaryOp for nDims = " << nDims << " not supported yet" )
    }
}

/* --------------------------------------------------------------------------- */
/*    struct ArrayKernels methods                                              */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASection::ArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register SectionKernel CUDA routines for GPU at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<SectionKernelTrait::assign<ValueType> >( assign, ctx, flag );
    KernelRegistry::set<SectionKernelTrait::assignScalar<ValueType> >( assignScalar, ctx, flag );
    KernelRegistry::set<SectionKernelTrait::UnaryOp<ValueType> >( UnaryOp, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void CUDASection::BinOpKernels<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;

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
