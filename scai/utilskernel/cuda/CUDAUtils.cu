/**
 * @file utilskernel/cuda/CUDAUtils.cu
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
 * @brief Implementation of CSR utilities with CUDA
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/utilskernel/cuda/CUDAUtils.hpp>
#include <scai/utilskernel/cuda/CUDASparseUtils.hpp>



// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/Math.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/SparseKernelTrait.hpp>

// thrust
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>


#include <complex.h>


using namespace scai::common;

namespace scai
{

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( CUDAUtils::logger, "CUDA.Utils" )

/* --------------------------------------------------------------------------- */
/*                                 UnaryOp kernels                               */
/* --------------------------------------------------------------------------- */

/** This kernel can be applied for any arbitrary UnaryOp operation */

template<typename ValueType>
__global__
void unaryOpKernel( ValueType out[], const UnaryOp op, const ValueType in[], IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = applyUnary( op, in[i] );
    }
}

/** The following kernel is a specialization of unaryOp with op == UnaryOp::CONJ
 *
 *  Note: latest CUDA compiler releases show no performance benefits but we keep it
 *        here to demonstrate how specific kernels might be used for optimization
 */

template<typename ValueType>
__global__
void conjKernel( ValueType out[], const ValueType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = Math::conj( in[i] );
    }
}

/** A specialization of conj kernel is required at Math::conj is not defined for IndexType  */

template<>
__global__
void conjKernel( IndexType out[], const IndexType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = in[i];
    }
}

/* abs */

template<typename ValueType>
__global__
void absKernel( ValueType out[], const ValueType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = Math::abs( in[i] );
    }
}

/* floor */

template<typename ValueType>
__global__
void floorKernel( ValueType out[], const ValueType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = Math::floor( in[i] );
    }
}

template<>
__global__
void floorKernel( IndexType out[], const IndexType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = in[i];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*   binaryOp kernels for two input arrays                                                                            */
/* ------------------------------------------------------------------------------------------------------------------ */

/** The following kernel can be used as binary operation on arrays for each operation.
 *  Due to inlining the code has the same efficiency than using any kernel were the operation
 *  is explicitly coded.
 */

template<typename ValueType>
__global__
void binOpKernel( ValueType out[], const ValueType in1[], const common::BinaryOp op, const ValueType in2[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = applyBinary( in1[i], op, in2[i] );
    }
}

/** Special version of binOpKernel with op == BinaryOp::MULT */

template<typename ValueType>
__global__
void multKernel( ValueType out[], const ValueType in1[], const ValueType in2[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = in1[i] * in2[i];
    }
}

/** Special version of binOpKernel with op == common::BinaryOp::ADD */

template<typename ValueType>
__global__
void addKernel( ValueType out[], const ValueType in1[], const ValueType in2[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = in1[i] + in2[i];
    }
}

template<typename ValueType>
__global__
void binOpScalar1Kernel( ValueType out[], const ValueType value, const common::BinaryOp op, const ValueType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = applyBinary( value, op, in[i] );
    }
}

template<typename ValueType>
__global__
void binOpScalar2Kernel( ValueType out[], const ValueType in[], const BinaryOp op, const ValueType value, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = applyBinary( in[i], op, value );
    }
}

template<typename ValueType>
__global__
void addScalarKernel( ValueType out[], const ValueType in[], const ValueType value, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = in[i] + value;
    }
}

template<typename ValueType>
__global__
void scaleVectorAddScalarKernel( ValueType array1[], const ValueType array2[], const ValueType alpha, const ValueType beta, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array1[i] = alpha * array2[i] + beta;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*   BinaryOp::MULT kernel                                                                                              */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void multScalarKernel( ValueType out[], const ValueType in[], const ValueType value, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = in[i] * value;
    }
}

template<typename ValueType>
__global__
void divScalarKernel( ValueType out[], const ValueType value, const ValueType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = value / in[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const BinaryOp op )
{
    SCAI_REGION( "CUDA.Utils.setVal" )

    using namespace thrust::placeholders;

    SCAI_LOG_INFO( logger, "setVal # array = " << array << ", n = " << n << ", val = " << val << ", op = " << op )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

    ValueType value = static_cast<ValueType>( val );

    switch ( op )
    {
        case BinaryOp::COPY:
            thrust::fill( data, data + n, value );
            break;

        case BinaryOp::ADD:
            thrust::for_each( data, data + n,  _1 += value );
            break;

        case BinaryOp::SUB:
            thrust::for_each( data, data + n,  _1 -= value );
            break;

        case BinaryOp::MULT:
        {
            if ( val == scai::common::Constants::ZERO )
            {
                thrust::fill( data, data + n, ValueType( 0 ) );
            }
            else
            {
                thrust::for_each( data, data + n,  _1 *= value );
            }
        }
        break;

        case BinaryOp::DIVIDE:
        {
            if ( val == scai::common::Constants::ZERO )
            {
                COMMON_THROWEXCEPTION( "Divide by ZERO" )
            }
            else
            {
                thrust::for_each( data, data + n,  _1 /= value );
            }
        }
        break;

        default:
            COMMON_THROWEXCEPTION( "unsupported binary op: " << op )
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

template<typename ValueType>
void CUDAUtils::scaleVectorAddScalar( ValueType array1[], const ValueType array2[], const IndexType n,
                                      const ValueType alpha, const ValueType beta )
{
    SCAI_REGION( "CUDA.Utils.scaleVectorAddScalar" )

    SCAI_LOG_INFO( logger, "scaleVectorAddScalar<" << TypeTraits<ValueType>::id() << ">( ..., n = " << n << ")" )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    scaleVectorAddScalarKernel<ValueType> <<< dimGrid, dimBlock>>>( array1, array2, alpha, beta, n );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::getValue( const ValueType* array, const IndexType i )
{
    SCAI_REGION( "CUDA.Utils.getValue" )
    SCAI_LOG_INFO( logger, "getValue # i = " << i )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> arrayPtr( const_cast<ValueType*>( array ) );
    thrust::host_vector<ValueType> arrayHost( arrayPtr + i, arrayPtr + i + 1 );
    return arrayHost[0];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void scatterVal_kernel( ValueType* out, const IndexType* indexes, const ValueType value, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[indexes[i]] = value;
    }
}

template<typename ValueType>
void CUDAUtils::scatterVal( ValueType out[], const IndexType indexes[], const ValueType val, const IndexType n  )
{
    SCAI_LOG_INFO( logger,
                   "scatterVal<" << TypeTraits<ValueType>::id() << ">( ..., n = " << n << ")" )

    if ( n > 0 )
    {
        SCAI_CHECK_CUDA_ACCESS
        const int blockSize = 256;
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( n, dimBlock.x );
        scatterVal_kernel <<< dimGrid, dimBlock>>>( out, indexes, val, n );
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::unaryOp( ValueType out[], const ValueType in[], const IndexType n, const UnaryOp op )
{
    SCAI_REGION( "CUDA.Utils.unaryOp" )

    SCAI_LOG_INFO( logger, "unaryOp<" << TypeTraits<ValueType>::id() << ">( ..., n = " << n << ")" )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    switch ( op )
    {
        case UnaryOp::CONJ :
        {
            conjKernel <<< dimGrid, dimBlock>>>( out, in, n );

            break;
        }

        case UnaryOp::MINUS :
        {
            multScalarKernel <<< dimGrid, dimBlock>>>( out, in, ValueType( -1 ), n );

            break;
        }

        case UnaryOp::ABS :
        {
            absKernel <<< dimGrid, dimBlock>>>( out, in, n );

            break;
        }

        case UnaryOp::FLOOR :
        {
            // the performance of this kernel should be compared to
            // unaryOp with op == UnaryOp::CEIL to decide whether it is worth

            floorKernel <<< dimGrid, dimBlock>>>( out, in, n );

            break;
        }

        default:

            unaryOpKernel <<< dimGrid, dimBlock>>>( out, op, in, n );
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::binaryOp(
    ValueType out[],
    const ValueType in1[],
    const ValueType in2[],
    const IndexType n,
    const BinaryOp op )
{
    SCAI_REGION( "CUDA.Utils.binOp" )

    SCAI_LOG_INFO( logger, "binaryOp<" << TypeTraits<ValueType>::id() << ">( ..., n = " << n << ")" )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    // Specific kernels for a given binary op does not give relevant performance gain

    switch ( op )
    {
        case BinaryOp::ADD :
        {
            addKernel<ValueType> <<< dimGrid, dimBlock>>>( out, in1, in2, n );
            break;
        }

        case BinaryOp::MULT :
        {
            multKernel<ValueType> <<< dimGrid, dimBlock>>>( out, in1, in2, n );
            break;
        }

        default:
        {
            binOpKernel<ValueType> <<< dimGrid, dimBlock>>>( out, in1, op, in2, n );
        }
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::binaryOpScalar(
    ValueType out[],
    const ValueType in[],
    const ValueType value,
    const IndexType n,
    const BinaryOp op,
    const bool swapScalar )
{
    SCAI_REGION( "CUDA.Utils.binOpScalar" )

    SCAI_LOG_INFO( logger, "binaryOp<" << TypeTraits<ValueType>::id() << ">( ..., n = " << n << ")" )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    switch ( op )
    {
        case BinaryOp::ADD :
        {
            addScalarKernel<ValueType> <<< dimGrid, dimBlock>>>( out, in, value, n );
            break;
        }

        case BinaryOp::MULT :
        {
            if ( value == scai::common::Constants::ZERO )
            {
                CUDAUtils::setVal( out, n, value, BinaryOp::COPY );
                return;
            }

            multScalarKernel<ValueType> <<< dimGrid, dimBlock>>>( out, in, value, n );

            break;
        }

        case BinaryOp::DIVIDE :
        {
            if ( swapScalar )
            {
                if ( value == scai::common::Constants::ZERO )
                {
                    CUDAUtils::setVal( out, n, value, BinaryOp::COPY );
                    return;
                }
    
                divScalarKernel<ValueType> <<< dimGrid, dimBlock>>>( out, value, in, n );
            }
            else
            {
                ValueType factor = ValueType( 1 ) / value;
                multScalarKernel<ValueType> <<< dimGrid, dimBlock>>>( out, in, factor, n );
            }

            break;
        }
        default:
        {
            if ( swapScalar )
            {
                binOpScalar1Kernel<ValueType> <<< dimGrid, dimBlock>>>( out, value, op, in, n );
            } 
            else
            {
                binOpScalar2Kernel<ValueType> <<< dimGrid, dimBlock>>>( out, in, op, value, n );
            }
        }
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "kernel for binary op with scalar" );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDAUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    // const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag << "]" )
    // we keep the registrations for IndexType as we do not need conversions
}

template<typename ValueType>
void CUDAUtils::RegArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;

    SCAI_LOG_DEBUG( logger, "registerV array UtilsKernel CUDA [" << flag
                    << "] --> ValueType = " << common::getScalarType<ValueType>() )

    // Note: these kernels will be instantiated for numeric types + IndexType

    KernelRegistry::set<UtilKernelTrait::getValue<ValueType> >( getValue, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setVal<ValueType> >( setVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scaleVectorAddScalar<ValueType> >( scaleVectorAddScalar, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scatterVal<ValueType> >( scatterVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::unaryOp<ValueType> >( unaryOp, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::binaryOp<ValueType> >( binaryOp, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::binaryOpScalar<ValueType> >( binaryOpScalar, ctx, flag );
}

template<typename ValueType, typename SourceValueType>
void CUDAUtils::RegistratorVO<ValueType, SourceValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    // const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "registerVO UtilsKernel CUDA [" << flag
                    << "] --> ValueType = " << common::getScalarType<ValueType>()
                    << ", SourceValueType = " << common::getScalarType<SourceValueType>() )
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDAUtils::CUDAUtils()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARRAY_TYPES_CUDA_LIST, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAUtils::~CUDAUtils()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARRAY_TYPES_CUDA_LIST, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAUtils CUDAUtils::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
