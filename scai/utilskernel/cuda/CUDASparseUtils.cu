/**
 * @file utilskernel/cuda/CUDASparseUtils.cu
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
 * @brief Implementation of utility kernels with CUDA, here sparse and set routines
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/utilskernel/cuda/CUDASparseUtils.hpp>


// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAUtils.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/Math.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/SparseKernelTrait.hpp>

// thrust
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/transform_reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>

#include <complex.h>


using namespace scai::common;

namespace scai
{

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( CUDASparseUtils::logger, "CUDA.SparseUtils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASparseUtils::setOrder( ValueType array[], const IndexType n )
{
    SCAI_REGION( "CUDA.Utils.setOrder" )
    SCAI_LOG_INFO( logger, "setOrder # array = " << array << ", n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> array_ptr( const_cast<ValueType*>( array ) );
    thrust::sequence( array_ptr, array_ptr + n );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDASparseUtils::setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n )
{
    SCAI_REGION( "CUDA.Utils.setSequence" )
    SCAI_LOG_INFO( logger, "setSequence # array = " << array << ", n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> array_ptr( const_cast<ValueType*>( array ) );
    thrust::sequence( thrust::device, array_ptr, array_ptr + n, startValue, inc );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  setInversePerm                                                    */
/* ------------------------------------------------------------------------------------------------------------------ */

void CUDASparseUtils::setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n )
{
    SCAI_REGION( "CUDA.Utils.invPerm" )

    SCAI_LOG_INFO( logger, "compute inverse perm, n = " << n )
    SCAI_CHECK_CUDA_ACCESS

    if ( n > 0 )
    {
        thrust::device_ptr<IndexType> inversePermPtr( const_cast<IndexType*>( inversePerm ) );
        thrust::device_ptr<IndexType> permPtr( const_cast<IndexType*>( perm ) );
        thrust::counting_iterator<IndexType> sequence( 0 );
        thrust::scatter( sequence, sequence + n, permPtr, inversePermPtr );
        SCAI_CHECK_CUDA_ERROR
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, typename SourceValueType>
__global__
void setKernelCopy( ValueType* out, const SourceValueType* in, const IndexType n )
{   
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    
    if ( i < n )
    {   
        out[i] = static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void setKernelAdd( ValueType* out, const SourceValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] += static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void setKernelSub( ValueType* out, const SourceValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] -= static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void setKernelMult( ValueType* out, const SourceValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] *= static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void setKernelDivide( ValueType* out, const SourceValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] /= static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void setKernelPow( ValueType* out, const SourceValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = Math::pow( out[i], in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
void CUDASparseUtils::set( ValueType out[], const SourceValueType in[], const IndexType n, const BinaryOp op )
{   
    SCAI_REGION( "CUDA.Utils.set" )
    
    SCAI_LOG_INFO( logger,
                   "set<" << TypeTraits<ValueType>::id() << "," << TypeTraits<SourceValueType>::id() << ">( ..., n = " << n << ")" )
    SCAI_LOG_DEBUG( logger, "out = " << out << ", in = " << in )
    
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
        case BinaryOp::COPY :
            setKernelCopy <<< dimGrid, dimBlock>>>( out, in, n );
            break;
        
        case BinaryOp::ADD :
            setKernelAdd <<< dimGrid, dimBlock>>>( out, in, n );
            break;
        
        case BinaryOp::SUB :
            setKernelSub <<< dimGrid, dimBlock>>>( out, in, n );
            break;
        
        case BinaryOp::MULT :
            setKernelMult <<< dimGrid, dimBlock>>>( out, in, n );
            break;
        
        case BinaryOp::DIVIDE :
            setKernelDivide <<< dimGrid, dimBlock>>>( out, in, n );
            break;
        
        default:
            COMMON_THROWEXCEPTION( "Unsupported binary op " << op )
    }
    
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "setKernel<Op>" );
}


/* --------------------------------------------------------------------------- */

template<typename ValueType, typename SourceValueType>
__global__
void setKernelCopySection( ValueType* out, const IndexType inc1,
                           const SourceValueType* in, const IndexType inc2, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i * inc1] = static_cast<ValueType>( in[i * inc2] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void setKernelAddSection( ValueType* out, const IndexType inc1,
                          const SourceValueType* in, const IndexType inc2, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i * inc1] += static_cast<ValueType>( in[i * inc2] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void setKernelSubSection( ValueType* out, const IndexType inc1,
                          const SourceValueType* in, const IndexType inc2, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i * inc1] -= static_cast<ValueType>( in[i * inc2] );
    }
}
template<typename ValueType, typename SourceValueType>
__global__
void setKernelMultSection( ValueType* out, const IndexType inc1,
                           const SourceValueType* in, const IndexType inc2, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i * inc1] *= static_cast<ValueType>( in[i * inc2] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void setKernelDivideSection( ValueType* out, const IndexType inc1,
                             const SourceValueType* in, const IndexType inc2, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i * inc1] /= static_cast<ValueType>( in[i * inc2] );
    }
}

template<typename ValueType, typename SourceValueType>
void CUDASparseUtils::setSection( ValueType out[], const IndexType inc1,
                                  const SourceValueType in[], const IndexType inc2,
                                  const IndexType n, const BinaryOp op )
{
    SCAI_REGION( "CUDA.Utils.setSection" )

    SCAI_LOG_INFO( logger,
                   "setSection<" << TypeTraits<ValueType>::id() << "," << TypeTraits<SourceValueType>::id()
                   << "> : out(:" << inc2 << ") <- in(:" << inc1 << "), n = " << n << ", op = " << op )

    SCAI_LOG_DEBUG( logger, "out = " << out << ", in = " << in )

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
        case BinaryOp::COPY :
            setKernelCopySection <<< dimGrid, dimBlock>>>( out, inc1, in, inc2, n );
            break;
        
        case BinaryOp::ADD :
            setKernelAddSection <<< dimGrid, dimBlock>>>( out, inc1, in, inc2, n );
            break;
        
        case BinaryOp::SUB :
            setKernelSubSection <<< dimGrid, dimBlock>>>( out, inc1, in, inc2, n );
            break;
        
        case BinaryOp::MULT :
            setKernelMultSection <<< dimGrid, dimBlock>>>( out, inc1, in, inc2, n );
            break;
        
        case BinaryOp::DIVIDE :
            setKernelDivideSection <<< dimGrid, dimBlock>>>( out, inc1, in, inc2, n );
            break;
        
        default:
            COMMON_THROWEXCEPTION( "Unsupported binary op " << op )
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */
/*    gather + gather kernels                                                  */
/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
__global__
void gatherKernel(
    ValueType1 out[],
    const ValueType2 in[],
    const IndexType indexes[],
    const BinaryOp op,
    const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = applyBinary( out[i], op, static_cast<ValueType1>( in[indexes[i]] ) );
    }
}

template<typename ValueType1, typename ValueType2>
__global__
void gatherCopyKernel( ValueType1 out[], const ValueType2 in[], const IndexType indexes[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = static_cast<ValueType1>( in[indexes[i]] );
    }
}

template<typename ValueType1, typename ValueType2>
void CUDASparseUtils::setGather(
    ValueType1 out[],
    const ValueType2 in[],
    const IndexType indexes[],
    const BinaryOp op,
    const IndexType n )
{   
    SCAI_REGION( "CUDA.Utils.setGather" )
    
    SCAI_LOG_INFO( logger, "setGather<" << TypeTraits<ValueType1>::id() << "," << TypeTraits<ValueType2>::id()
                   << ">, n = " << n << ", op = " << op )
    
    SCAI_CHECK_CUDA_ACCESS
    
    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );
    
    switch ( op )
    {   
        case BinaryOp::COPY :
            gatherCopyKernel <<< dimGrid, dimBlock>>>( out, in, indexes, n );
            break;
        
        default:
            gatherKernel <<< dimGrid, dimBlock>>>( out, in, indexes, op, n );
    }
    
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "setGather::kernel" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, typename SourceValueType>
__global__
void scatter_kernel( ValueType* out, const IndexType* indexes, const SourceValueType* in, const IndexType n )
{   
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    
    if ( i < n )
    {   
        out[indexes[i]] = static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void scatter_add_kernel( ValueType* out, const IndexType* indexes, const SourceValueType* in, const IndexType n )
{   
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    
    if ( i < n )
    {   
        out[indexes[i]] += static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void scatter_add1_kernel( ValueType* out, const IndexType* indexes, const SourceValueType* in, const IndexType n )
{   
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    
    if ( i < n )
    {   
        common::CUDAUtils::atomicAdd( &out[indexes[i]], static_cast<ValueType>( in[i] ) );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void scatter_sub_kernel( ValueType out[], const IndexType indexes[], const SourceValueType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[indexes[i]] -= static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void scatter_sub1_kernel( ValueType out[], const IndexType indexes[], const SourceValueType in[], const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        common::CUDAUtils::atomicAdd( &out[indexes[i]], -static_cast<ValueType>( in[i] ) );
    }
}

template<typename ValueType, typename SourceValueType>
__global__
void scatter_op_kernel( ValueType* out, const IndexType* indexes, const SourceValueType* in, const IndexType n, const BinaryOp op )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[indexes[i]] = applyBinary( out[indexes[i]], op, static_cast<ValueType>( in[i] ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void CUDASparseUtils::setScatter(
    ValueType1 out[],
    const IndexType indexes[],
    const bool unique,
    const ValueType2 in[],
    const BinaryOp op,
    const IndexType n )
{
    SCAI_REGION( "CUDA.Utils.setScatter" )

    SCAI_LOG_INFO( logger,
                   "setScatter<" << TypeTraits<ValueType1>::id() << "," << TypeTraits<ValueType2>::id() 
                   << ">( ..., n = " << n << "), op = " << op << ", unique = " << unique )

    if ( n > 0 )
    {
        SCAI_CHECK_CUDA_ACCESS
        const int blockSize = 256;
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( n, dimBlock.x );

        if ( op == BinaryOp::COPY )
        {
            // unique does not matter, result can depend on race conditions

            scatter_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n );
        }
        else if ( op == BinaryOp::ADD )
        {
            if ( unique )
            {
                scatter_add_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n );
            }
            else
            {
                scatter_add1_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n );
            }
        }
        else if ( op == BinaryOp::SUB )
        {
            if ( unique )
            {
                scatter_sub_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n );
            }
            else
            {
                scatter_sub1_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n );
            }
        }
        else if ( unique )
        {
            scatter_op_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n, op );
        }
        else
        {
            COMMON_THROWEXCEPTION( "setScatter with atomic update not supported for op = " << op )
        }

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    }
}

/* ---------------------------------------------------------------------------- */

template<typename ValueType>
struct nonZero
{
    // RealType is required otherwise operator > might be undefined 

    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    const ValueType mZero;
    const RealType mEps;

    nonZero( ValueType zero, ValueType eps ) :
        mZero( zero ),
        mEps( Math::abs( eps ) )
    {
    }

    __host__ __device__
    bool operator()( ValueType x )
    {
        RealType absX = Math::abs( x - mZero );
        return absX > mEps;
    }
};

template<typename ValueType>
IndexType CUDASparseUtils::countNonZeros( const ValueType denseArray[], const IndexType n, const ValueType zero, const ValueType eps )
{
    SCAI_REGION( "CUDA.Utils.countNZ" )

    SCAI_LOG_INFO( logger, "countNonZeros of array[ " << n << " ]" )

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> array_ptr( const_cast<ValueType*>( denseArray ) );

    IndexType nonZeros = thrust::transform_reduce( array_ptr,
                         array_ptr + n,
                         nonZero<ValueType>( zero, eps ),
                         0,
                         thrust::plus<IndexType>() );
    return nonZeros;
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*     compress                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
struct invalidateZeros
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    const ValueType mZero;
    const RealType mEps;
    const IndexType mZeroIndex;

    invalidateZeros( ValueType zero, ValueType eps, IndexType zeroIndex ) :

        mZero( zero ),
        mEps( eps ),
        mZeroIndex( zeroIndex )
    {
    }

    __host__ __device__
    IndexType operator()( const ValueType& value, const IndexType& index )
    {
        RealType tmp = common::Math::abs( value - mZero );

        if ( tmp > mEps )       // true -> value is considered as non-zero
        {
            return index;
        }
        else
        {
            return mZeroIndex;
        }
    }
};

struct notEqual
{
    const IndexType mOutOfRange;

    notEqual( const IndexType val ) : mOutOfRange( val )
    {
    }

    __host__ __device__
    bool operator()( const IndexType x )
    {
        return x != mOutOfRange;
    }
};

template<typename ValueType, typename SourceType>
IndexType CUDASparseUtils::compress(
    ValueType sparseArray[],
    IndexType sparseIndexes[],
    const SourceType denseArray[],
    const IndexType n,
    const SourceType zero,
    const SourceType eps )
{
    SCAI_REGION( "CUDA.Utils.compress" )

    SCAI_LOG_INFO( logger, "compress<" << TypeTraits<ValueType>::id() << ", "
                           << TypeTraits<SourceType>::id() << "> array[ " << n << " ]" )

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<SourceType> dense_ptr( const_cast<SourceType*>( denseArray ) );
    thrust::counting_iterator<IndexType> sequence( IndexType( 0 ) );
    thrust::device_vector<IndexType> tmp( n );

    // Create device ptr and help variables

    thrust::transform( dense_ptr, dense_ptr + n, sequence, tmp.begin(), invalidateZeros<ValueType>( zero, eps, invalidIndex ) );

    // transform array, replace in sequence all non-zero entries with -1
    // e.g. sizes = [ 0, 2, 0, 1, 3 ], sequence = [ 0, 1, 2, 3, 4 ] -> [ invalidIndex, 1, invalidIndex, 3, 4 ]

    thrust::device_ptr<IndexType> sparseIndexes_ptr( sparseIndexes );

    // now compact all indexes with values not equal invalidIndex

    IndexType cnt = thrust::copy_if( tmp.begin(),
                                     tmp.end(),
                                     sparseIndexes_ptr,
                                     notEqual( invalidIndex ) ) - sparseIndexes_ptr;

    if ( sparseArray )
    {
        const int blockSize = 256;
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( cnt, dimBlock.x );

        gatherCopyKernel <<< dimGrid, dimBlock>>>( sparseArray, denseArray, sparseIndexes, cnt );
    }

    return cnt;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDASparseUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<UtilKernelTrait::setInversePerm>( setInversePerm, ctx, flag );
}

template<typename ValueType>
void CUDASparseUtils::RegArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;

    SCAI_LOG_DEBUG( logger, "registerV array UtilsKernel CUDA [" << flag
                    << "] --> ValueType = " << common::getScalarType<ValueType>() )

    // Note: these kernels will be instantiated for numeric types + IndexType

    KernelRegistry::set<UtilKernelTrait::setOrder<ValueType> >( setOrder, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setSequence<ValueType> >( setSequence, ctx, flag );
    KernelRegistry::set<SparseKernelTrait::countNonZeros<ValueType> >( countNonZeros, ctx, flag );
}

template<typename ValueType, typename SourceValueType>
void CUDASparseUtils::RegistratorVO<ValueType, SourceValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    KernelRegistry::set<SparseKernelTrait::compress<ValueType, SourceValueType> >( compress, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setScatter<ValueType, SourceValueType> >( setScatter, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<ValueType, SourceValueType> >( setGather, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::set<ValueType, SourceValueType> >( set, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setSection<ValueType, SourceValueType> >( setSection, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDASparseUtils::CUDASparseUtils()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARRAY_TYPES_CUDA_LIST, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDASparseUtils::~CUDASparseUtils()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARRAY_TYPES_CUDA_LIST, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDASparseUtils CUDASparseUtils::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
