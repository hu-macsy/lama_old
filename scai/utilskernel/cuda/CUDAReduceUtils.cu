/**
 * @file utilskernel/cuda/CUDAReduceUtils.cu
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
#include <scai/utilskernel/cuda/CUDAReduceUtils.hpp>
#include <scai/utilskernel/cuda/CUDAUtils.hpp>

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
#include <thrust/reduce.h>
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

SCAI_LOG_DEF_LOGGER( CUDAReduceUtils::logger, "CUDA.ReduceUtils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAReduceUtils::reduceSum( const ValueType array[], const IndexType n, const ValueType zero )
{
    SCAI_REGION( "CUDA.Utils.reduceSum" )
    SCAI_LOG_INFO( logger, "sum # array = " << array << ", n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

    ValueType result = thrust::reduce( data, data + n, zero, thrust::plus<ValueType>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "sum of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAReduceUtils::reduceMaxVal( const ValueType array[], const IndexType n, const ValueType zero )
{
    SCAI_REGION( "CUDA.Utils.reduceMax" )
    SCAI_LOG_INFO( logger, "maxval for " << n << " elements " )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
    ValueType result = thrust::reduce( data, data + n, zero, thrust::maximum<ValueType>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "max of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAReduceUtils::reduceMinVal( const ValueType array[], const IndexType n, const ValueType zero )
{
    SCAI_REGION( "CUDA.Utils.reduceMin" )
    SCAI_LOG_INFO( logger, "minval for " << n << " elements " )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
    ValueType result = thrust::reduce( data, data + n, zero, thrust::minimum<ValueType>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "min of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

// Be careful: template<ResultType, ArgumentType>, but unary_function<ArgumentType, ResultType> 

template<typename RealType, typename ValueType>
struct absolute_value: public thrust::unary_function<ValueType, RealType>
{
    __host__ __device__
    RealType operator()( const ValueType& x ) const
    {
        return static_cast<RealType>( Math::abs( x ) );
    }
};

template<typename ValueType>
ValueType CUDAReduceUtils::reduceAbsMaxVal( const ValueType array[], const IndexType n, const ValueType zero )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    SCAI_REGION( "CUDA.Utils.reduceAbsMax" )

    SCAI_LOG_INFO( logger, "absMaxVal for " << n << " elements " )

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

    RealType result = thrust::transform_reduce(
                           data,
                           data + n,
                           absolute_value<RealType, ValueType>(),
                           RealType( zero ),
                           thrust::maximum<RealType>() );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "abs max of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAReduceUtils::reduce( const ValueType array[], const IndexType n, const ValueType zero, BinaryOp op )
{
    SCAI_LOG_INFO ( logger, "reduce # array = " << array << ", n = " << n << ", op = " << op )

    ValueType result;

    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    RealType redResult;
    RealType redZero = zero;
    const RealType* redArray = reinterpret_cast<const RealType*>( array );

    switch ( op )
    {
        case BinaryOp::ADD :
            result = reduceSum( array, n, zero );
            break;

        case BinaryOp::MAX :
            // SCAI_ASSERT_EQ_ERROR( common::TypeTraits<RealType>::stype, common::TypeTraits<ValueType>::stype, "MAX not supported for complex" )
            redResult = reduceMaxVal( redArray, n, redZero );
            result = redResult;
            break;

        case BinaryOp::MIN :
            redResult = reduceMinVal( redArray, n, redZero );
            result = redResult;
            break;

        case BinaryOp::ABS_MAX :
            result = reduceAbsMaxVal( array, n, zero );
            break;

        default:
            COMMON_THROWEXCEPTION( "Unsupported reduce op " << op )
    }

    return result;
}

/* --------------------------------------------------------------------------- */

template <typename ValueType>
ValueType CUDAReduceUtils::reduce2(
    const ValueType array1[],
    const ValueType array2[],
    const IndexType n,
    const BinaryOp binOp,
    const ValueType zero,
    const BinaryOp redOp )
{
    SCAI_REGION( "CUDA.Utils.reduce2" )

    SCAI_LOG_DEBUG( logger, "reduce2<" << TypeTraits<ValueType>::id() << ">, n = " << n )

    // Currently on CUDA: reduce operator requires temporary array in any case

    thrust::device_vector<ValueType> temp( n );

    ValueType* tmpData  = temp.data().get();

    CUDAUtils::binaryOp( tmpData, array1, array2, n, binOp );

    ValueType result = reduce( tmpData, n, zero, redOp );

    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAReduceUtils::scan( ValueType array[], const IndexType n, ValueType first, const bool exclusive, const bool append )
{
    using namespace thrust::placeholders;

    SCAI_REGION( "CUDA.Utils.scan" )

    SCAI_LOG_INFO( logger, "scan<" << TypeTraits<ValueType>::id() <<  ">, #n = " << n 
                            << ", first = " << first << ", exclusive = " << exclusive )

    SCAI_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> array_ptr( array );

    ValueType result = 0;

    if ( exclusive )
    {
        if ( append )
        {
            thrust::exclusive_scan( array_ptr, array_ptr + n + 1, array_ptr );

            if ( first != common::Constants::ZERO )
            {
                thrust::for_each( array_ptr, array_ptr + n + 1, _1 += first );
            }

            thrust::host_vector<ValueType> numValues( array_ptr + n, array_ptr + n + 1 );
            result = numValues[0];
        }
        else
        {
            COMMON_THROWEXCEPTION( "exclusive scan, no append is unsupported" )
        }
    }
    else
    {
        if ( append )
        {
            COMMON_THROWEXCEPTION( "inclusive scan, append is unsupported" )
        }
        else if ( n > 0 )
        {
            thrust::inclusive_scan( array_ptr, array_ptr + n, array_ptr );

            if ( first != common::Constants::ZERO )
            {
		        SCAI_LOG_INFO( logger, "now add first = " << first << ", n = " << n )
                thrust::for_each( array_ptr, array_ptr + n, _1 += first );
            }

            thrust::host_vector<ValueType> numValues( array_ptr + n - 1, array_ptr + n );
            result = numValues[0];
        }
    }

    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
__global__
void isSortedKernel( bool* result, const IndexType numValues, const ValueType* values, const CompareOp op )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numValues )
    {
        result[i] = compare( values[i], op, values[i + 1] );
    }
}

template<typename ValueType>
bool CUDAReduceUtils::isSorted( const ValueType array[], const IndexType n, const CompareOp op )
{
    SCAI_REGION( "CUDA.Utils.isSorted" )

    SCAI_LOG_INFO( logger, "isSorted<" << TypeTraits<ValueType>::id() << ">, n = " << n
                           << ", op = " << op )

    SCAI_CHECK_CUDA_ACCESS

    if ( n < 2 )
    {
        return true; // 0 or 1 element is always sorted
    }

    // create a tempory bool array on device with n-1 entries
    thrust::device_ptr<bool> resultPtr = thrust::device_malloc<bool>( n - 1 );
    bool* resultRawPtr = thrust::raw_pointer_cast( resultPtr );
    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n - 1, dimBlock.x );

    isSortedKernel<ValueType> <<< dimGrid, dimBlock>>> ( resultRawPtr, n - 1, array, op );

    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
    return thrust::reduce( resultPtr, resultPtr + n - 1, true, thrust::logical_and<bool>() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
struct InvalidIndex
{
    const ValueType size; //!< size of array for which index is checked

    InvalidIndex( ValueType _size ) : size( _size )
    {}

    __host__ __device__
    bool operator()( ValueType y )
    {
        return ! common::Utils::validIndex( y, size );
    }
};

/* --------------------------------------------------------------------------- */

bool CUDAReduceUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
    SCAI_REGION( "CUDA.Utils.validIndexes" )

    SCAI_LOG_DEBUG( logger, "validIndexes: array[" << n << "], size " << size )

    bool validFlag = true;

    if ( n > 0 )
    {
        SCAI_CHECK_CUDA_ACCESS
        thrust::device_ptr<IndexType> arrayPtr( const_cast<IndexType*> ( array ) );
        bool error = false;
        // count invalid indexes
        error = thrust::transform_reduce( arrayPtr,
                                          arrayPtr + n,
                                          InvalidIndex<IndexType>( size ),
                                          false,
                                          thrust::logical_or<bool>() );

        if ( error )
        {
            validFlag = false;
        }
    }

    return validFlag;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDAReduceUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag << "]" )
    // we keep the registrations for IndexType as we do not need conversions
    KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, ctx, flag );
}

template<typename ValueType>
void CUDAReduceUtils::RegArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;

    SCAI_LOG_DEBUG( logger, "registerV array UtilsKernel CUDA [" << flag
                    << "] --> ValueType = " << common::getScalarType<ValueType>() )

    // Note: these kernels will be instantiated for numeric types + IndexType

    KernelRegistry::set<UtilKernelTrait::reduce<ValueType> >( reduce, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::reduce2<ValueType> >( reduce2, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scan<ValueType> >( scan, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<ValueType> >( isSorted, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDAReduceUtils::CUDAReduceUtils()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAReduceUtils::~CUDAReduceUtils()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAReduceUtils CUDAReduceUtils::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
