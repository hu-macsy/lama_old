/**
 * @file utilskernel/cuda/CUDAUtils.cu
 *
 * @license
 * Copyright (c) 2009-2016
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

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>


// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/cuda/CUDASettings.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/Math.hpp>

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

SCAI_LOG_DEF_LOGGER( CUDAUtils::logger, "CUDA.Utils" )

/* --------------------------------------------------------------------------- */
/*                                 elementwise kernel                          */
/* --------------------------------------------------------------------------- */

/* invert / reciprocal */

template<typename ValueType>
__global__
void invertVectorComponents_kernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );
    ValueType one = 1.0;

    if ( i < n )
    {
        array[i] = one / array[i];
    }
}

/* conj */

template<typename ValueType>
__global__
void conjKernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array[i] = Math::conj( array[i] );
    }
}

/* exp */

template<typename ValueType>
__global__
void expKernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array[i] = Math::exp( array[i] );
    }
}

/* sqrt */

template<typename ValueType>
__global__
void sqrtKernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array[i] = Math::sqrt( array[i] );
    }
}

/* sin */

template<typename ValueType>
__global__
void sinKernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array[i] = Math::sin( array[i] );
    }
}

/* cos */

template<typename ValueType>
__global__
void cosKernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array[i] = Math::cos( array[i] );
    }
}

/* tan */

template<typename ValueType>
__global__
void tanKernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array[i] = Math::tan( array[i] );
    }
}

/* atan */

template<typename ValueType>
__global__
void atanKernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array[i] = Math::atan( array[i] );
    }
}

/* log */

template<typename ValueType>
__global__
void logKernel( ValueType* array, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        array[i] = Math::log( array[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                            copysign                                                                */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void copysignKernel( ValueType* result, const ValueType* x, const ValueType* y, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        result[i] = copysign( x[i], y[i] ); // copysign not from Math.hpp, CUDA has own implementation
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                            vectorScale                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
__global__
void vectorScaleKernel( ValueType* result, const ValueType* x, const ValueType* y, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        result[i] = x[i] * y[i];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*   Kernel used for scale, set, setScale                                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, typename OtherValueType>
__global__
void setScaleKernel( ValueType* out, const ValueType beta, const OtherValueType* in, const IndexType n )
{
    // Note: out == in does not harm, also not for performance
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = beta * static_cast<ValueType>( in[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  scale                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDAUtils::scale( ValueType* values, const ValueType scale, const IndexType n )
{
    SCAI_LOG_INFO( logger, "scale, #n = " << n << ", scale = " << scale )

    if ( n == 0 )
    {
        return;
    }

    if ( scale == scai::common::constants::ZERO )
    {
        setVal( values, n, scale );
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );
    // there is no performance loss in using same kernel as setScale
    // kernel fits well even if in and out are aliased
    setScaleKernel <<< dimGrid, dimBlock>>>( values, scale, values, n );
}

template<typename ValueType>
void CUDAUtils::vectorScale( ValueType* result, const ValueType* x, const ValueType* y, const IndexType n )
{
    SCAI_LOG_INFO( logger, "vectorScale, #n = " << n )

    if ( n == 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );
    vectorScaleKernel <<< dimGrid, dimBlock>>>( result, x, y, n );
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

bool CUDAUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
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

template<typename ValueType>
ValueType CUDAUtils::reduceSum( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "sum # array = " << array << ", n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
    ValueType zero = ValueType( 0 );
    ValueType result = thrust::reduce( data, data + n, zero, thrust::plus<ValueType>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "sum of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::reduceMaxVal( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "maxval for " << n << " elements " )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
    ValueType zero( TypeTraits<ValueType>::getMin() );
    ValueType result = thrust::reduce( data, data + n, zero, thrust::maximum<ValueType>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "max of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::reduceMinVal( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "minval for " << n << " elements " )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
    ValueType zero( TypeTraits<ValueType>::getMax() );
    ValueType result = thrust::reduce( data, data + n, zero, thrust::minimum<ValueType>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "min of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
struct absolute_value: public thrust::unary_function<ValueType, ValueType>
{
    __host__ __device__
    ValueType operator()( const ValueType& x ) const
    {
        // return x < ValueType( 0 ) ? -x : x;
        return Math::abs( x );
    }
};

template<typename ValueType>
ValueType CUDAUtils::reduceAbsMaxVal( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "absMaxVal for " << n << " elements " )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
    ValueType zero( 0 );
    ValueType result = thrust::transform_reduce( data, data + n, absolute_value<ValueType>(), zero,
                       thrust::maximum<ValueType>() );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_LOG_INFO( logger, "abs max of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::reduce( const ValueType array[], const IndexType n, reduction::ReductionOp op )
{
    SCAI_LOG_INFO ( logger, "reduce # array = " << array << ", n = " << n << ", op = " << op )
    ValueType result;

    switch ( op )
    {
        case reduction::ADD :
            result = reduceSum( array, n );
            break;

        case reduction::MAX :
            result = reduceMaxVal( array, n );
            break;

        case reduction::MIN :
            result = reduceMinVal( array, n );
            break;

        case reduction::ABS_MAX :
            result = reduceAbsMaxVal( array, n );
            break;

        default:
            COMMON_THROWEXCEPTION( "Unsupported reduce op " << op )
    }

    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const reduction::ReductionOp op )
{
    using namespace thrust::placeholders;
    SCAI_LOG_INFO( logger, "setVal # array = " << array << ", n = " << n << ", val = " << val << ", op = " << op )
    SCAI_CHECK_CUDA_ACCESS

    if ( n > 0 )
    {
        thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
        ValueType value = static_cast<ValueType>( val );

        switch ( op )
        {
            case reduction::COPY:
                thrust::fill( data, data + n, value );
                break;

            case reduction::ADD:
                thrust::for_each( data, data + n,  _1 += value );
                break;

            case reduction::SUB:
                thrust::for_each( data, data + n,  _1 -= value );
                break;

            case reduction::MULT:
            {
                if ( val == scai::common::constants::ZERO )
                {
                    thrust::fill( data, data + n, ValueType( 0 ) );
                }
                else
                {
                    thrust::for_each( data, data + n,  _1 *= value );
                }
            }
            break;

            case reduction::DIVIDE:
            {
                if ( val == scai::common::constants::ZERO )
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
                COMMON_THROWEXCEPTION( "unsupported reduction op: " << op )
        }

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::setOrder( ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "setOrder # array = " << array << ", n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> array_ptr( const_cast<ValueType*>( array ) );
    thrust::sequence( array_ptr, array_ptr + n );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n )
{
    SCAI_LOG_INFO( logger, "setSequence # array = " << array << ", n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> array_ptr( const_cast<ValueType*>( array ) );
    thrust::sequence( thrust::device, array_ptr, array_ptr + n, startValue, inc );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::getValue( const ValueType* array, const IndexType i )
{
    SCAI_LOG_INFO( logger, "getValue # i = " << i )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> arrayPtr( const_cast<ValueType*>( array ) );
    thrust::host_vector<ValueType> arrayHost( arrayPtr + i, arrayPtr + i + 1 );
    return arrayHost[0];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "absMaxDiffVal for " << n << " elements " )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> data1( const_cast<ValueType*>( array1 ) );
    thrust::device_ptr<ValueType> data2( const_cast<ValueType*>( array2 ) );
    thrust::device_vector<ValueType> temp( n );
    // compute temp =  array1 - array2
    thrust::transform( data1, data1 + n, data2, temp.begin(), thrust::minus<ValueType>() );
    ValueType result = thrust::transform_reduce( temp.begin(), temp.end(), absolute_value<ValueType>(), static_cast<ValueType>( 0.0 ),
                       thrust::maximum<ValueType>() );
    /* Not available, but would be useful:

     ValueType result = thrust::transform_reduce( data1, data1 + n,
     data2,
     thrust::minus<ValueType>(),
     zero,
     thrust::maximum<ValueType>());
     */
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" )
    SCAI_LOG_INFO( logger, "abs max diff of " << n << " values = " << result )
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::copysign( ValueType result[], const ValueType x[], const ValueType y[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "copysign<" << TypeTraits<ValueType>::id() << ">( ..., n = " << n << ")" )
    SCAI_LOG_DEBUG( logger, "result = " << result << ", x = " << x << ", y = " << y )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );
    copysignKernel<ValueType, true> <<< dimGrid, dimBlock>>> ( result, x, y, n );
}

/* --------------------------------------------------------------------------- */

// template argument ascending: make two instantiations of kernel to avoid bool test
template<typename ValueType, bool ascending>
__global__
void isSortedKernel( bool* result, const IndexType numValues, const ValueType* values )
{
    const int i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < numValues )
    {
        if ( ascending )
        {
            // not possible, <= not defined on complex
            // ToDo: warp divergence possible?
//            result[i] = values[i] <= values[i + 1];
            result[i] = values[i] < values[i + 1] || values[i] == values[i + 1];
        }
        else
        {
//            result[i] = values[i] >= values[i + 1];
            result[i] = values[i] > values[i + 1] || values[i] == values[i + 1];
        }
    }
}

template<typename ValueType>
bool CUDAUtils::isSorted( const ValueType array[], const IndexType n, bool ascending )
{
    SCAI_LOG_INFO( logger, "isSorted<" << TypeTraits<ValueType>::id()
                   << ">, n = " << n << ", ascending = " << ascending )
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

    if ( ascending )
    {
        isSortedKernel<ValueType, true> <<< dimGrid, dimBlock>>> ( resultRawPtr, n - 1, array );
    }
    else
    {
        isSortedKernel<ValueType, false> <<< dimGrid, dimBlock>>> ( resultRawPtr, n - 1, array );
    }

    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
    return thrust::reduce( resultPtr, resultPtr + n - 1, true, thrust::logical_and<bool>() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
__global__
void gatherKernel( ValueType1* out, const ValueType2* in, const IndexType* indexes, const IndexType n )
{
    // Kernel also supports implicit type conversions
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = static_cast<ValueType1>( in[indexes[i]] );
    }
}

template<typename ValueType1, typename ValueType2>
void CUDAUtils::setGather( ValueType1 out[], const ValueType2 in[], const IndexType indexes[], const IndexType n )
{
    SCAI_LOG_INFO( logger,
                   "setGather<" << TypeTraits<ValueType1>::id() << "," << TypeTraits<ValueType2>::id() << ">( ..., n = " << n << ")" )
    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );
    gatherKernel <<< dimGrid, dimBlock>>>( out, in, indexes, n );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, typename OtherValueType>
__global__
void scatter_kernel( ValueType* out, const IndexType* indexes, const OtherValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[indexes[i]] = static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename OtherValueType>
__global__
void scatter_add_kernel( ValueType* out, const IndexType* indexes, const OtherValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[indexes[i]] += static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType1, typename ValueType2>
void CUDAUtils::setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const reduction::ReductionOp op, const IndexType n )
{
    SCAI_LOG_INFO( logger,
                   "setScatter<" << TypeTraits<ValueType1>::id() << "," << TypeTraits<ValueType2>::id() << ">( ..., n = " << n << ")" )

    if ( n > 0 )
    {
        SCAI_CHECK_CUDA_ACCESS
        const int blockSize = 256;
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( n, dimBlock.x );
     
        if ( op == reduction::COPY )
        {
            scatter_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n );
        }
        else if ( op == reduction::ADD )
        {
            scatter_add_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n );
        }

        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    }
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

template<typename ValueType, typename OtherValueType>
__global__
void setKernelCopy( ValueType* out, const OtherValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename OtherValueType>
__global__
void setKernelAdd( ValueType* out, const OtherValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] += static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename OtherValueType>
__global__
void setKernelSub( ValueType* out, const OtherValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] -= static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename OtherValueType>
__global__
void setKernelMult( ValueType* out, const OtherValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] *= static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename OtherValueType>
__global__
void setKernelDivide( ValueType* out, const OtherValueType* in, const IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] /= static_cast<ValueType>( in[i] );
    }
}

template<typename ValueType, typename OtherValueType>
void CUDAUtils::set( ValueType out[], const OtherValueType in[], const IndexType n, const reduction::ReductionOp op )
{
    SCAI_LOG_INFO( logger,
                   "set<" << TypeTraits<ValueType>::id() << "," << TypeTraits<OtherValueType>::id() << ">( ..., n = " << n << ")" )
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
        case reduction::COPY :
            setKernelCopy <<< dimGrid, dimBlock>>>( out, in, n );
            break;

        case reduction::ADD :
            setKernelAdd <<< dimGrid, dimBlock>>>( out, in, n );
            break;

        case reduction::SUB :
            setKernelSub <<< dimGrid, dimBlock>>>( out, in, n );
            break;

        case reduction::MULT :
            setKernelMult <<< dimGrid, dimBlock>>>( out, in, n );
            break;

        case reduction::DIVIDE :
            setKernelDivide <<< dimGrid, dimBlock>>>( out, in, n );
            break;

        default:
            COMMON_THROWEXCEPTION( "Unsupported reduction op " << op )
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::execElementwise( ValueType array[], const IndexType n, const elementwise::ElementwiseOp op )
{
    SCAI_LOG_INFO( logger, "execElementwise<" << TypeTraits<ValueType>::id() << ">( ..., n = " << n << ")" )
    SCAI_LOG_DEBUG( logger, "array = " << array )

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
        case elementwise::INVERT :
        {
            invertVectorComponents_kernel <<< dimGrid, dimBlock>>>( array, n );

            break;
        }

        case elementwise::CONJ :
        {
            conjKernel <<< dimGrid, dimBlock>>>( array, n );
        
            break;
        }

        case elementwise::EXP :
        {
            expKernel <<< dimGrid, dimBlock>>>( array, n );
        
            break;
        }

        case elementwise::SQRT :
        {
            sqrtKernel <<< dimGrid, dimBlock>>>( array, n );
        
            break;
        }

        case elementwise::SIN :
        {
            sinKernel <<< dimGrid, dimBlock>>>( array, n );
        
            break;
        }

        case elementwise::COS :
        {
            cosKernel <<< dimGrid, dimBlock>>>( array, n );
        
            break;
        }

        case elementwise::TAN :
        {
            tanKernel <<< dimGrid, dimBlock>>>( array, n );
        
            break;
        }

        case elementwise::ATAN :
        {
            atanKernel <<< dimGrid, dimBlock>>>( array, n );
        
            break;
        }

        case elementwise::LOG :
        {
            logKernel <<< dimGrid, dimBlock>>>( array, n );
        
            break;
        }

        default:
            COMMON_THROWEXCEPTION( "Unsupported reduction op " << op )
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    SCAI_CHECK_CUDA_ERROR
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, typename OtherValueType>
void CUDAUtils::setScale( ValueType out[],
                          const ValueType beta,
                          const OtherValueType in[], const IndexType n )
{
    SCAI_LOG_INFO( logger,
                   "set<" << TypeTraits<ValueType>::id() << "," << TypeTraits<OtherValueType>::id() << ">( ..., n = " << n << ")" )
    SCAI_LOG_DEBUG( logger, "out = " << out << ", in = " << in )

    if ( n <= 0 )
    {
        return;
    }

    if ( beta == scai::common::constants::ZERO )
    {
        // in array might be undefined
        setVal( out, n, beta, reduction::COPY );
        return;
    }

    SCAI_CHECK_CUDA_ACCESS
    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );
    setScaleKernel <<< dimGrid, dimBlock>>>( out, beta, in, n );
    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::scan( ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "scan<" << TypeTraits<ValueType>::id() <<  ">, #n = " << n )
    SCAI_CHECK_CUDA_ACCESS
    thrust::device_ptr<ValueType> array_ptr( array );
    thrust::exclusive_scan( array_ptr, array_ptr + n + 1, array_ptr );
    thrust::host_vector<ValueType> numValues( array_ptr + n, array_ptr + n + 1 );
    return numValues[0];
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CUDAUtils::sort( ValueType array[], IndexType perm[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "sort " << n << " values" )

    if ( n > 1 )
    {
        SCAI_CHECK_CUDA_ACCESS
        thrust::device_ptr<ValueType> array_d( array );
        thrust::device_ptr<IndexType> perm_d( perm );
        thrust::sequence( perm_d, perm_d + n );
        // stable sort, descending order, so override default comparison
        thrust::stable_sort_by_key( array_d, array_d + n, perm_d, thrust::greater<ValueType>() );
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "Utils: synchronize for sort FAILED" )
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDAUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag << "]" )
    // we keep the registrations for IndexType as we do not need conversions
    KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, ctx, flag );
}

template<typename ValueType>
void CUDAUtils::RegArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;

    SCAI_LOG_DEBUG( logger, "registerV array UtilsKernel CUDA [" << flag
                   << "] --> ValueType = " << common::getScalarType<ValueType>() )

    // Note: these kernels will be instantiated for numeric types + IndexType

    KernelRegistry::set<UtilKernelTrait::reduce<ValueType> >( reduce, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setOrder<ValueType> >( setOrder, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setSequence<ValueType> >( setSequence, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<ValueType> >( getValue, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<ValueType> >( absMaxDiffVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<ValueType> >( isSorted, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setVal<ValueType> >( setVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scan<ValueType> >( scan, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::sort<ValueType> >( sort, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::vectorScale<ValueType> >( vectorScale, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scatterVal<ValueType> >( scatterVal, ctx, flag );
}

template<typename ValueType>
void CUDAUtils::RegNumericKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;

    SCAI_LOG_DEBUG( logger, "registerV numeric UtilsKernel CUDA [" << flag
                   << "] --> ValueType = " << common::getScalarType<ValueType>() )

    KernelRegistry::set<UtilKernelTrait::execElementwise<ValueType> >( execElementwise, ctx, flag );
    //KernelRegistry::set<UtilKernelTrait::pow<ValueType> >( pow, ctx, flag );
    //KernelRegistry::set<UtilKernelTrait::powBase<ValueType> >( powBase, ctx, flag );
    //KernelRegistry::set<UtilKernelTrait::powExp<ValueType> >( powExp, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void CUDAUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_DEBUG( logger, "registerVO UtilsKernel CUDA [" << flag 
                     << "] --> ValueType = " << common::getScalarType<ValueType>() 
                     << ", OtherValueType = " << common::getScalarType<OtherValueType>() )
    KernelRegistry::set<UtilKernelTrait::setScale<ValueType, OtherValueType> >( setScale, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<ValueType, OtherValueType> >( setGather, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setScatter<ValueType, OtherValueType> >( setScatter, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::set<ValueType, OtherValueType> >( set, ctx, flag );
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
    kregistry::mepr::RegistratorV<RegNumericKernels, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARRAY_TYPES_CUDA_LIST, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAUtils::~CUDAUtils()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel CUDA version" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegNumericKernels, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARRAY_TYPES_CUDA_LIST, SCAI_ARRAY_TYPES_CUDA_LIST>::registerKernels( flag );
}

CUDAUtils CUDAUtils::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
