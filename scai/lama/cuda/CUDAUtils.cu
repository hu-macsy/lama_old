/**
 * @file CUDAUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Implementation of CSR utilities with CUDA
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/cuda/CUDAUtils.hpp>

// local library
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/lama/UtilKernelTrait.hpp>

#include <scai/lama/cuda/CUDASettings.hpp>

// internal scai libraries
#include <scai/common/macros/assert.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/Constants.hpp>

// thrust
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>

// boost
#include <boost/preprocessor.hpp>

using namespace scai::common;

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( CUDAUtils::logger, "CUDA.Utils" )

/* ------------------------------------------------------------------------------------------------------------------ */
/*   Kernel used for scale, set, setScale                                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename T1, typename T2>
__global__
void setScaleKernel( T1* out, const T1 beta, const T2* in, IndexType n )
{
    // Note: out == in does not harm, also not for performance

    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = static_cast<T1>( in[i] * beta );
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
        return y >= size || y < 0;
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

    ValueType zero( - TypeTraits<ValueType>::getMax() );

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
        return x < ValueType( 0 ) ? -x : x;
        // return abs( x );
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
ValueType CUDAUtils::reduce( const ValueType array[], const IndexType n, common::reduction::ReductionOp op )
{
    SCAI_LOG_INFO ( logger, "reduce # array = " << array << ", n = " << n << ", op = " << op )

    ValueType result;

    switch ( op )
    {
        case common::reduction::ADD :
            result = reduceSum( array, n );
            break;
        case common::reduction::MAX :
            result = reduceMaxVal( array, n );
            break;
        case common::reduction::MIN :
            result = reduceMinVal( array, n );
            break;
        case common::reduction::ABS_MAX :
            result = reduceAbsMaxVal( array, n );
            break;
        default:
            COMMON_THROWEXCEPTION( "Unsupported reduce op " << op )
    }

    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const common::reduction::ReductionOp op )
{
    using namespace thrust::placeholders;

    SCAI_LOG_INFO( logger, "setVal # array = " << array << ", n = " << n << ", val = " << val << ", op = " << op )

    SCAI_CHECK_CUDA_ACCESS

    if ( n > 0 )
    {
        thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

        switch ( op ) 
        {
            case common::reduction::COPY:
                thrust::fill( data, data + n, val );
                break;
            case common::reduction::ADD:
                thrust::for_each( data, data + n,  _1 += val);
                break;
            case common::reduction::MULT:
                {
                    if ( val == scai::common::constants::ZERO )
                    {
                        thrust::fill( data, data + n, ValueType( 0 ) );
                    }
                    else
                    {
                        thrust::for_each( data, data + n,  _1 *= val);
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
            result[i] = values[i] <= values[i + 1];
        }
        else
        {
            result[i] = values[i] >= values[i + 1];
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
void gatherKernel( ValueType1* out, const ValueType2* in, const IndexType* indexes, IndexType n )
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

template<typename T1, typename T2>
__global__
void scatter_kernel( T1* out, const IndexType* indexes, const T2* in, IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[indexes[i]] = in[i];
    }
}

template<typename ValueType1, typename ValueType2>
void CUDAUtils::setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const IndexType n )
{
    SCAI_LOG_INFO( logger,
                   "setScatter<" << TypeTraits<ValueType1>::id() << "," << TypeTraits<ValueType2>::id() << ">( ..., n = " << n << ")" )

    if ( n > 0 )
    {
        SCAI_CHECK_CUDA_ACCESS

        const int blockSize = 256;
        dim3 dimBlock( blockSize, 1, 1 );
        dim3 dimGrid = makeGrid( n, dimBlock.x );
    
        scatter_kernel <<< dimGrid, dimBlock>>>( out, indexes, in, n );
    
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    }
}

/* --------------------------------------------------------------------------- */

template<typename T1, typename T2>
__global__
void setKernelCopy( T1* out, const T2* in, IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = static_cast<T1>( in[i] );
    }
}

template<typename T1, typename T2>
__global__
void setKernelAdd( T1* out, const T2* in, IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] += static_cast<T1>( in[i] );
    }
}

template<typename T1, typename T2>
__global__
void setKernelMult( T1* out, const T2* in, IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] *= static_cast<T1>( in[i] );
    }
}

template<typename ValueType1, typename ValueType2>
void CUDAUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n, const common::reduction::ReductionOp op )
{
    SCAI_LOG_INFO( logger,
                   "set<" << TypeTraits<ValueType1>::id() << "," << TypeTraits<ValueType2>::id() << ">( ..., n = " << n << ")" )

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
        case common::reduction::COPY :
            setKernelCopy <<< dimGrid, dimBlock>>>( out, in, n );
            break;
        case common::reduction::ADD :
            setKernelAdd <<< dimGrid, dimBlock>>>( out, in, n );
            break;
        case common::reduction::MULT :
            setKernelMult <<< dimGrid, dimBlock>>>( out, in, n );
            break;
         default:
            COMMON_THROWEXCEPTION( "Unsupported reduction op " << op )
    }

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void CUDAUtils::setScale( ValueType1 out[],
                          const ValueType1 beta,
                          const ValueType2 in[], const IndexType n )
{
    SCAI_LOG_INFO( logger,
                   "set<" << TypeTraits<ValueType1>::id() << "," << TypeTraits<ValueType2>::id() << ">( ..., n = " << n << ")" )

    SCAI_LOG_DEBUG( logger, "out = " << out << ", in = " << in )

    if ( n <= 0 )
    {
        return;
    }

    if ( beta == scai::common::constants::ZERO )
    {
        // in array might be undefined

        setVal( out, n, beta, common::reduction::COPY );
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
__global__
void invertVectorComponents_kernel( ValueType* array, IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    ValueType one = 1.0;

    if ( i < n )
    {
        array[i] = one / array[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::invert( ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger,
                   "invert Vector components for vector of type " << TypeTraits<ValueType>::id() << " and size n = " << n << "." )

    if ( n <= 0 )
    {
        return;
    }

    SCAI_CHECK_CUDA_ACCESS

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    invertVectorComponents_kernel <<< dimGrid, dimBlock>>>( array, n );
    cudaStreamSynchronize( 0 );
    SCAI_CHECK_CUDA_ERROR
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDAUtils::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;
    using common::context::CUDA;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // lower priority

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

    SCAI_LOG_INFO( logger, "set general utilty routines for CUDA in Interface" )

    KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, CUDA, flag );
    KernelRegistry::set<UtilKernelTrait::reduce<IndexType> >( reduce, CUDA, flag );

    KernelRegistry::set<UtilKernelTrait::setVal<IndexType> >( setVal, CUDA, flag );
    KernelRegistry::set<UtilKernelTrait::setOrder<IndexType> >( setOrder, CUDA, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<IndexType> >( getValue, CUDA, flag );

    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<IndexType> >( absMaxDiffVal, CUDA, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<IndexType> >( isSorted, CUDA, flag );

    KernelRegistry::set<UtilKernelTrait::setScatter<IndexType, IndexType> >( setScatter, CUDA, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<IndexType, IndexType> >( setGather, CUDA, flag );
    KernelRegistry::set<UtilKernelTrait::set<IndexType, IndexType> >( set, CUDA, flag );

#define LAMA_UTILS2_REGISTER(z, J, TYPE )                                                                        \
    KernelRegistry::set<UtilKernelTrait::setScale<TYPE, ARITHMETIC_CUDA_TYPE_##J> >( setScale, CUDA, flag );     \
    KernelRegistry::set<UtilKernelTrait::setGather<TYPE, ARITHMETIC_CUDA_TYPE_##J> >( setGather, CUDA, flag );   \
    KernelRegistry::set<UtilKernelTrait::setScatter<TYPE, ARITHMETIC_CUDA_TYPE_##J> >( setScatter, CUDA, flag ); \
    KernelRegistry::set<UtilKernelTrait::set<TYPE, ARITHMETIC_CUDA_TYPE_##J> >( set, CUDA, flag );               \
     
#define LAMA_UTILS_REGISTER(z, I, _)                                                                             \
    KernelRegistry::set<UtilKernelTrait::reduce<ARITHMETIC_CUDA_TYPE_##I> >( reduce, CUDA, flag );               \
    KernelRegistry::set<UtilKernelTrait::setVal<ARITHMETIC_CUDA_TYPE_##I> >( setVal, CUDA, flag );               \
    KernelRegistry::set<UtilKernelTrait::setOrder<ARITHMETIC_CUDA_TYPE_##I> >( setOrder, CUDA, flag );           \
    KernelRegistry::set<UtilKernelTrait::getValue<ARITHMETIC_CUDA_TYPE_##I> >( getValue, CUDA, flag );           \
    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<ARITHMETIC_CUDA_TYPE_##I> >( absMaxDiffVal, CUDA, flag ); \
    KernelRegistry::set<UtilKernelTrait::isSorted<ARITHMETIC_CUDA_TYPE_##I> >( isSorted, CUDA, flag );           \
    KernelRegistry::set<UtilKernelTrait::invert<ARITHMETIC_CUDA_TYPE_##I> >( invert, CUDA, flag );               \
    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT,                                                                   \
                     LAMA_UTILS2_REGISTER,                                                                       \
                     ARITHMETIC_CUDA_TYPE_##I )                                                                  \
     
    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT, LAMA_UTILS_REGISTER, _ )

#undef LAMA_UTILS_REGISTER
#undef LAMA_UTILS2_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDAUtils::CUDAUtils()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

CUDAUtils::~CUDAUtils()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

CUDAUtils CUDAUtils::guard;    // guard variable for registration

} /* end namespace lama */

} /* end namespace scai */
