/**
 * @file CUDAUtils.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAUtils.hpp>

#include <lama/exception/LAMAAssert.hpp>

// thrust
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>

// others
#include <typeinfo>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDAUtils::logger, "CUDA.Utils" )

/* ------------------------------------------------------------------------------------------------------------------ */
/*                                                  scale                                                             */
/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void CUDAUtils::scale( ValueType *mValues, const IndexType n, const OtherValueType value )
{
    LAMA_LOG_INFO( logger, "scale, #n = " << n << ", value = " << value )

    LAMA_CHECK_CUDA_ACCESS

    ValueType castedValue = static_cast<ValueType>( value );

    thrust::device_ptr<ValueType> mValuesPtr( const_cast<ValueType*>( mValues ) );
    thrust::constant_iterator<ValueType> valueItr( castedValue );

    thrust::transform( mValuesPtr, mValuesPtr + n, valueItr, mValuesPtr, thrust::multiplies<ValueType>() );
}

/* --------------------------------------------------------------------------- */

template<typename T>
struct InvalidIndex
{
    const T size;  //!< size of array for which index is checked

    InvalidIndex( T _size ) : size( _size ) {}

    __host__ __device__
    bool operator()( T y )
    {
        return y >= size || y < 0;
    }
};

/* --------------------------------------------------------------------------- */

bool CUDAUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
    LAMA_LOG_DEBUG( logger, "validIndexes: array[" << n << "], size " << size )

    bool validFlag = true;

    if ( n > 0 ) 
    {
        LAMA_CHECK_CUDA_ACCESS

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
ValueType CUDAUtils::sum( const ValueType array[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "sum # array = " << array << ", n = " << n )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

    ValueType zero = static_cast<ValueType>( 0 );

    ValueType result = thrust::reduce( data, data + n, zero, thrust::plus<ValueType>() );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );

    LAMA_LOG_INFO( logger, "sum of " << n << " values = " << result )

    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::setVal( ValueType array[], const IndexType n, const ValueType val )
{
    LAMA_LOG_INFO( logger, "setVal # array = " << array << ", n = " << n << ", val = " << val )

    LAMA_CHECK_CUDA_ACCESS

    if ( n > 0 )
    {
        thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
        thrust::fill( data, data + n, val );

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDAUtils::setOrder( ValueType array[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "setOrder # array = " << array << ", n = " << n )
    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> array_ptr( const_cast<ValueType*>( array ) );
    thrust::sequence( array_ptr, array_ptr + n );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::getValue( const ValueType* array, const IndexType i )
{
    LAMA_LOG_INFO( logger, "getValue # i = " << i )
    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> arrayPtr( const_cast<ValueType*>( array ) );
    thrust::host_vector<ValueType> arrayHost( arrayPtr + i, arrayPtr + i + 1 );

    return arrayHost[0];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::maxval( const ValueType array[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "maxval for " << n << " elements " )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );
    ValueType zero = static_cast<ValueType>( 0 );
    ValueType result = thrust::reduce( data, data + n, zero, thrust::maximum<ValueType>() );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );

    LAMA_LOG_INFO( logger, "max of " << n << " values = " << result )

    return result;
}

/* --------------------------------------------------------------------------- */

template<typename T>
struct absolute_value: public thrust::unary_function<T,T>
{
    __host__ __device__
    T operator()( const T &x ) const
    {
        return x < T( 0 ) ? -x : x;
    }
};

template<typename ValueType>
ValueType CUDAUtils::absMaxVal( const ValueType array[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "absMaxVal for " << n << " elements " )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

    ValueType zero = static_cast<ValueType>( 0 );

    ValueType result = thrust::transform_reduce( data, data + n, absolute_value<ValueType>(), zero,
                       thrust::maximum<ValueType>() );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );

    LAMA_LOG_INFO( logger, "abs max of " << n << " values = " << result )

    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CUDAUtils::absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "absMaxDiffVal for " << n << " elements " )

    LAMA_CHECK_CUDA_ACCESS

    thrust::device_ptr<ValueType> data1( const_cast<ValueType*>( array1 ) );
    thrust::device_ptr<ValueType> data2( const_cast<ValueType*>( array2 ) );

    thrust::device_vector<ValueType> temp( n );

    // compute temp =  array1 - array2

    thrust::transform( data1, data1 + n, data2, temp.begin(), thrust::minus<ValueType>() );

    ValueType zero = static_cast<ValueType>( 0 );

    ValueType result = thrust::transform_reduce( temp.begin(), temp.end(), absolute_value<ValueType>(), zero,
                       thrust::maximum<ValueType>() );

    /* Not available, but would be useful:

     ValueType result = thrust::transform_reduce( data1, data1 + n,
     data2,
     thrust::minus<ValueType>(),
     zero,
     thrust::maximum<ValueType>());
     */

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" )

    LAMA_LOG_INFO( logger, "abs max diff of " << n << " values = " << result )

    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
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

template<typename ValueType1,typename ValueType2>
void CUDAUtils::setGather( ValueType1 out[], const ValueType2 in[], const IndexType indexes[], const IndexType n )
{
    LAMA_LOG_INFO( logger,
                   "setGather<" << typeid(ValueType1).name() << "," << typeid(ValueType2).name() << ">( ..., n = " << n << ")" )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    gatherKernel<<<dimGrid,dimBlock>>>( out, in, indexes, n );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename T1,typename T2>
__global__
void scatter_kernel( T1* out, const IndexType* indexes, const T2* in, IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[indexes[i]] = in[i];
    }
}

template<typename ValueType1,typename ValueType2>
void CUDAUtils::setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const IndexType n )
{
    LAMA_LOG_INFO( logger,
                   "setScatter<" << typeid(ValueType1).name() << "," << typeid(ValueType2).name() << ">( ..., n = " << n << ")" )

    LAMA_CHECK_CUDA_ACCESS

    const int block_size = 256;
    dim3 dimBlock( block_size, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    scatter_kernel<<<dimGrid,dimBlock>>>( out, indexes, in, n );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename T1,typename T2>
__global__
void setKernel( T1* out, const T2* in, IndexType n )
{
    const IndexType i = threadId( gridDim, blockIdx, blockDim, threadIdx );

    if ( i < n )
    {
        out[i] = static_cast<T1>( in[i] );
    }
}

template<typename ValueType1,typename ValueType2>
void CUDAUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n )
{
    LAMA_LOG_INFO( logger,
                   "set<" << typeid(ValueType1).name() << "," << typeid(ValueType2).name() << ">( ..., n = " << n << ")" )

    LAMA_LOG_DEBUG( logger, "out = " << out << ", in = " << in )

    LAMA_CHECK_CUDA_ACCESS

    if ( n > 0 )
    {
        const int block_size = 256;
        dim3 dimBlock( block_size, 1, 1 );
        dim3 dimGrid = makeGrid( n, dimBlock.x );

        setKernel<<<dimGrid,dimBlock>>>( out, in, n );

        LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
    }
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
    LAMA_LOG_INFO( logger,
                   "invert Vector components for vector of type " << typeid(ValueType).name() << " and size n = " << n << "." )

    LAMA_CHECK_CUDA_ACCESS

    if ( n > 0 )
    {
        const int block_size = 256;
        dim3 dimBlock( block_size, 1, 1 );
        dim3 dimGrid = makeGrid( n, dimBlock.x );

        invertVectorComponents_kernel<<<dimGrid,dimBlock>>>( array, n );
        cudaStreamSynchronize( 0 );
        LAMA_CHECK_CUDA_ERROR
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDAUtils::setInterface( UtilsInterface& Utils )
{
    LAMA_INTERFACE_REGISTER( Utils, validIndexes )

    LAMA_INTERFACE_REGISTER_TT( Utils, scale, float, float )
    LAMA_INTERFACE_REGISTER_TT( Utils, scale, double, float )
    LAMA_INTERFACE_REGISTER_TT( Utils, scale, float, double )
    LAMA_INTERFACE_REGISTER_TT( Utils, scale, double, double )

    LAMA_INTERFACE_REGISTER_T( Utils, sum, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, sum, float )
    LAMA_INTERFACE_REGISTER_T( Utils, sum, double )

    LAMA_INTERFACE_REGISTER_T( Utils, setVal, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, setVal, float )
    LAMA_INTERFACE_REGISTER_T( Utils, setVal, double )

    LAMA_INTERFACE_REGISTER_T( Utils, setOrder, IndexType )

    LAMA_INTERFACE_REGISTER_T( Utils, getValue, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, getValue, float )
    LAMA_INTERFACE_REGISTER_T( Utils, getValue, double )

    LAMA_INTERFACE_REGISTER_T( Utils, maxval, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, maxval, float )
    LAMA_INTERFACE_REGISTER_T( Utils, maxval, double )

    LAMA_INTERFACE_REGISTER_T( Utils, absMaxVal, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxVal, float )
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxVal, double )

    LAMA_INTERFACE_REGISTER_T( Utils, absMaxDiffVal, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxDiffVal, float )
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxDiffVal, double )

    LAMA_INTERFACE_REGISTER_TT( Utils, set, int, int )
    LAMA_INTERFACE_REGISTER_TT( Utils, set, float, float )
    LAMA_INTERFACE_REGISTER_TT( Utils, set, float, double )
    LAMA_INTERFACE_REGISTER_TT( Utils, set, double, float )
    LAMA_INTERFACE_REGISTER_TT( Utils, set, double, double )

    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, int, int )
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, float, float )
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, float, double )
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, double, float )
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, double, double )

    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, int, int )
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, float, float )
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, float, double )
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, double, float )
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, double, double )

    LAMA_INTERFACE_REGISTER_T( Utils, invert, float )
    LAMA_INTERFACE_REGISTER_T( Utils, invert, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool CUDAUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::CUDA );
    setInterface( interface.Utils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool CUDAUtils::initialized = registerInterface();

} // namespace lama
