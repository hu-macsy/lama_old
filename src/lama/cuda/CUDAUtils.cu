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

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/cuda/utils.cu.h>
#include <lama/cuda/CUDAError.hpp>
#include <lama/cuda/CUDAUtils.hpp>
#include <lama/cuda/CUDASettings.hpp>

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

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDAUtils::logger, "CUDA.Utils" )

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
void CUDAUtils::scale( ValueType *values, const ValueType scale, const IndexType n )
{
    LAMA_LOG_INFO( logger, "scale, #n = " << n << ", scale = " << scale )

    if ( n == 0 )
    {
        return;
    }

    if ( scale == static_cast<ValueType>( 0 ) )
    {
        setVal( values, n, scale );  
    }

    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    // there is no performance loss in using same kernel as setScale
    // kernel fits well even if in and out are aliased

    setScaleKernel<<<dimGrid, dimBlock>>>( values, scale, values, n );
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
        // return x < T( 0 ) ? -x : x;
        return abs( x );
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
            result[i] = values[i] <= values[i+1];
        }
        else
        {
            result[i] = values[i] >= values[i+1];
        }
    }
}

template<typename ValueType>
bool CUDAUtils::isSorted( const ValueType array[], const IndexType n, bool ascending )
{
    LAMA_LOG_INFO( logger, "isSorted<" << Scalar::getType<ValueType>() 
                           << ">, n = " << n << ", ascending = " << ascending )

    LAMA_CHECK_CUDA_ACCESS

    if ( n < 2 )
    {
        return true;   // 0 or 1 element is always sorted
    }

    // create a tempory bool array on device with n-1 entries

    thrust::device_ptr<bool> resultPtr = thrust::device_malloc<bool>( n - 1 );

    bool* resultRawPtr = thrust::raw_pointer_cast( resultPtr );

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n - 1, dimBlock.x );

    if ( ascending )
    {
        isSortedKernel<ValueType, true><<<dimGrid, dimBlock>>> ( resultRawPtr, n - 1, array );
    }
    else
    {
        isSortedKernel<ValueType, false><<<dimGrid, dimBlock>>> ( resultRawPtr, n - 1, array );
    }

    cudaStreamSynchronize( 0 );

    LAMA_CHECK_CUDA_ERROR

    return thrust::reduce( resultPtr, resultPtr + n - 1, true, thrust::logical_and<bool>() );
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
                   "setGather<" << Scalar::getType<ValueType1>() << "," << Scalar::getType<ValueType2>() << ">( ..., n = " << n << ")" )

    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    gatherKernel<<<dimGrid, dimBlock>>>( out, in, indexes, n );

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
                   "setScatter<" << Scalar::getType<ValueType1>() << "," << Scalar::getType<ValueType2>() << ">( ..., n = " << n << ")" )

    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    scatter_kernel<<<dimGrid, dimBlock>>>( out, indexes, in, n );

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
                   "set<" << Scalar::getType<ValueType1>() << "," << Scalar::getType<ValueType2>() << ">( ..., n = " << n << ")" )

    LAMA_LOG_DEBUG( logger, "out = " << out << ", in = " << in )

    if ( n <= 0 )
    {
        return;
    }

    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    setKernel<<<dimGrid, dimBlock>>>( out, in, n );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void CUDAUtils::setScale( ValueType1 out[], 
                          const ValueType1 beta,
                          const ValueType2 in[], const IndexType n )
{
    LAMA_LOG_INFO( logger,
                   "set<" << Scalar::getType<ValueType1>() << "," << Scalar::getType<ValueType2>() << ">( ..., n = " << n << ")" )

    LAMA_LOG_DEBUG( logger, "out = " << out << ", in = " << in )

    if ( n <= 0 )
    {
        return;
    }

    if ( beta == static_cast<ValueType1>( 0 ) )
    {
        // in array might be undefined

        setVal( out, n, beta );
        return;
    }

    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = CUDASettings::getBlockSize( n );
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    setScaleKernel<<<dimGrid, dimBlock>>>( out, beta, in, n );

    LAMA_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );
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
                   "invert Vector components for vector of type " << Scalar::getType<ValueType>() << " and size n = " << n << "." )

    if ( n <= 0 )
    {
        return;
    }

    LAMA_CHECK_CUDA_ACCESS

    const int blockSize = 256;
    dim3 dimBlock( blockSize, 1, 1 );
    dim3 dimGrid = makeGrid( n, dimBlock.x );

    invertVectorComponents_kernel<<<dimGrid, dimBlock>>>( array, n );
    cudaStreamSynchronize( 0 );
    LAMA_CHECK_CUDA_ERROR
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDAUtils::setInterface( UtilsInterface& Utils )
{
    LAMA_LOG_INFO( logger, "set general utilty routines for CUDA in Interface" )

    LAMA_INTERFACE_REGISTER( Utils, validIndexes )

    LAMA_INTERFACE_REGISTER_T( Utils, sum, IndexType )

    LAMA_INTERFACE_REGISTER_T( Utils, setVal, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, setOrder, IndexType )

    LAMA_INTERFACE_REGISTER_T( Utils, getValue, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, maxval, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxVal, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxDiffVal, IndexType )

    LAMA_INTERFACE_REGISTER_T( Utils, isSorted, IndexType )

    LAMA_INTERFACE_REGISTER_TT( Utils, set, int, int )

    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, int, int )
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, int, int )

#define LAMA_UTILS2_REGISTER(z, J, TYPE )                                            \
    LAMA_INTERFACE_REGISTER_TT( Utils, setScale, TYPE, ARITHMETIC_TYPE##J )          \
    LAMA_INTERFACE_REGISTER_TT( Utils, set, TYPE, ARITHMETIC_TYPE##J )               \
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, TYPE, ARITHMETIC_TYPE##J )        \
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, TYPE, ARITHMETIC_TYPE##J )         \

#define LAMA_UTILS_REGISTER(z, I, _)                                                 \
    LAMA_INTERFACE_REGISTER_T( Utils, invert, ARITHMETIC_TYPE##I )                   \
    LAMA_INTERFACE_REGISTER_T( Utils, isSorted, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxDiffVal, ARITHMETIC_TYPE##I )            \
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxVal, ARITHMETIC_TYPE##I )                \
    LAMA_INTERFACE_REGISTER_T( Utils, maxval, ARITHMETIC_TYPE##I )                   \
    LAMA_INTERFACE_REGISTER_T( Utils, sum, ARITHMETIC_TYPE##I )                      \
    LAMA_INTERFACE_REGISTER_T( Utils, setVal, ARITHMETIC_TYPE##I )                   \
    LAMA_INTERFACE_REGISTER_T( Utils, getValue, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( Utils, scale, ARITHMETIC_TYPE##I )                    \
                                                                                     \
    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT,                                            \
                     LAMA_UTILS2_REGISTER,                                           \
                     ARITHMETIC_TYPE##I )                                            \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_UTILS_REGISTER, _ )

#undef LAMA_UTILS_REGISTER
#undef LAMA_UTILS2_REGISTER

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
