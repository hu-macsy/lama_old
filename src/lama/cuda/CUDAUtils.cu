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

///* sparse dot */
//
////TODO: Implement the Texture Cache thing here like it is done
////      for level2 and lapack stuff.
///** The texture holding the data of the vector. */
//texture<float, 1> tex_Y_Data_ref_d;
//
///**
// * @brief Returns the value of the vector at position i.
// * @param[in] i Position.
// * @return The value of the vector at position i.
// */
//template<typename T>
//__inline__ __device__
//T fetch_Y_Data( const int i )
//{
//    return tex1Dfetch( tex_Y_Data_ref_d, i );
//}
//
///**
// * @brief Adds the stride of the given array to the value in the array at
// *        position of the threadId, if the threadId is smaller than the stride.
// * @param[in,out] sum   The array holding the values to be added.
// * @param[in] threadId  The threadId.
// * @param[in] stride    The stride to add.
// */
//template<typename T>
//__inline__ __device__
//void stridedAdd( T* sum, const unsigned int threadId, const unsigned int stride )
//{
//    if ( threadId < stride )
//    {
//        sum[threadId] += sum[threadId + stride];
//    }
//}
//
///**
// * @brief multiplies SparseVector x and DenseVector y, each block reduces its
// * product to one value.
// *
// * multiplies SparseVector x and DenseVector y, each block reduces its
// * product to one value. So the size of the global array has to be just the size
// * of the number of blocks.
// *
// * @param[in]  nz length of SparseVector.
// * @param[in]  x_d The array holding the data of the sparseVector.
// * @param[in]  indx_d The array holding the rows of the sparseVector.
// * @param[out] product_d global array with size == gridDim.x*gridDim.y
// *
// * TODO: optimization.
// */
//template<typename ValueType, unsigned int blocksize>
//__global__
//void sdoti_kernel(
//    const int nz,
//    const ValueType* const x_d,
//    const int* const indx_d,
//    ValueType* product_d )
//{
//    extern __shared__ float localSharedProduct[ ];
//    localSharedProduct[threadIdx.x] = 0.0;
//    {
//        /* Block to make thisColumn local */
//        unsigned const int thisColumn = threadId( gridDim, blockIdx, blockDim, threadIdx );
//
//        if ( thisColumn < nz )
//        {
//            localSharedProduct[threadIdx.x] =
//                x_d[thisColumn] * fetch_Y_Data<float>( indx_d[thisColumn] );
//        }
//    } /* Block to make thisColumn local */
//    __syncthreads( );
//
//    if ( blocksize >= 512 )
//    {
//        stridedAdd( localSharedProduct, threadIdx.x, 256 );
//        __syncthreads( );
//    }
//
//    if ( blocksize >= 256 )
//    {
//        stridedAdd( localSharedProduct, threadIdx.x, 128 );
//        __syncthreads( );
//    }
//
//    if ( blocksize >= 128 )
//    {
//        stridedAdd( localSharedProduct, threadIdx.x, 64 );
//        __syncthreads( );
//    }
//
//    /* unrolling the last six steps of the loop. This is possible, because a  */
//    /* warp is always doing the same at a time and 32 threads are a warp.     */
//    if ( threadIdx.x < 32 )
//    {
//        if ( blocksize > 32 )
//        {
//            localSharedProduct[threadIdx.x] += localSharedProduct[threadIdx.x + 32];
//        }
//
//        localSharedProduct[threadIdx.x] += localSharedProduct[threadIdx.x + 16];
//        localSharedProduct[threadIdx.x] += localSharedProduct[threadIdx.x + 8];
//        localSharedProduct[threadIdx.x] += localSharedProduct[threadIdx.x + 4];
//        localSharedProduct[threadIdx.x] += localSharedProduct[threadIdx.x + 2];
//        localSharedProduct[threadIdx.x] += localSharedProduct[threadIdx.x + 1];
//    }
//
//    /* write sum of this block to global array */
//    if ( threadIdx.x == 0 )
//    {
//        product_d[blockId( gridDim, blockIdx )] = localSharedProduct[0];
//    }
//} /* cspblas_sdoti_kernel */
//
///**
// * @brief Reduces values of global array to 1/threadsPerBlock
// *
// * Each block sums up threadsPerBlock values of the global array, so each time,
// * the kernel is invoked, it reduces the number of values to 1/threadsPerBlock,
// * but at least 1.
// *
// * @param[in,out] product_d global array, containing values to reduce
// * @param[in]     sizeOfProduct size of product_d
// *
// * TODO: optimization.
// */
//template<typename ValueType, unsigned int blocksize>
//__global__
//void sum_reduction_kernel(
//    ValueType* product_d,
//    int sizeOfProduct )
//{
//    extern __shared__ float sum[ ];
//    {
//        /* Block to make thisColumn and interval local */
//        unsigned const int thisColumn =
//            threadId( gridDim, blockIdx, blockDim, threadIdx );
//        unsigned const int interval = halve( sizeOfProduct );
//
//        if ( thisColumn + interval < sizeOfProduct )
//        {
//            sum[threadIdx.x] = product_d[thisColumn]
//                               + product_d[thisColumn + interval];
//        }
//        /* thisColumn is the middle of the vector                             */
//        /* example:                                                           */
//        /* Input ( 0 1 2 3 4 5 6 7 8 ) => sizeOfProduct == 9                  */
//        /* => interval == 5                                                   */
//        /* we need to handle the case thisColumn == 4, because 4 + 5 >= 9     */
//        /* and the element at position 4 will not be handled by threadIdx.x==0*/
//        /* as it would be with an vector of even length                       */
//        else if ( thisColumn + interval == sizeOfProduct )
//        {
//            /* if the vector is odd, the midst value must be saved. */
//            sum[threadIdx.x] = ( ( sizeOfProduct & 1 ) ? product_d[thisColumn] : 0.0 );
//        }
//        else
//        {
//            sum[threadIdx.x] = 0.0;
//        }
//
//        __syncthreads( );
//    } /* Block to make thisColumn and interval local */
//
//    if ( blocksize >= 512 )
//    {
//        stridedAdd( sum, threadIdx.x, 256 );
//        __syncthreads( );
//    }
//
//    if ( blocksize >= 256 )
//    {
//        stridedAdd( sum, threadIdx.x, 128 );
//        __syncthreads( );
//    }
//
//    if ( blocksize >= 128 )
//    {
//        stridedAdd( sum, threadIdx.x, 64 );
//        __syncthreads( );
//    }
//
//    /* unrolling the last six steps of the loop. This is possible, because a  */
//    /* warp is always doing the same at a time and 32 threads are a warp.     */
//    if ( threadIdx.x < 32 )
//    {
//        if ( blocksize > 32 )
//        {
//            sum[threadIdx.x] += sum[threadIdx.x + 32];
//        }
//
//        sum[threadIdx.x] += sum[threadIdx.x + 16];
//        sum[threadIdx.x] += sum[threadIdx.x + 8];
//        sum[threadIdx.x] += sum[threadIdx.x + 4];
//        sum[threadIdx.x] += sum[threadIdx.x + 2];
//        sum[threadIdx.x] += sum[threadIdx.x + 1];
//    }
//
//    /* write sum of this block to global array. */
//    if ( threadIdx.x == 0 )
//    {
//        product_d[blockId( gridDim, blockIdx )] = sum[0];
//    }
//} /* sum_reduction_kernel */
//
//template<typename ValueType>
//ValueType CUDAUtils::sdoti(
//    const IndexType nz,
//    const ValueType* const x_d,
//    const IndexType* const indx_d,
//    const ValueType* const y_d,
//    ContextPtr context )
//{
//    LAMA_CHECK_CUDA_ACCESS;
//
//    if ( nz <= 0 )
//    {
//        return 0.0;
//    }
//
//    /* The getenv may be used to check  */
//    /* the performance with different   */
//    /* block sizes. The function then   */
//    /* can be called easily from an     */
//    /* external script.                 */
//    const IndexType SDOTI_numberOfThreads = 256; //atoi( getenv("SDOTI_BLOCK_SIZE") );
//    const IndexType RED_numberOfThreads   = 128; //atoi( getenv("RED_BLOCK_SIZE") );
//    /* define size of grid and blocks */
//    IndexType numberOfThreads = SDOTI_numberOfThreads;
//    dim3 numberOfBlocks = makeGrid( nz, numberOfThreads );
//    /* allocate global array for kernel */
//    IndexType sizeOfProduct = numberOfBlocks.x * numberOfBlocks.y;
//    ValueType* product_d;
//    IndexType error;
//    lama_Status lama_error;
//    //TODO
////    lama_Status lama_error = lama_deviceAlloc_cuda(
////        ( void** )&product_d,
////        sizeOfProduct*sizeof( float ) );
////
////    /* check for failures */
////    if( LAMA_STATUS_SUCCESS != lama_error )
////    {
////        lama_setLastError( lama_error );
////        return 0.0;
////    }
//    /* bind textures. */
//    LAMA_CUDA_SETLASTERROR( cudaBindTexture( NULL, tex_Y_Data_ref_d, y_d ),     \
//                            LAMA_STATUS_CUDA_BINDTEX_FAILED );
//
//    /* invoke kernel for mutliplication. */
//    // TODO switch is not used, because number of threads is fixed.
//
//    switch ( numberOfThreads )
//    {
//        case 512:
//            sdoti_kernel<512>
//            <<< numberOfBlocks, numberOfThreads, numberOfThreads * sizeof( ValueType )>>>(
//                nz,
//                x_d,
//                indx_d,
//                product_d );
//            break;
//        case 256:
//            sdoti_kernel<256>
//            <<< numberOfBlocks, numberOfThreads, numberOfThreads* sizeof( ValueType )>>>(
//                nz,
//                x_d,
//                indx_d,
//                product_d );
//            break;
//        case 128:
//            sdoti_kernel<128>
//            <<< numberOfBlocks, numberOfThreads, numberOfThreads* sizeof( ValueType )>>>(
//                nz,
//                x_d,
//                indx_d,
//                product_d );
//            break;
//        case  64:
//            sdoti_kernel< 64>
//            <<< numberOfBlocks, numberOfThreads, numberOfThreads* sizeof( ValueType )>>>(
//                nz,
//                x_d,
//                indx_d,
//                product_d );
//            break;
//        case  32:
//            sdoti_kernel< 32>
//            <<< numberOfBlocks, numberOfThreads, numberOfThreads* sizeof( ValueType )>>>(
//                nz,
//                x_d,
//                indx_d,
//                product_d );
//            break;
//    }
//
//    /* check for failures. */
//    error = cudaGetLastError( );
//    LAMA_CUDA_SETLASTERROR( error, LAMA_STATUS_SDOTI_CUDAKERNEL_FAILED );
//
//    if ( cudaSuccess != error )
//    {
//        //TODO
////        lama_deviceFree_cuda( product_d );
//        return 0.0;
//    }
//
//    error = cudaStreamSynchronize( 0 );
//    LAMA_CUDA_SETLASTERROR( error, LAMA_STATUS_CUDA_THREADSYNCHRONIZE_FAILED );
//
//    if ( cudaSuccess != error )
//    {
//        //TODO
////        lama_deviceFree_cuda( product_d );
//        return 0.0;
//    }
//
//    numberOfThreads = RED_numberOfThreads;
//
//    /* invokes kernel for reduction of global array in a loop, until one */
//    /* value is left.                                                    */
//    do
//    {
//        numberOfBlocks = makeGrid( sizeOfProduct, numberOfThreads );
//
//        /* invoke kernel for sum reduction. */
//        // TODO switch is not used, because number of threads is fixed.
//
//        switch ( numberOfThreads )
//        {
//            case 512:
//                sum_reduction_kernel<512>
//                <<< numberOfBlocks, numberOfThreads, numberOfThreads * sizeof( ValueType )>>>(
//                    product_d,
//                    sizeOfProduct );
//                break;
//            case 256:
//                sum_reduction_kernel<256>
//                <<< numberOfBlocks, numberOfThreads, numberOfThreads* sizeof( ValueType )>>>(
//                    product_d,
//                    sizeOfProduct );
//                break;
//            case 128:
//                sum_reduction_kernel<128>
//                <<< numberOfBlocks, numberOfThreads, numberOfThreads* sizeof( ValueType )>>>(
//                    product_d,
//                    sizeOfProduct );
//                break;
//            case  64:
//                sum_reduction_kernel< 64>
//                <<< numberOfBlocks, numberOfThreads, numberOfThreads* sizeof( ValueType )>>>(
//                    product_d,
//                    sizeOfProduct );
//                break;
//            case  32:
//                sum_reduction_kernel< 32>
//                <<< numberOfBlocks, numberOfThreads, numberOfThreads* sizeof( ValueType )>>>(
//                    product_d,
//                    sizeOfProduct );
//                break;
//        }
//
//        /* check for failures. */
//        error = cudaGetLastError( );
//        LAMA_CUDA_SETLASTERROR( error, LAMA_STATUS_SDOTI_SUMREDUCTION_CUDAKERNEL_FAILED );
//
//        if ( cudaSuccess != error )
//        {
//            //TODO
////            lama_deviceFree_cuda( product_d );
//            return 0.0;
//        }
//
//        /* the size of the relevant values in the array is the number of just */
//        /* invoked blocks, because each block has written out exact one value.*/
//        sizeOfProduct = numberOfBlocks.x * numberOfBlocks.y;
//    }
//    while ( numberOfBlocks.x > 1 );
//
//    error = cudaStreamSynchronize( 0 );
//    LAMA_CUDA_SETLASTERROR( error, LAMA_STATUS_CUDA_THREADSYNCHRONIZE_FAILED );
//
//    if ( cudaSuccess != error )
//    {
//        //TODO
////        lama_deviceFree_cuda( product_d );
//        return 0.0;
//    }
//
//    /* copy final scalar value to host. */
//    ValueType scalar = 0.0;
//    //TODO
////    lama_error = lama_memcpyToHost_cuda(
////            ( void* )&scalar,
////            product_d,
////            sizeof( float ) );
//
//    /* check for failures */
//    if ( LAMA_STATUS_SUCCESS != lama_error )
//    {
//        lama_setLastError( lama_error );
//        //TODO
////        lama_deviceFree_cuda( product_d );
//        return 0.0;
//    }
//
//    //TODO
////    lama_error = lama_deviceFree_cuda( ( void* )product_d );
//
//    /* check for failures */
//    if ( LAMA_STATUS_SUCCESS != lama_error )
//    {
//        lama_setLastError( lama_error );
//        return 0.0;
//    }
//
//    /* unbind textures. */
//    LAMA_CUDA_SETLASTERROR( cudaUnbindTexture( tex_Y_Data_ref_d ),            \
//                            LAMA_STATUS_CUDA_UNBINDTEX_FAILED );
//    lama_setLastError( LAMA_STATUS_SUCCESS );
//    return scalar;
//} /* lama_CSPBLAS_sdoti_launcher */


/* --------------------------------------------------------------------------- */

template<typename T>
struct InvalidIndex
{
    const T size;  //!< size of array for which index is checked

    InvalidIndex( T _size ) : size( _size ) {}

    __host__ __device__
    int operator()( T y )
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
                                          0,
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

} // namespace lama
