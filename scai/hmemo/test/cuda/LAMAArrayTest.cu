/**
 * @file cudamem/test/HArrayTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Basic tests for LAMA arrays with context/memory at CUDA devices
 * @author: Thomas Brandes
 * @date 08.07.2015
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/hmemo.hpp>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#include <scai/common/cuda/CUDAError.hpp>

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HArrayTest )

using namespace scai;
using namespace scai::hmemo;
using namespace scai::common;

/* --------------------------------------------------------------------- */

template <typename T>
void initArray( HArray<T>& array, int N, T val )
{
    WriteOnlyAccess<T> write( array, N );
    T* data = write.get();

    for ( int i = 0; i < N; ++i )
    {
        data[i] = val + static_cast<T>( i );
    }
}

/* --------------------------------------------------------------------- */

template<typename ValueType>
ValueType sum( const ValueType array[], const IndexType n )
{
    thrust::device_ptr<ValueType> data( const_cast<ValueType*>( array ) );

    ValueType zero = static_cast<ValueType>( 0 );

    ValueType result = thrust::reduce( data, data + n, zero, thrust::plus<ValueType>() );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( 0 ), "cudaStreamSynchronize( 0 )" );

    return result;
}

template <typename ValueType>
ValueType cudaSum( HArray<ValueType>& array, ContextPtr cuda )
{
    IndexType n = array.size();
    WriteAccess<ValueType> write( array, cuda );
    return sum ( write.get(), n );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    ContextPtr host = Context::getContextPtr( context::Host );
    ContextPtr cuda = Context::getContextPtr( context::CUDA );

    HArray<float> array( cuda );

    initArray<float>( array, 128, 1.0 );

    std::cout << "array : " << array << std::endl;

    BOOST_ASSERT( array.isValid( host ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( PrefetchTest )
{
    ContextPtr host = Context::getContextPtr( context::Host );
    ContextPtr cuda = Context::getContextPtr( context::CUDA );

    IndexType N = 128;

    HArray<float> array( cuda );

    initArray<float>( array, N, 1.0 );

    std::cout << "array : " << array << std::endl;

    BOOST_ASSERT( array.isValid( host ) );
    BOOST_ASSERT( !array.isValid( cuda ) );

    array.prefetch( cuda );

    // Note: array is already valid at cuda context even if not finished

    std::cout << "array : " << array << std::endl;

    BOOST_ASSERT( array.isValid( cuda ) );

    // wait: would also be done before any new access

    float expected = N * ( N - 1 ) / 2 + N;
    float result = cudaSum( array, cuda );

    std::cout << "sum (on CUDA device) is " << result << ", expected = " << expected << std::endl;
    BOOST_CHECK_EQUAL( result, expected );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CopyTest )
{
    ContextPtr host = Context::getContextPtr( context::Host );
    ContextPtr cuda = Context::getContextPtr( context::CUDA );

    IndexType N = 128;

    HArray<float> array( cuda );

    initArray<float>( array, N, 1.0 );
    // array is valid @ host
    array.prefetch( cuda );
    // now array is valid @ host, cuda

    // Copy contructor, copies all valid data at their locations
    HArray<float> array1( array );

    BOOST_ASSERT( array1.isValid( host ) );
    BOOST_ASSERT( array1.isValid( cuda ) );

    float expected = cudaSum( array, cuda );
    float result = cudaSum( array1, cuda );

    std::cout << "HArray, copy: sum (on CUDA device) is " << result 
              << ", expected = " << expected << std::endl;

    BOOST_CHECK_EQUAL( result, expected );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
