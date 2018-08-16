/**
 * @file hmemo/test/cuda/HArrayTest.cu
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
 * @brief Basic tests for LAMA arrays with context/memory at CUDA devices
 * @author Thomas Brandes
 * @date 08.07.2015
 */

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
    ContextPtr host = Context::getContextPtr( ContextType::Host );
    ContextPtr cuda = Context::getContextPtr( ContextType::CUDA );
    HArray<float> array( cuda );
    initArray<float>( array, 128, 1.0 );
    //std::cout << "array : " << array << std::endl;
    BOOST_ASSERT( array.isValid( host ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( PrefetchTest )
{
    ContextPtr host = Context::getContextPtr( ContextType::Host );
    ContextPtr cuda = Context::getContextPtr( ContextType::CUDA );
    IndexType N = 128;
    HArray<float> array( cuda );
    initArray<float>( array, N, 1.0 );
    //std::cout << "array : " << array << std::endl;
    BOOST_ASSERT( array.isValid( host ) );
    BOOST_ASSERT( !array.isValid( cuda ) );
    array.prefetch( cuda );
    // Note: array is already valid at cuda context even if not finished
    //std::cout << "array : " << array << std::endl;
    BOOST_ASSERT( array.isValid( cuda ) );
    // wait: would also be done before any new access
    float expected = N * ( N - 1 ) / 2 + N;
    float result = cudaSum( array, cuda );
    //std::cout << "sum (on CUDA device) is " << result << ", expected = " << expected << std::endl;
    BOOST_CHECK_EQUAL( result, expected );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CopyTest )
{
    ContextPtr host = Context::getContextPtr( ContextType::Host );
    ContextPtr cuda = Context::getContextPtr( ContextType::CUDA );
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
    //std::cout << "HArray, copy: sum (on CUDA device) is " << result << ", expected = " << expected << std::endl;
    BOOST_CHECK_EQUAL( result, expected );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
