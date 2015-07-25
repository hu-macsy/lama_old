/**
 * @file CUDA_LAMAArrayTest.cpp
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
 * @brief Contains the implementation of the class CUDA_LAMAArrayTest.
 * @author: Alexander Büchel, Lauretta Schubert
 * @date 18.04.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/scoped_array.hpp>

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAArray.hpp>
#include <lama/HostWriteAccess.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/ContextAccess.hpp>

#include <lama/exception/LAMAAssert.hpp>

#include <test/cuda/CUDAContext.hpp>

using namespace lama;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDA_LAMAArrayTest );

LAMA_LOG_DEF_LOGGER( logger, "Test.CUDA_LAMAArrayTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( baseTest )
{
    const IndexType n = 10;
    const IndexType val = 5;
    LAMAArray<IndexType> array1( n ); // will already be allocated on host
    LAMAArray<IndexType> array2; // no size
    ContextPtr cuda = ContextFactory::getContext( Context::CUDA );
    {
        UtilsInterface::Setter<IndexType>::setVal setVal = cuda->getInterface().Utils.setVal<IndexType>();
        BOOST_REQUIRE( setVal );
        {
            WriteAccess<IndexType> wArray1( array1, cuda );
            WriteOnlyAccess<IndexType> wArray2( array2, cuda, n );
            LAMA_CONTEXT_ACCESS( cuda );
            setVal( wArray1.get(), n, val );
            setVal( wArray2.get(), n, val + 1 );
        }
    }
    HostReadAccess<IndexType> rArray1( array1 );
    HostReadAccess<IndexType> rArray2( array2 );

    for ( IndexType i = 0; i < n; i++ )
    {
        BOOST_CHECK_EQUAL( rArray1[i] + 1, rArray2[i] );
    }
}

/* -------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( refTest1 )
{
    typedef float ValueType;
    ContextPtr cuda = ContextFactory::getContext( Context::CUDA );
    const IndexType n = 10;
    const ValueType value = 3.5;
    ValueType myData[n] =
    { 1, 2, 3, 4, 5, 5, 4, 3, 2, 1 };
    {
        // use the existing host data as input for LAMA array
        // LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType );
        UtilsInterface::Setter<ValueType>::setVal setVal = cuda->getInterface().Utils.setVal<ValueType>();
        BOOST_REQUIRE( setVal );
        LAMAArrayRef<ValueType> lamaArray( myData, n );
        {
            WriteAccess<ValueType> lamaArrayWAccess( lamaArray, cuda );
            LAMA_CONTEXT_ACCESS( cuda );
            setVal( lamaArrayWAccess.get(), n, value );
        }
        // This should become redundant by the destructor
        {
            HostWriteAccess<ValueType> lamaArrayWAccess( lamaArray );
        }
        // destructor of LAMAArrayRef will write back data to the valid location
    }

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value, myData[i] );
    }
}
//BOOST_AUTO_TEST_CASE( accessTestonCUDA )
//{
//    typedef float  ValueType;
//    const IndexType n = 10;
//    const ValueType value = 1.0;
//    const ValueType value2 = 2.0;
//    const ValueType value3 = 3.0;
//
//    ContextPtr cuda = lama_test::CUDAContext::getContext();
//
//    LAMAArray<ValueType> lamaArray( n, value );
//    {
//        ReadAccess<ValueType> cudaReadAccess( lamaArray, cuda );
//        boost::scoped_array<ValueType> tmpArray( new ValueType[n]);
//
////        TODO:
////        LAMA_CALL( lama_memcpyToHost_cuda( tmpArray.get(), cudaReadAccess.get(), n*sizeof( ValueType ) ),
////                   "lama_memcpyToHost_cuda failed" );
//
//        for ( IndexType i = 0; i < cudaReadAccess.size(); i++ )
//        {
//            tmpArray[i] = cudaReadAccess;
//        }
//
//        for ( IndexType i = 0; i < n; ++i )
//        {
//            BOOST_CHECK_EQUAL( value2, tmpArray[i] );
//        }
//        cudaReadAccess.release();
//        BOOST_CHECK_THROW( cudaReadAccess.get(), Exception );
//
//        WriteAccess<ValueType> cudaWriteAccess( lamaArray, cuda );
//        for ( IndexType i = 0; i < n; ++i )
//        {
//            tmpArray[i] = value3;
//        }
//
////        TODO:
////        LAMA_CALL( lama_memcpyToDevice_cuda(  cudaWriteAccess.get(), tmpArray.get(), n*sizeof( ValueType ) ),
////                   "lama_memcpyToDevice_cuda failed" );
//
//
//        cudaWriteAccess.release();
//        HostReadAccess<ValueType> lamaArrayRAccess(lamaArray);
//        for ( IndexType i = 0; i < n; ++i )
//        {
//            BOOST_CHECK_EQUAL( value3, lamaArrayRAccess[i] );
//        }
//    }
//}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

