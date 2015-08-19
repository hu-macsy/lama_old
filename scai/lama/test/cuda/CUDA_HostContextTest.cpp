/**
 * @file CUDA_HostContextTest.cpp
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
 * @brief Contains the implementation of the class CUDA_HostContextTest.
 * @author: Alexander Büchel, Lauretta Schubert
 * @date 26.04.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DefaultHostContextManager.hpp>
#include <scai/lama/cuda/CUDAHostContextManager.hpp>
#include <scai/lama/cuda/CUDAHostContext.hpp>

#include <scai/lama/HostReadAccess.hpp>
#include <scai/lama/HostWriteAccess.hpp>
#include <scai/lama/LAMAArray.hpp>

#include <scai/lama/ContextAccess.hpp>

#include <test/cuda/CUDAContext.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDA_HostContextTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.CUDA_HostContextTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getContextTest )
{
    DefaultHostContextManager::setAsCurrent();
    ContextPtr defaultHostContext = ContextFactory::getContext( Context::Host );
    ContextPtr cudaHostContext = ContextFactory::getContext( Context::Host );
    // Two queries for the same context should deliver same pointer
    BOOST_CHECK( defaultHostContext.get() == cudaHostContext.get() );
    // Test will take the default CUDA device
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    CUDAHostContextManager::setAsCurrent( cudaContext );
    cudaHostContext = ContextFactory::getContext( Context::Host );
    SCAI_LOG_INFO( logger, "defaultHostContext = " << *defaultHostContext );
    SCAI_LOG_INFO( logger, "cudaHostContext = " << *cudaHostContext );
    BOOST_CHECK( defaultHostContext.get() != cudaHostContext.get() );
    // Note: the two shared pointers will be freed at the end of the subroutine
    //       so the CUDA context is freed now
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTest )
{
    // Problem: CUDA must be initialized before we can allocate memory via CUDAHostContext
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    CUDAHostContextManager::setAsCurrent( cudaContext );
    LAMAArray<IndexType> arrContext;
    {
        HostWriteAccess<IndexType> arr( arrContext );
        arr.resize( 5 );

        for ( IndexType i = 0; i < 5; i++ )
        {
            arr[i] = i;
        }
    }
    // Note: CUDA initialization done after call of lama_cudaHost_alloc
    {
        // transfer to CUDA context invalidates host instance
        WriteAccess<IndexType> arr( arrContext, cudaContext );
    }
    {
        HostWriteAccess<IndexType> arr( arrContext );

        for ( IndexType i = 0; i < 5; i++ )
        {
            BOOST_CHECK_EQUAL( i, arr[i] );
        }

        for ( IndexType i = 0; i < 5; i++ )
        {
            arr[i] = 0;
        }
    }
    arrContext.prefetch( cudaContext );
    {
        WriteAccess<IndexType> arr( arrContext, cudaContext );
    }
    arrContext.prefetch( ContextFactory::getContext( Context::Host ) );
    {
        HostReadAccess<IndexType> arr( arrContext );

        for ( IndexType i = 0; i < 5; i++ )
        {
            BOOST_CHECK_EQUAL( 0, arr[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

//BOOST_AUTO_TEST_CASE( prefetchTest )
//{
//    typedef double ValueType;
//
//    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
//
//    CUDAHostContextManager::setAsCurrent( cudaContext );
//
//    ContextPtr hostContext = ContextFactory::getContext( Context::Host );
//
//    const IndexType n = 100;
//
//    const ValueType value1 = 1.0;
//
//    LAMAArray<ValueType> vector1( n );
//    LAMAArray<ValueType> vector2( n );
//
//    {
//        HostWriteAccess<ValueType> v1(vector1);
//
//        for ( IndexType i = 0; i < n; ++i )
//        {
//            v1[i] = value1;
//        }
//    }
//
//    vector1.prefetch( cudaContext );
//
//    {
//        ReadAccess<ValueType> v1( vector1, cudaContext );
//
//        {
//            SCAI_CONTEXT_ACCESS( cudaContext );
//            double norm = lama_DNRM2_cuda( n, v1.get(), 1 );
//            double expNorm = sqrt( n * value1 );
//            BOOST_CHECK_EQUAL( expNorm, norm );
//        }
//
//        // vector2 = vector1   via copy
//
//        WriteAccess<ValueType> v2( vector2, cudaContext );
//
//        {
//            SCAI_CONTEXT_ACCESS( cudaContext );
//            lama_DCOPY_cuda( n, v1.get(), 1, v2.get(), 1 );
//        }
//    }
//
//    vector2.prefetch( hostContext );
//
//    {
//        HostReadAccess<ValueType> v2(vector2);
//
//        for ( IndexType i = 0; i < n; ++i )
//        {
//            BOOST_CHECK_EQUAL( value1, v2[i] );
//        }
//    }
//}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( multiPrefetchTest, ValueType, test_types )
{
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    CUDAHostContextManager::setAsCurrent( cudaContext );
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );
    const IndexType n = 128 * 1024;
    const ValueType value1 = 1.0;
    LAMAArray<ValueType> vector1( n, value1 );
    LAMAArray<ValueType> vector2( n, value1 );
    // Run two transfers Host -> CUDA in parallel
    vector1.prefetch( cudaContext );
    vector2.prefetch( cudaContext );
    {
        // invalidate vector data on host
        WriteAccess<ValueType> v1( vector1, cudaContext );
        WriteAccess<ValueType> v2( vector2, cudaContext );
    }
    // Run two transfers CUDA -> Host in parallel
    vector1.prefetch( hostContext );
    vector2.prefetch( hostContext );
    {
        WriteAccess<ValueType> v1( vector1, hostContext );
        WriteAccess<ValueType> v2( vector2, cudaContext );
    }
    // Run two transfers parallel for different directions
    vector1.prefetch( cudaContext );
    vector2.prefetch( hostContext );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
