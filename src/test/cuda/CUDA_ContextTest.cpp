/**
 * @file CUDA_ContextTest.cpp
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
 * @brief Contains the implementation of the class CUDA_ContextTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 18.04.2012
 * $
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/scoped_array.hpp>

#include <lama/ContextFactory.hpp>
#include <lama/WriteAccess.hpp>
#include <lama/LAMAArray.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>
#include <lama/ContextAccess.hpp>

#include <lama/exception/LAMAAssert.hpp>

#include <lama/task/Task.hpp>

#include <lama/LAMAInterface.hpp>
#include <lama/cuda/CUDAStreamSyncToken.hpp>
#include <lama/cuda/CUDAContext.hpp>
#include <lama/cuda/CUDAError.hpp>

#include <lama/cuda/CUDABLAS1.hpp>
#include <lama/openmp/OpenMPBLAS1.hpp>

#include <test/cuda/CUDAContext.hpp>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<double,float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDA_ContextTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.CUDA_ContextTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getContextTest )
{
    // Test will take the default CUDA device

    ContextPtr cudaContext1 = lama_test::CUDAContext::getContext();
    ContextPtr cudaContext2 = lama_test::CUDAContext::getContext();

    // Two queries for the same context should deliver same pointer

    BOOST_CHECK( cudaContext1.get() == cudaContext2.get() );

    cudaContext1 = ContextPtr();

    cudaContext1 = lama_test::CUDAContext::getContext();

    BOOST_CHECK( cudaContext1.get() == cudaContext2.get() );

    // Note: the two shared pointers will be freed at the end of the subroutine
    //       so the CUDA context is freed now
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTest )
{
    LAMAArray<int> ctxArray; // default, not allocated at all

    ContextPtr cudaContext = lama_test::CUDAContext::getContext();

    WriteAccess<int> array( ctxArray, cudaContext );

    array.resize( 10 );

    // Destructors will be called for array, ctxArray, cudaContext
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( releaseTest )
{
    LAMAArray<IndexType> ctxArray; // default, not allocated at all

    HostReadAccess<IndexType> readTestAccess( ctxArray );

    readTestAccess.release();

    HostWriteAccess<IndexType> writeAccess( ctxArray );

    writeAccess.resize( 10 );

    for ( IndexType i = 0; i < 10; i++ )
    {
        writeAccess[i] = 3;
    }

    writeAccess.release();

    //Should throw exception, because of already released writeAccess
    BOOST_CHECK_THROW( { writeAccess.resize( 20 ); }, Exception );

    // This is not checked:  writeAccess[0] = 5.0; -> crashes

    HostReadAccess<IndexType> readAccess( ctxArray );

    for ( IndexType i = 0; i < 5; i++ )
    {
        BOOST_CHECK_EQUAL( 3, readAccess[i] );
    }

    readAccess.release();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( resizeTest )
{
    LAMAArray<IndexType> ctxArray; // default, not allocated at all
    {
        HostWriteAccess<IndexType> writeAccess( ctxArray );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( 10 );

        for ( IndexType i = 0; i < 10; i++ )
        {
            writeAccess[i] = 3;
        }
    }

    ctxArray.purge();

    {
        HostWriteAccess<IndexType> writeAccess( ctxArray );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( 10 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( prefetchTest )
{
    typedef double ValueType;

    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );

    const IndexType n = 100;

    const ValueType value1 = 1.4;

    LAMAArray<ValueType> vector1( n );
    LAMAArray<ValueType> vector2( n );

    {
        HostWriteAccess<ValueType> v1( vector1 );

        for ( IndexType i = 0; i < n; ++i )
        {
            v1[i] = value1;
        }
    }

    vector1.prefetch( cudaContext );

    {
        LAMA_INTERFACE_FN_t(  nrm2, cudaContext, BLAS, BLAS1, ValueType );
        LAMA_INTERFACE_FN_t(  copy, cudaContext, BLAS, BLAS1, ValueType );

        ReadAccess<ValueType> v1( vector1, cudaContext );
        cudaContext->enable( __FILE__, __LINE__ );
        double norm = nrm2( n, v1.get(), 1, cudaContext->getSyncToken().get() );
        double expNormSqr = n * ( value1 * value1 );
        cudaContext->disable( __FILE__, __LINE__ );
        BOOST_CHECK_CLOSE( expNormSqr, norm * norm, 1 );

        // vector2 = vector1   via copy
        WriteAccess<ValueType> v2( vector2, cudaContext );
        cudaContext->enable( __FILE__, __LINE__ );
        copy( n, v1.get(), 1, v2.get(), 1, cudaContext->getSyncToken().get() );
        cudaContext->disable( __FILE__, __LINE__ );
    }

    vector2.prefetch( hostContext );
    vector1.prefetch( hostContext );

    {
        HostReadAccess<ValueType> v2( vector2 );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value1, v2[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( asyncTest )
{
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );

    const IndexType n = 100;

    const float value = 1.4;
    const float alpha = 0.5;

    LAMAArray<float> vector( n, value );

    CUDAStreamSyncTokenPtr token;

    std::auto_ptr<WriteAccess<float> > cudaV( new WriteAccess<float>( vector, cudaContext ) );

    //CUDAContext* cuda = dynamic_cast<const CUDAContext*>( cudaContext.get() );
    const CUDAContext* cuda = (CUDAContext*) cudaContext.get();

    LAMA_ASSERT_ERROR( cuda, "dynamic cast for CUDAContext failed" );

    token = cuda->getComputeSyncToken();

    {
        LAMA_INTERFACE_FN_t( scal, cudaContext, BLAS, BLAS1, float )

        LAMA_CONTEXT_ACCESS( cudaContext )

        // test: should be async!!!
        scal( n, alpha, cudaV->get(), 1, cudaContext->getSyncToken().get() );

        LAMA_CHECK_CUDA_ERROR
        ;
    }

    token->pushAccess( std::auto_ptr<BaseAccess>( cudaV ) );

    token->wait(); // frees/releases the pushed write access

    HostReadAccess<float> hostV( vector );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
    }
}

//TODO: To benchmarks?
//BOOST_AUTO_TEST_CASE( benchTest )
//{
//    ContextPtr cudaContext = ContextFactory::getContext( Context::CUDA );
//    ContextPtr hostContext = ContextFactory::getContext( Context::Host );
//
//    const int n = 8 * 1024;
//
//    const float value = 1.4;
//    const float alpha = 1.0;
//
//    LAMAArray<float> vector( n, value );
//
//    const int ITER = 10000;
//
//    // Avoid measuring overhead of first call
//
//    {
//        WriteAccess<float> cudaV( vector, cudaContext ) ;
//        LAMA_CONTEXT_ACCESS( cudaContext );
//        lama_SSCAL_cuda(n, alpha, cudaV.get(), 1 );
//        LAMA_CHECK_ERROR( "lama_SSCAL_cuda" );
//    }
//
//    double time_pure = 0.0;
//    double time_access = 0.0;
//    double time_full = 0.0;
//
//    {
//        double start = timestamp();
//        WriteAccess<float> cudaV( vector, cudaContext ) ;
//        LAMA_CONTEXT_ACCESS( cudaContext );
//        for ( int i = 0; i < ITER; i++)
//        {
//            lama_SSCAL_cuda(n, alpha, cudaV.get(), 1 );
//            LAMA_CHECK_ERROR( "lama_SSCAL_cuda" );
//        }
//        double stop = timestamp();
//        time_pure = (stop - start) / double(ITER);
//    }
//
//    {
//        double start = timestamp();
//        LAMA_CONTEXT_ACCESS( cudaContext );
//        for ( int i = 0; i < ITER; i++)
//        {
//            WriteAccess<float> cudaV( vector, cudaContext ) ;
//            lama_SSCAL_cuda(n, alpha, cudaV.get(), 1 );
//            LAMA_CHECK_ERROR( "lama_SSCAL_cuda" );
//        }
//        double stop = timestamp();
//        time_access = (stop - start) / double(ITER);
//    }
//
//    {
//        double start = timestamp();
//        for ( int i = 0; i < ITER; i++)
//        {
//            WriteAccess<float> cudaV( vector, cudaContext ) ;
//            LAMA_CONTEXT_ACCESS( cudaContext );
//            lama_SSCAL_cuda(n, alpha, cudaV.get(), 1 );
//            LAMA_CHECK_ERROR( "lama_SSCAL_cuda" );
//        }
//        double stop = timestamp();
//        time_full = (stop - start) / double(ITER);
//    }
//
//    cout << "Execution times in milliseconds for 1 call (average over "
//              << ITER << " iterations)" << endl;
//
//    cout << "   pure    = " << ( time_pure * 1000.0 ) << " ms" << endl;
//    cout << "   access  = " << ( time_access * 1000.0 ) << " ms" << endl;
//    cout << "   full    = " << ( time_full * 1000.0 ) << " ms" << endl;
//
//    double overhead1 = ( time_access - time_pure ) / time_pure * 100;
//    double overhead2 = ( time_full - time_pure ) / time_pure * 100;
//
//    cout << "Overhead access = " << overhead1 << " %" << endl;
//    cout << "Overhead full   = " << overhead2 << " %" << endl;
//}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( syncTest )
{
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );

    const IndexType n = 100;

    const float value = 1.4;
    const float alpha = 0.5;

    LAMAArray<float> vector( n, value );

    {
        LAMA_INTERFACE_FN_t( scal, cudaContext, BLAS, BLAS1, float );

        WriteAccess<float> cudaV( vector, cudaContext );
        LAMA_CONTEXT_ACCESS( cudaContext );

        scal( n, alpha, cudaV.get(), 1, cudaContext->getSyncToken().get() );
    }

    HostReadAccess<float> hostV( vector );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
    }
}

/* --------------------------------------------------------------------- */

namespace
{

static void callSSCAL( LAMAArray<float>& vector, const float alpha, ContextPtr context )
{
    // get routine for context, must be available, otherwise Exception

    LAMA_INTERFACE_FN_t( scal, context, BLAS, BLAS1, float );

    WriteAccess<float> vectorAccess( vector, context );

    LAMA_CONTEXT_ACCESS( context );

    scal( vector.size(), alpha, vectorAccess.get(), 1, context->getSyncToken().get() );
}

}

BOOST_AUTO_TEST_CASE( threadTest )
{
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );

    const IndexType n = 100;

    const float value = 1.4;
    const float alpha = 0.5;

    const LAMAArray<float> vectorOrig( n, value );
    LAMAArray<float> vector( n, value );
    LAMAArray<float> vector2( n, value );

    //CUDA Synchronous
    {
        callSSCAL( vector, alpha, cudaContext );
        HostReadAccess<float> hostV( vector );
        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
        }
    }

    vector = vectorOrig;

    //Host Synchronous
    {
        callSSCAL( vector, alpha, hostContext );
        HostReadAccess<float> hostV( vector );
        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
        }
    }

    vector = vectorOrig;
    boost::function<void()> sscalCUDA = boost::bind( callSSCAL, ref( vector ), alpha, cudaContext );
    boost::function<void()> sscalHOST = boost::bind( callSSCAL, ref( vector ), alpha, hostContext );

    //CUDA Asynchronous Simple
    {
        Task asyncCall( sscalCUDA );
        asyncCall.synchronize();
        HostReadAccess<float> hostV( vector );
        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
        }
    }

    vector = vectorOrig;

    //Host Asynchronous Simple
    {
        Task asyncCall( sscalHOST );
        asyncCall.synchronize();
        HostReadAccess<float> hostV( vector );
        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
        }
    }

    vector2 = vectorOrig;
    vector = vectorOrig;

    //CUDA Asynchronous
    {
        Task asyncCall( sscalCUDA );
        callSSCAL( vector2, alpha, cudaContext );
        asyncCall.synchronize();

        HostReadAccess<float> hostV( vector );
        HostReadAccess<float> hostV2( vector2 );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
            BOOST_CHECK_EQUAL( value * alpha, hostV2[i] );
        }
    }

    vector2 = vectorOrig;
    vector = vectorOrig;

    //Host Asynchronous
    {
        Task asyncCall( sscalHOST );
        callSSCAL( vector2, alpha, hostContext );
        asyncCall.synchronize();

        HostReadAccess<float> hostV( vector );
        HostReadAccess<float> hostV2( vector2 );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
            BOOST_CHECK_EQUAL( value * alpha, hostV2[i] );
        }
    }

    vector2 = vectorOrig;
    vector = vectorOrig;

    //CUDA Synchronous
    {
        callSSCAL( vector, alpha, cudaContext );
        HostReadAccess<float> hostV( vector );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
        }
    }
}
/* ------------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();

