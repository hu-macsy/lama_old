/**
 * @file CUDA_ContextTest.cpp
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
 * @brief Contains the implementation of the class CUDA_ContextTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 18.04.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <scai/common/bind.hpp>

#include <scai/lama/ContextFactory.hpp>
#include <scai/lama/WriteAccess.hpp>
#include <scai/lama/HArray.hpp>
#include <scai/lama/HostReadAccess.hpp>
#include <scai/lama/HostWriteAccess.hpp>
#include <scai/lama/ContextAccess.hpp>

#include <scai/lama/exception/LAMAmacros/assert.hpp>

#include <scai/tasking/Task.hpp>

#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/cuda/CUDAStreamSyncToken.hpp>
#include <scai/lama/cuda/CUDAContext.hpp>
#include <scai/lama/cuda/CUDAError.hpp>

#include <scai/lama/cuda/CUDABLAS1.hpp>
#include <scai/lama/test/cuda/CUDAContext.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using scai::tasking::Task;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDA_ContextTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.CUDA_ContextTest" );

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
    HArray<int> ctxArray; // default, not allocated at all
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    WriteAccess<int> array( ctxArray, cudaContext );
    array.resize( 10 );
    // Destructors will be called for array, ctxArray, cudaContext
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( releaseTest )
{
    HArray<IndexType> ctxArray; // default, not allocated at all
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
    HArray<IndexType> ctxArray; // default, not allocated at all
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
//BOOST_AUTO_TEST_CASE_TEMPLATE( prefetchTest, ValueType, test_types )
{
    typedef double ValueType;
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );
    const IndexType n = 100;
    const ValueType value1 = 1.4;
    HArray<ValueType> vector1( n );
    HArray<ValueType> vector2( n );
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
        SyncToken* token = cudaContext->getSyncToken();
        ValueType norm = nrm2( n, v1.get(), 1, token );
        delete token;
        ValueType expNormSqr = n * ( value1 * value1 );
        cudaContext->disable( __FILE__, __LINE__ );
        BOOST_CHECK_CLOSE( expNormSqr, norm * norm, 1 );
        // vector2 = vector1   via copy
        WriteAccess<ValueType> v2( vector2, cudaContext );
        cudaContext->enable( __FILE__, __LINE__ );
        token = cudaContext->getSyncToken();
        copy( n, v1.get(), 1, v2.get(), 1, token );
        delete token;
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
    HArray<float> vector( n, value );
    common::shared_ptr<WriteAccess<float> > cudaV( new WriteAccess<float>( vector, cudaContext ) );
    //CUDAContext* cuda = dynamic_cast<const CUDAContext*>( cudaContext.get() );
    const CUDAContext* cuda = ( CUDAContext* ) cudaContext.get();
    SCAI_ASSERT_ERROR( cuda, "dynamic cast for CUDAContext failed" );
    common::shared_ptr<CUDAStreamSyncToken> token( cuda->getComputeSyncToken() );
    {
        LAMA_INTERFACE_FN_t( scal, cudaContext, BLAS, BLAS1, float )
        SCAI_CONTEXT_ACCESS( cudaContext )
        // test: should be async!!!
        common::unique_ptr<SyncToken> token( cudaContext->getSyncToken() );
        scal( n, alpha, cudaV->get(), 1, token.get() );
        token->pushAccess( cudaV );
        token->wait();
        SCAI_CHECK_CUDA_ERROR;
    }
    cudaV.reset();  // give free the write access
    HostReadAccess<float> hostV( vector );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( syncTest )
{
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );
    const IndexType n = 100;
    const float value = 1.4;
    const float alpha = 0.5;
    HArray<float> vector( n, value );
    {
        LAMA_INTERFACE_FN_t( scal, cudaContext, BLAS, BLAS1, float );
        WriteAccess<float> cudaV( vector, cudaContext );
        SCAI_CONTEXT_ACCESS( cudaContext );
        common::unique_ptr<SyncToken> token( cudaContext->getSyncToken() );
        scal( n, alpha, cudaV.get(), 1, token.get() );
        // synchronize on token at end of this scope
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

static void callSSCAL( HArray<float>& vector, const float alpha, ContextPtr context )
{
    // get routine for context, must be available, otherwise Exception
    LAMA_INTERFACE_FN_t( scal, context, BLAS, BLAS1, float );
    WriteAccess<float> vectorAccess( vector, context );
    SCAI_CONTEXT_ACCESS( context );
    common::unique_ptr<SyncToken> token( context->getSyncToken() );
    scal( vector.size(), alpha, vectorAccess.get(), 1, token.get() );
    // syncronize on token is done here at end of scope where token will be destroyed
}

}

BOOST_AUTO_TEST_CASE( threadTest )
{
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );
    const IndexType n = 100;
    const float value = 1.4;
    const float alpha = 0.5;
    const HArray<float> vectorOrig( n, value );
    HArray<float> vector( n, value );
    HArray<float> vector2( n, value );
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
    common::function<void()> sscalCUDA = common::bind( callSSCAL, ref( vector ), alpha, cudaContext );
    common::function<void()> sscalHOST = common::bind( callSSCAL, ref( vector ), alpha, hostContext );
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

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

