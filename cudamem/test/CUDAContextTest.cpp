/**
 * @file CUDAContextTest.cpp
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
 * @brief Extensive testing of CUDAContext.
 * @author: Thomas Brandes
 * @date 18.07.2015
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <memory/memory.hpp>

#include <tasking/Task.hpp>

#include <cudamem/CUDAStreamSyncToken.hpp>
#include <cudamem/CUDAContext.hpp>
#include <cudamem/CUDAError.hpp>

#include <common/shared_ptr.hpp>
#include <common/function.hpp>
#include <common/bind.hpp>

using tasking::Task;
using tasking::SyncToken;
using tasking::CUDAStreamSyncToken;

using namespace memory;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

extern cublasHandle_t CUDAContext_cublasHandle;

static void scal( int n, float alpha, float* x_d, int inc_x, SyncToken* syncToken )
{
    LAMA_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        COMMON_ASSERT( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    std::cout << "scal( n = " << n << ", alpha = " << alpha << ", x[], inc_x = " << inc_x << std::endl;

    LAMA_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "scal set stream" );
    LAMA_CUBLAS_CALL( cublasSscal( CUDAContext_cublasHandle, n, &alpha, x_d, inc_x ), cublasSscal );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDAContextTest );

LAMA_LOG_DEF_LOGGER( logger, "Test.CUDAContextTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getContextTest )
{
    // Test will take the default CUDA device
    ContextPtr cudaContext1 = Context::getContext( context::CUDA );
    ContextPtr cudaContext2 = Context::getContext( context::CUDA );
    // Two queries for the same context should deliver same pointer
    BOOST_CHECK( cudaContext1.get() == cudaContext2.get() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTest )
{
    LAMAArray<int> ctxArray; // default, not allocated at all
    ContextPtr cudaContext = Context::getContext( context::CUDA );
    {
        WriteAccess<int> array( ctxArray, cudaContext );
        array.resize( 10 );
        // Destructors will be called for array, ctxArray, cudaContext
    }
    BOOST_CHECK_EQUAL( 10, ctxArray.size() );
    BOOST_CHECK( ctxArray.capacity( cudaContext ) >= 10 );
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
    BOOST_CHECK_THROW( { writeAccess.resize( 20 ); }, common::Exception );
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

BOOST_AUTO_TEST_CASE( asyncTest )
{
    ContextPtr cudaContext = Context::getContext( context::CUDA );
    const IndexType n = 100;
    const float value = 1.4;
    const float alpha = 0.5;
    LAMAArray<float> vector( n, value );
    common::shared_ptr<WriteAccess<float> > cudaV( new WriteAccess<float>( vector, cudaContext ) );

    common::shared_ptr<SyncToken> token( cudaContext->getSyncToken() );

    {
        LAMA_CONTEXT_ACCESS( cudaContext )
        scal( n, alpha, cudaV->get(), 1, token.get() );
        token->pushToken( cudaV );
    }

    cudaV.reset();  // give free the write access, but ownership also in token

    token->wait();  

    HostReadAccess<float> hostV( vector );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( syncTest )
{
    ContextPtr cudaContext = Context::getContext( context::CUDA );

    const IndexType n = 100;
    const float value = 1.4;
    const float alpha = 0.5;
    LAMAArray<float> vector( n, value );
    {
        WriteAccess<float> cudaV( vector, cudaContext );
        LAMA_CONTEXT_ACCESS( cudaContext );
        std::auto_ptr<SyncToken> token( cudaContext->getSyncToken() );
        scal( n, alpha, cudaV.get(), 1, token.get() );
        // synchronize on token at end of this scope
    }
    HostReadAccess<float> hostV( vector );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV[i] );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

