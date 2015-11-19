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

#include <scai/hmemo.hpp>

#include <scai/tasking/Task.hpp>

#include <scai/hmemo/cuda/CUDAStreamSyncToken.hpp>
#include <scai/hmemo/cuda/CUDAContext.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/shared_ptr.hpp>
#include <scai/common/function.hpp>
#include <scai/common/bind.hpp>

using namespace scai;

using tasking::Task;
using tasking::SyncToken;
using tasking::CUDAStreamSyncToken;

using namespace hmemo;
namespace context = common::context;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */
namespace scai
{
	extern cublasHandle_t CUDAContext_cublasHandle;
} /* end namespace scai */

static void scal( int n, float alpha, float* x_d, int inc_x, SyncToken* syncToken )
{
    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = NULL;

    if ( syncToken )
    {
        CUDAStreamSyncToken* cudaStreamSyncToken = dynamic_cast<CUDAStreamSyncToken*>( syncToken );
        SCAI_ASSERT( cudaStreamSyncToken, "no cuda stream sync token provided" )
        stream = cudaStreamSyncToken->getCUDAStream();
    }

    std::cout << "scal( n = " << n << ", alpha = " << alpha << ", x[], inc_x = " << inc_x << std::endl;

    SCAI_CUBLAS_CALL( cublasSetStream( CUDAContext_cublasHandle, stream ), "scal set stream" );
    SCAI_CUBLAS_CALL( cublasSscal( CUDAContext_cublasHandle, n, &alpha, x_d, inc_x ), "cublasSscal" );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CUDAContextTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.CUDAContextTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getContextTest )
{
    // Test will take the default CUDA device
    ContextPtr cudaContext1 = Context::getContextPtr( context::CUDA );
    ContextPtr cudaContext2 = Context::getContextPtr( context::CUDA );
    // Two queries for the same context should deliver same pointer
    BOOST_CHECK( cudaContext1.get() == cudaContext2.get() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTest )
{
    HArray<int> ctxArray; // default, not allocated at all
    ContextPtr cudaContext = Context::getContextPtr( context::CUDA );
    {
        WriteAccess<int> array( ctxArray, cudaContext );
        array.resize( 10 );
        // Destructors will be called for array, ctxArray, cudaContext
    }
    BOOST_CHECK_EQUAL( static_cast<IndexType>( 10 ), ctxArray.size() );
    BOOST_CHECK( ctxArray.capacity( cudaContext ) >= 10 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( releaseTest )
{
    ContextPtr contextPtr = Context::getContextPtr( context::Host );

    HArray<IndexType> ctxArray; // default, not allocated at all
    ReadAccess<IndexType> readTestAccess( ctxArray, contextPtr );
    readTestAccess.release();
    WriteAccess<IndexType> writeAccess( ctxArray, contextPtr );
    writeAccess.resize( 10 );

    for ( IndexType i = 0; i < 10; i++ )
    {
        writeAccess.get()[i] = 3;
    }

    writeAccess.release();
    //Should throw exception, because of already released writeAccess
    BOOST_CHECK_THROW( { writeAccess.resize( 20 ); }, common::Exception );
    // This is not checked:  writeAccess[0] = 5.0; -> crashes
    ReadAccess<IndexType> readAccess( ctxArray, contextPtr );

    for ( IndexType i = 0; i < 5; i++ )
    {
        BOOST_CHECK_EQUAL( static_cast<IndexType>( 3 ), readAccess.get()[i] );
    }

    readAccess.release();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( resizeTest )
{
    ContextPtr contextPtr = Context::getContextPtr( context::Host );

    HArray<IndexType> ctxArray; // default, not allocated at all
    {
        WriteAccess<IndexType> writeAccess( ctxArray, contextPtr );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( 10 );

        for ( IndexType i = 0; i < 10; i++ )
        {
            writeAccess.get()[i] = 3;
        }
    }
    ctxArray.purge();
    {
        WriteAccess<IndexType> writeAccess( ctxArray, contextPtr );
        // Possible problem: fetch from any location not possible
        writeAccess.resize( 10 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( asyncTest )
{
    ContextPtr hostContext = Context::getContextPtr( context::Host );
    ContextPtr cudaContext = Context::getContextPtr( context::CUDA );
    const IndexType n = 100;
    const float value = 1.4;
    const float alpha = 0.5;
    HArray<float> vector( n, value );
    common::shared_ptr<WriteAccess<float> > cudaV( new WriteAccess<float>( vector, cudaContext ) );

    common::shared_ptr<SyncToken> token( cudaContext->getSyncToken() );

    {
        SCAI_CONTEXT_ACCESS( cudaContext )
        scal( n, alpha, cudaV->get(), 1, token.get() );
        token->pushToken( cudaV );
    }

    cudaV.reset();  // give free the write access, but ownership also in token

    token->wait();  

    ReadAccess<float> hostV( vector, hostContext );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV.get()[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( syncTest )
{
    ContextPtr hostContext = Context::getContextPtr( context::Host );
    ContextPtr cudaContext = Context::getContextPtr( context::CUDA );

    const IndexType n = 100;
    const float value = 1.4;
    const float alpha = 0.5;
    HArray<float> vector( n, value );
    {
        WriteAccess<float> cudaV( vector, cudaContext );
        SCAI_CONTEXT_ACCESS( cudaContext );
        common::shared_ptr<SyncToken> token( cudaContext->getSyncToken() );
        scal( n, alpha, cudaV.get(), 1, token.get() );
        // synchronize on token at end of this scope
    }
    ReadAccess<float> hostV( vector, hostContext );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( value * alpha, hostV.get()[i] );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

