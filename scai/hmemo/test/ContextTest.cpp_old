/**
 * @file ContextTest.cpp
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
 * @brief Contains the implementation of the class ContextTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 01.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/thread.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/lama/ContextManager.hpp>
#include <scai/lama/ContextFactory.hpp>
#include <scai/lama/HostReadAccess.hpp>
#include <scai/lama/HostWriteAccess.hpp>
#include <scai/lama/task/TaskSyncToken.hpp>
#include <scai/lama/exception/LAMAAssert.hpp>

#include <scai/common/unique_ptr.hpp>
#include <scai/common/bind.hpp>

#include <scai/tasking/SyncToken.hpp>
#include <scai/tasking/TaskSyncToken.hpp>

#include <omp.h>

using namespace scai::lama;
using namespace scai::hmemo;
using scai::tasking::SyncToken;
using scai::tasking::TaskSyncToken;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE ( ContextTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.ContextTest" )

/* --------------------------------------------------------------------- */

class MockContext: public Context
{
    friend class MockContextManager;

public:

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "MockContext";
    }

    virtual bool canUseData( const Context& other ) const
    {
        return other.getType() == getType();
    }

    virtual ContextType getType() const
    {
        return context::UserContext;
    }

    virtual void* allocate( const size_t size ) const
    {
        return malloc( size );
    }

    virtual void free( void* pointer, const size_t ) const
    {
        ::free( pointer );
    }

    virtual void allocate( ContextData& contextData, const size_t size ) const
    {
        contextData.pointer = allocate( size );
    }

    virtual void free( ContextData& contextData ) const
    {
        SCAI_ASSERT_EQUAL_ERROR( contextData.context->getType(), getType() );
        BOOST_CHECK_EQUAL( contextData.context->getType(), getType() );
        contextData.free();
    }

    virtual void memcpy( void* target, const void* source, const size_t size ) const
    {
        ::memcpy( target, source, size );
    }

    static scai::common::unique_ptr<SyncToken> theMemcpyAsync( void* dst, const void* src, const size_t size )
    {
        return scai::common::unique_ptr < SyncToken > ( new TaskSyncToken( scai::common::bind( &::memcpy, dst, src, size ) ) );
    }

    virtual scai::common::unique_ptr<SyncToken> memcpyAsync( void* dst, const void* src, const size_t size ) const
    {
        return scai::common::unique_ptr < SyncToken > ( new TaskSyncToken( scai::common::bind( &::memcpy, dst, src, size ) ) );
    }

    virtual bool cancpy( const ContextData& dst, const ContextData& src ) const
    {
        return ( dst.context->getType() == getType() && src.context->getType() == getType() )
               || ( dst.context->getType() == Context::Host && src.context->getType() == getType() )
               || ( dst.context->getType() == getType() && src.context->getType() == Context::Host )
               || ( dst.context->getType() == Context::Host && src.context->getType() == Context::Host );
    }

    virtual void memcpy( ContextData& dst, const ContextData& src, const size_t size ) const
    {
        SCAI_ASSERT_ERROR( cancpy( dst, src ), "Can not copy from " << * ( src.context ) << " to " << * ( dst.context ) );
        memcpy( dst.pointer, src.pointer, size );
    }

    virtual scai::common::unique_ptr<SyncToken> memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const
    {
        SCAI_ASSERT_ERROR( cancpy( dst, src ), "Can not copy from " << * ( src.context ) << " to " << * ( dst.context ) );
        return memcpyAsync( dst.pointer, src.pointer, size );
    }

    virtual scai::common::unique_ptr<SyncToken> getSyncToken() const
    {
        return scai::common::unique_ptr < SyncToken > ( new TaskSyncToken() );
    }

private:

    // MockContext uses the type NewContext as its type

    MockContext()
        : Context( context::UserContext )
    {
    }
};

/* ---------------------------------------------------------------------- */

class MockContextManager: public ContextManager
{
public:

    ContextPtr getContext( int /* context */ )
    {
        return getInstance();
    }

    ~MockContextManager()
    {
    }

    static ContextPtr getInstance()
    {
        if ( !mMockContext )
        {
            mMockContext = ContextPtr( new MockContext() );
        }

        return mMockContext;
    }

private:

    MockContextManager()
        : ContextManager( Context::NewContext )
    {
        registerFactory();
    }

    static ContextPtr mMockContext; // the only one context managed

    static MockContextManager theInstance; // singleton instance of this class
};

ContextPtr MockContextManager::mMockContext = ContextPtr(); // will be allocated later

MockContextManager MockContextManager::theInstance;

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( releaseTest )
{
    LAMAArray<IndexType> lamaArray; // default, not allocated at all
    // read access on empty array should work even if not useful
    HostReadAccess<IndexType> readTestAccess( lamaArray );
    // release read on empty array
    readTestAccess.release();
    // get write access on empty array
    HostWriteAccess<IndexType> writeAccess( lamaArray );
    writeAccess.resize( 10 );

    for ( IndexType i = 0; i < 10; i++ )
    {
        writeAccess[i] = 3;
    }

    writeAccess.release();
    // working on a released array should give an exception
    SCAI_CHECK_THROW( { writeAccess.resize( 20 ); }, Exception );
    // This is not checked:  writeAccess[0] = 5.0; -> crashes
    HostReadAccess<IndexType> readAccess( lamaArray );

    for ( IndexType i = 0; i < 5; i++ )
    {
        BOOST_CHECK_EQUAL( 3, readAccess[i] );
    }

    readAccess.release();
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTest )
{
    LAMAArray<IndexType> data; // default, not allocated at all
    ContextPtr context = ContextFactory::getContext( Context::NewContext );
    {
        WriteAccess<IndexType> arr( data, context );
        arr.resize( 10 );
        IndexType* idata = arr.get();

        for ( IndexType i = 0; i < 10; i++ )
        {
            idata[i] = 23;
        }
    }
    {
        HostWriteAccess<IndexType> arr( data );
        BOOST_CHECK_EQUAL( 23, arr[8] );
        arr[5] = 13;
        arr.resize( 20 );
        arr[17] = 34;
    }
    {
        ReadAccess<IndexType> arr( data, context );
        const IndexType* idata = arr.get();
        BOOST_CHECK_EQUAL( 13, idata[5] );
        BOOST_CHECK_EQUAL( 23, idata[8] );
        BOOST_CHECK_EQUAL( 34, idata[17] );
    }
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    LAMAArray<IndexType> data; // default, not allocated at all
    ContextPtr context = ContextFactory::getContext( Context::NewContext );
    WriteAccess<IndexType> arr( data, context );
    arr.resize( 10 );
    IndexType* idata = arr.get();

    for ( IndexType i = 0; i < 10; i++ )
    {
        idata[i] = 23;
    }

    // the copy operator should also work on a write-locked array
    LAMAArray<IndexType> data1( data );
    BOOST_CHECK_EQUAL( data.size(), data1.size() );
    // check that copied array has valid data
    {
        HostReadAccess<IndexType> arr( data1 );
        BOOST_CHECK_EQUAL( 23, arr[8] );
    }
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( prefetchTest )
{
    LAMAArray<IndexType> data;
    ContextPtr context = ContextFactory::getContext( Context::NewContext );
    {
        HostWriteAccess<IndexType> arr( data );
        arr.resize( 10 );

        for ( IndexType i = 0; i < 10; i++ )
        {
            arr[i] = 10 + i;
        }
    }
    data.prefetch( context );
    {
        ReadAccess<IndexType> arr( data, context );
        const IndexType* data = arr.get();
        BOOST_CHECK_EQUAL( 12, data[2] );
        BOOST_CHECK_EQUAL( 19, data[9] );
    }
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( swapTest )
{
    LAMAArray<IndexType> data1;
    LAMAArray<IndexType> data2;
    ContextPtr context = ContextFactory::getContext( Context::NewContext );
    {
        WriteAccess<IndexType> arr( data1, context );
        arr.resize( 10 );
        IndexType* idata = arr.get(); // no indexing allowed for arbitrary write context

        for ( IndexType i = 0; i < 10; i++ )
        {
            idata[i] = 10 + i;
        }
    }
    {
        HostWriteOnlyAccess<IndexType> arr( data2, 10 );

        for ( IndexType i = 0; i < 10; i++ )
        {
            arr[i] = 100 - i;
        }

        // swapping of arrays not allowed with any access
        SCAI_CHECK_THROW( data1.swap( data2 ), Exception );
    }
    data1.swap( data2 );
    {
        // data1 contains now original data2
        HostReadAccess<IndexType> arr( data1 );
        BOOST_CHECK_EQUAL( 100, arr[0] );
        BOOST_CHECK_EQUAL( 91, arr[9] );
    }
    {
        HostReadAccess<IndexType> arr( data2 );
        BOOST_CHECK_EQUAL( 12, arr[2] );
        BOOST_CHECK_EQUAL( 19, arr[9] );
    }
}

/* ---------------------------------------------------------------------- */

const IndexType N = 300;
const IndexType ITER = 10;

void sumit( IndexType& sum, const LAMAArray<IndexType>& data )
{
    for ( IndexType i = 0; i < N; i++ )
    {
        HostReadAccess<IndexType> arr( data );
        sum += arr[i];
    }
}

BOOST_AUTO_TEST_CASE( threadSafetyTest )
{
    LAMAArray<IndexType> data;
    // define the array with some data
    {
        HostWriteAccess<IndexType> arr( data );
        arr.resize( N );

        for ( IndexType i = 0; i < N; i++ )
        {
            arr[i] = 1;
        }
    }
    IndexType sum1 = 0;
    IndexType sum2 = 0;
    IndexType sum3 = 0;

    for ( IndexType k = 0; k < ITER; k++ )
    {
        // run this process + two additional threads
        boost::thread t1( common::bind( sumit, common::ref( sum1 ), common::ref( data ) ) );
        boost::thread t2( common::bind( sumit, common::ref( sum2 ), common::ref( data ) ) );
        sumit( sum3, data );
        t1.join();
        t2.join();
    }

    BOOST_CHECK_EQUAL( sum1, N * ITER );
    BOOST_CHECK_EQUAL( sum2, N * ITER );
    BOOST_CHECK_EQUAL( sum3, N * ITER );
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ompSafetyTest )
{
    LAMAArray<IndexType> data;
    LAMAArray<IndexType> data1;
    // define the array with some data
    {
        HostWriteAccess<IndexType> arr( data );
        arr.resize( N );

        for ( IndexType i = 0; i < N; i++ )
        {
            arr[i] = 1;
        }
    }
    // now create some OpenMP threads that take read access
    omp_set_num_threads( 2 );
    #pragma omp parallel
    {
        IndexType sum = 0;

        for ( IndexType i = 0; i < N; i++ )
        {
            HostReadAccess<IndexType> arr( data );
            sum = sum + arr[i];
        }

        BOOST_CHECK_EQUAL( N, sum );
    }
    // no more access, so we can swap
    data.swap( data1 );
    #pragma omp parallel
    {
        IndexType sum = 0;

        for ( IndexType i = 0; i < N; i++ )
        {
            HostReadAccess<IndexType> arr( data1 );
            sum = sum + arr[i];
        }

        BOOST_CHECK_EQUAL( N, sum );
    }
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    ContextPtr context = ContextFactory::getContext( Context::Host );
    LAMA_WRITEAT_PTR_TEST( context );
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( EqualityTest )
{
    ContextPtr contextA = ContextFactory::getContext( Context::Host );
    ContextPtr contextB = ContextFactory::getContext( Context::Host );
    ContextPtr contextC = ContextFactory::getContext( Context::NewContext );
    BOOST_CHECK( contextA == contextB );
    BOOST_CHECK( contextA != contextC );
}
/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
