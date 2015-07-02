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

#include <memory/Context.hpp>
#include <memory/HostReadAccess.hpp>
#include <memory/HostWriteAccess.hpp>
#include <common/Exception.hpp>
#include <memory/TaskSyncToken.hpp>
#include <memory/Scalar.hpp>

#include <boost/weak_ptr.hpp>
#include <boost/bind.hpp>

using namespace memory;

/* --------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( logger, "ContextTest" )

/* --------------------------------------------------------------------- */

using memory::SyncToken;
using memory::TaskSyncToken;

/** Exampes of a new context class that implements all relevant routines. */

class MockContext: public Context
{
private: 

    // Member variables

    int mDeviceNr;     // MockContext with different device numbers are not equal

public:

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "Context( kind = User, dev = " << mDeviceNr << ")";
    }

    virtual bool canUseData( const Context& other ) const
    {
        if ( other.getType() != Context::UserContext ) 
        {
            return false;
        }

        const MockContext* otherMock = dynamic_cast<const MockContext*>( &other );
 
        COMMON_ASSERT( otherMock, "dynamic_cast<const MockContext*> failed" )

        return otherMock->mDeviceNr == mDeviceNr;   // only on same device
    }

    virtual ContextType getType() const
    {
        return Context::UserContext;
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
        COMMON_ASSERT_EQUAL( contextData.context->getType(), getType(), "type mismatch" )
        contextData.free();
    }

    virtual void memcpy( void* target, const void* source, const size_t size ) const
    {
        ::memcpy( target, source, size );
    }

    static SyncToken* theMemcpyAsync( void* dst, const void* src, const size_t size )
    {
        return new TaskSyncToken( boost::bind( &::memcpy, dst, src, size ) );
    }

    virtual SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const
    {
        return new TaskSyncToken( boost::bind( &::memcpy, dst, src, size ) );
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
        COMMON_ASSERT( cancpy( dst, src ), "Can not copy from " << * ( src.context ) << " to " << * ( dst.context ) );
        memcpy( dst.pointer, src.pointer, size );
    }

    virtual bool canCopyFrom( const Context& other ) const
    {
        // copy from host to this context should always be supported

        return other.getType() == Context::Host;
    }

    virtual bool canCopyTo( const Context& other ) const
    {
        // copy from this context to host should always be supported

        return other.getType() == Context::Host;
    }

    virtual void memcpyFrom( void* dst, const Context& srcContext, const void* src, size_t size ) const 
    {
        if ( srcContext.getType() == Context::Host )
        {
            ::memcpy( dst, src, size );
        }
        else
        {
            COMMON_THROWEXCEPTION( "copy from " << srcContext << " to " << *this << " not supported" )
        }
    }

    virtual void memcpyTo( const Context& dstContext, void* dst, const void* src, size_t size ) const 
    {
        if ( dstContext.getType() == Context::Host )
        {
            ::memcpy( dst, src, size );
        }
        else
        {
            COMMON_THROWEXCEPTION( "copy to " << dstContext << " from " << *this << " not supported" )
        }
    }

    virtual SyncToken* memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const
    {
        COMMON_ASSERT( cancpy( dst, src ), "Can not copy from " << * ( src.context ) << " to " << * ( dst.context ) );
        return memcpyAsync( dst.pointer, src.pointer, size );
    }

    virtual TaskSyncToken* getSyncToken() const
    {
        return new TaskSyncToken();
    }

    static ContextPtr getContext( int deviceNr );

private:

    // MockContext uses the type UserContext as its type

    MockContext( int deviceNr )
        : Context( Context::UserContext )
    {
        mDeviceNr = deviceNr;
    }

    static bool init();

    static bool initialized;   // initialization will register getContext as creator for Context
};

bool MockContext::init()
{
    Context::addCreator( Context::UserContext, &MockContext::getContext );
}

bool MockContext::initialized = MockContext::init();

static std::vector<boost::weak_ptr<class MockContext> > contextInstances( 6 );

ContextPtr MockContext::getContext( int deviceNr )
{
    boost::shared_ptr<MockContext> context;

    COMMON_ASSERT( deviceNr < 6, "number of instances limited" )

    // use the last contextInstance if it is still valid

    if( contextInstances[deviceNr].expired() )
    {
        // create a new instance of MockContext and keep it for further uses

        context = boost::shared_ptr<MockContext>( new MockContext( deviceNr ) );

        contextInstances[deviceNr] = context;
    }
    else
    {
        // the last context instance is still valid, so we return new shared pointer to it

        context = contextInstances[deviceNr].lock();
    }

    return context;
}

void releaseTest ()
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
    writeAccess.resize( 20 ); 
    // This is not checked:  writeAccess[0] = 5.0; -> crashes
    HostReadAccess<IndexType> readAccess( lamaArray );

    for ( IndexType i = 0; i < 5; i++ )
    {
        COMMON_ASSERT_EQUAL( 3, readAccess[i], "check" );
    }

    readAccess.release();
}

int main()
{
    // releaseTest();

    ContextPtr userContext  = Context::getContext( Context::UserContext, 1 );
    ContextPtr userContext2 = Context::getContext( Context::UserContext, 2 );
    ContextPtr hostContext  = Context::getContext( Context::Host );

    LAMA_LOG_INFO( logger, "userContext = " << *userContext );

    LAMAArray<double> X( 10, 5.0 );

    {
        WriteAccess<double> write( X, userContext );  
    }

    // read @ userContext2: valid data is transfered from userContext to here

    ReadAccess<double> read( X, userContext2 );

    const double* vals = read.get();

    for ( int i = 0; i < 10; ++i )
    {
        COMMON_ASSERT_EQUAL( vals[i], 5.0, "check" )
    }

    // Now make some checks

    std::cout << "X @ " << *userContext << ", valid = " << X.isValid( userContext )
              << ", capacity = " << X.capacity( userContext ) << std::endl;

    std::cout << "X @ " << *userContext2 << ", valid = " << X.isValid( userContext2 )
              << ", capacity = " << X.capacity( userContext2 ) << std::endl;

    std::cout << "X @ " << *hostContext << ", valid = " << X.isValid( hostContext )
              << ", capacity = " << X.capacity( hostContext ) << std::endl;

    LAMAArray<double> Y( X );

    // valid should be the same for Y, capacity should be 0 if not valid

    std::cout << "Y @ " << *userContext << ", valid = " << Y.isValid( userContext )
              << ", capacity = " << Y.capacity( userContext ) << std::endl;

    std::cout << "Y @ " << *userContext2 << ", valid = " << Y.isValid( userContext2 )
              << ", capacity = " << Y.capacity( userContext2 ) << std::endl;

    std::cout << "Y @ " << *hostContext << ", valid = " << Y.isValid( hostContext )
              << ", capacity = " << Y.capacity( hostContext ) << std::endl;

    Y.clear();

    std::cout << "Y cleared now" << std::endl;

    // valid should be the same for Y, capacity should be 0 if not valid

    std::cout << "Y @ " << *userContext << ", valid = " << Y.isValid( userContext )
              << ", capacity = " << Y.capacity( userContext ) << std::endl;

    std::cout << "Y @ " << *userContext2 << ", valid = " << Y.isValid( userContext2 )
              << ", capacity = " << Y.capacity( userContext2 ) << std::endl;

    std::cout << "Y @ " << *hostContext << ", valid = " << Y.isValid( hostContext )
              << ", capacity = " << Y.capacity( hostContext ) << std::endl;

    Y.purge();

    std::cout << "Y purged now" << std::endl;

    // valid should be the same for Y, capacity should be 0 if not valid

    std::cout << "Y @ " << *userContext << ", valid = " << Y.isValid( userContext )
              << ", capacity = " << Y.capacity( userContext ) << std::endl;

    std::cout << "Y @ " << *userContext2 << ", valid = " << Y.isValid( userContext2 )
              << ", capacity = " << Y.capacity( userContext2 ) << std::endl;

    std::cout << "Y @ " << *hostContext << ", valid = " << Y.isValid( hostContext )
              << ", capacity = " << Y.capacity( hostContext ) << std::endl;

    int values[] = { 1, 2, 3, 4 };

    LAMAArray<float> v ( 4, values );   // implicit type conversion allowed
}

