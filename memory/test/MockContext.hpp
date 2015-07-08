/**
 * @file MockContext.hpp
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
#include <memory/TaskSyncToken.hpp>

#include <boost/bind.hpp>
#include <boost/weak_ptr.hpp>

using namespace memory;

/** Exampes of a new context class that implements all relevant routines. */

class MockContext: 

     public Context, 
     public Context::Register<MockContext>

{
private: 

    // Member variables

    int mDeviceNr;     // MockContext with different device numbers are not equal

public:

    ~MockContext()
    {
        LAMA_LOG_DEBUG( logger, "~MockContext: " << *this )
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "MockContext( dev = " << mDeviceNr << " )";
    }

    virtual bool canUseData( const Context& other ) const
    {
        if ( other.getType() != context::UserContext ) 
        {
            return false;
        }

        const MockContext* otherMock = dynamic_cast<const MockContext*>( &other );
 
        COMMON_ASSERT( otherMock, "dynamic_cast<const MockContext*> failed" )

        return otherMock->mDeviceNr == mDeviceNr;   // only on same device
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

    virtual bool canCopyFrom( const Context& other ) const
    {
        // copy from host to this context should always be supported

        return other.getType() == context::Host;
    }

    virtual bool canCopyTo( const Context& other ) const
    {
        // copy from this context to host should always be supported

        return other.getType() == context::Host;
    }

    virtual void memcpyFrom( void* dst, const Context& srcContext, const void* src, size_t size ) const 
    {
        if ( srcContext.getType() == context::Host )
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
        if ( dstContext.getType() == context::Host )
        {
            ::memcpy( dst, src, size );
        }
        else
        {
            COMMON_THROWEXCEPTION( "copy to " << dstContext << " from " << *this << " not supported" )
        }
    }

    virtual TaskSyncToken* getSyncToken() const
    {
        return new TaskSyncToken();
    }

    /** Static method that delivers a MockContext for a certain device.
     *
     *  During the initialization this function will be registered at the base class Context
     */

    static ContextPtr create( int deviceNr );

    static ContextType createValue() 
    { 
        return context::UserContext; 
    }

private:

    // MockContext uses the type UserContext as its type

    MockContext( int deviceNr ) : Context( context::UserContext )
    {
        mDeviceNr = deviceNr;
    }
};

/* --------------------------------------------------------------------- */

static std::vector<boost::weak_ptr<class MockContext> > contextInstances( 6 );

/* --------------------------------------------------------------------- */

ContextPtr MockContext::create( int deviceNr )
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

