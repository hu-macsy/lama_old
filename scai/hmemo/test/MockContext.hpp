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
 * @brief Definition of a Context class used for mock objects during tests.
 * @author: Thomas Brandes
 * @date 05.07.2015
 **/

#include "MockMemory.hpp"

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/Memory.hpp>
 
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/weak_ptr.hpp>
#include <scai/common/macros/assert.hpp>

/** Exampes of a new context class that implements all relevant routines. */

class MockContext: public scai::hmemo::Context, public scai::hmemo::Context::Register<MockContext>

{
private: 

    // Member variables

    int mDeviceNr;     // MockContext with different device numbers are not equal

    mutable scai::common::weak_ptr<scai::hmemo::Memory> mMemory;

public:

    ~MockContext()
    {
        SCAI_LOG_DEBUG( logger, "~MockContext: " << *this )
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "MockContext( dev = " << mDeviceNr << " )";
    }

    virtual scai::hmemo::MemoryPtr getMemoryPtr() const
    {
        scai::hmemo::MemoryPtr memory;

        if ( mMemory.expired() )
        {
            memory.reset( new MockMemory( mDeviceNr ) );
            mMemory = memory;
        }
        else
        {
            memory = mMemory.lock();
        }

        return memory;
    }

    virtual bool canUseMemory( const scai::hmemo::Memory& memory ) const
    {
        SCAI_ASSERT( &memory, "NULL memory" )

        if ( memory.getType() != scai::hmemo::memtype::UserMemory )
        {
            return false;
        }

        const MockMemory* mockMemory = dynamic_cast<const MockMemory*>( &memory );

        SCAI_ASSERT( mockMemory, "NULL mock memory" )

        return mockMemory->getDeviceNr() == mDeviceNr;
        
        // return &memory == getMemory().get();
    }

    virtual scai::common::context::ContextType getType() const
    {
        return scai::common::context::UserContext;
    }

    virtual scai::tasking::TaskSyncToken* getSyncToken() const
    {
        return new scai::tasking::TaskSyncToken();
    }

    /** Static method that delivers a MockContext for a certain device.
     *
     *  During the initialization this function will be registered at the base class Context
     */

    static scai::hmemo::ContextPtr create( int deviceNr );

    static scai::common::context::ContextType createValue()
    { 
        return scai::common::context::UserContext;
    }

private:

    // MockContext uses the type UserContext as its type

    MockContext( int deviceNr ) : scai::hmemo::Context( scai::common::context::UserContext )
    {
        mDeviceNr = deviceNr;
    }
};

/* --------------------------------------------------------------------- */

static std::vector<scai::common::weak_ptr<class MockContext> > contextInstances( 6 );

/* --------------------------------------------------------------------- */

inline scai::hmemo::ContextPtr MockContext::create( int deviceNr )
{
    scai::common::shared_ptr<MockContext> context;

    SCAI_ASSERT( deviceNr < 6, "number of instances limited" )

    // use the last contextInstance if it is still valid

    if( contextInstances[deviceNr].expired() )
    {
        // create a new instance of MockContext and keep it for further uses

        context = scai::common::shared_ptr<MockContext>( new MockContext( deviceNr ) );

        contextInstances[deviceNr] = context;
    }
    else
    {
        // the last context instance is still valid, so we return new shared pointer to it

        context = contextInstances[deviceNr].lock();
    }

    return context;
}

