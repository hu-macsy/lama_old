/**
 * @file MockContext.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Definition of a Context class used for mock objects during tests.
 * @author Thomas Brandes
 * @date 05.07.2015
 */

#include "MockMemory.hpp"

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/Memory.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/macros/assert.hpp>

#include <memory>
#include <functional>

/** Exampes of a new context class that implements all relevant routines. */

class MockContext: public scai::hmemo::Context, public scai::hmemo::Context::Register<MockContext>

{
private:

    // Member variables

    int mDeviceNr;     // MockContext with different device numbers are not equal

    mutable std::weak_ptr<scai::hmemo::Memory> mMemory;

public:

    ~MockContext()
    {
        SCAI_LOG_DEBUG( logger, "~MockContext: " << *this )
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "MockContext( dev = " << mDeviceNr << " )";
    }

    virtual scai::hmemo::MemoryPtr getLocalMemoryPtr() const
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
        if ( memory.getType() != scai::hmemo::MemoryType::UserMemory )
        {
            return false;
        }

        const MockMemory* mockMemory = dynamic_cast<const MockMemory*>( &memory );

        SCAI_ASSERT( mockMemory, "NULL mock memory" )
        return mockMemory->getDeviceNr() == mDeviceNr;

        // return &memory == getMemory().get();
    }

    virtual scai::common::ContextType getType() const
    {
        return scai::common::ContextType::UserContext;
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

    static scai::common::ContextType createValue()
    {
        return scai::common::ContextType::UserContext;
    }

    bool isEqual( const Context& other ) const
    {
        if ( &other == this )
        {
            return true;
        }
    
        if ( other.getType() != scai::common::ContextType::UserContext )
        {
            return false;
        }
    
        const MockContext& otherMock = static_cast<const MockContext&>( other );

        return mDeviceNr == otherMock.mDeviceNr;
    }

private:

    // MockContext uses the type UserContext as its type

    MockContext( int deviceNr ) : scai::hmemo::Context( scai::common::ContextType::UserContext )
    {
        mDeviceNr = deviceNr;
    }
};

/* --------------------------------------------------------------------- */

static std::vector<std::weak_ptr<class MockContext> > contextInstances( 6 );

/* --------------------------------------------------------------------- */

inline scai::hmemo::ContextPtr MockContext::create( int deviceNr )
{
    std::shared_ptr<MockContext> context;
    SCAI_ASSERT( deviceNr < 6, "number of instances limited" )

    // use the last contextInstance if it is still valid

    if ( contextInstances[deviceNr].expired() )
    {
        // create a new instance of MockContext and keep it for further uses
        context = std::shared_ptr<MockContext>( new MockContext( deviceNr ) );
        contextInstances[deviceNr] = context;
    }
    else
    {
        // the last context instance is still valid, so we return new shared pointer to it
        context = contextInstances[deviceNr].lock();
    }

    return context;
}

