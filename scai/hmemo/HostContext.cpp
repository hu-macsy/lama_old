/**
 * @file HostContext.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of methods for class HostContext.
 * @author Thomas Brandes
 * @date 11.07.2011
 */

// hpp
#include <scai/hmemo/HostContext.hpp>

// local library
#include <scai/hmemo/HostMemory.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/weak_ptr.hpp>

namespace scai
{

using  common::shared_ptr;
using  common::weak_ptr;

namespace hmemo
{

SCAI_LOG_DEF_LOGGER( HostContext::logger, "Context.HostContext" )

HostContext::HostContext() : Context( common::context::Host )
{
    SCAI_LOG_INFO( logger, "HostContext created" )
}

HostContext::~HostContext()
{
    SCAI_LOG_INFO( logger, "~HostContext" )
}

void HostContext::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    int nThreads = 1;

#pragma omp parallel
    {
#pragma omp master
        {
            nThreads = omp_get_num_threads();
        }
    }

    stream << "HostContext( #Threads = " << nThreads << " )";
}

static weak_ptr<class HostContext> contextInstance;

ContextPtr HostContext::create( int deviceNr )
{
    shared_ptr<HostContext> context;

    if ( deviceNr >= 0 )
    {
        SCAI_LOG_WARN( logger, "Context number ignored for HostContext, deviceNr = " << deviceNr )
    }

    // use the last contextInstance if it is still valid

    if ( contextInstance.expired() )
    {
        // create a new instance of HostContext and keep it for further uses

        context.reset( new HostContext() );

        contextInstance = context;
    }
    else
    {
        // the last context instance is still valid, so we return new shared pointer to it

        context = contextInstance.lock();
    }

    return context;
}

/* ------------------------------------------------------------------------- */

bool HostContext::canUseMemory( const Memory& other ) const
{
    bool canUseIt = false;

    if ( other.getType() == memtype::HostMemory )
    {
        canUseIt = true;
    }
    else if ( other.getContext().getType() == common::context::Host )
    {
        // If other memory can be used on Host it is okay

        canUseIt = true;
    }

    return canUseIt;
}

/* ------------------------------------------------------------------------- */

MemoryPtr HostContext::getLocalMemoryPtr() const
{
    return HostMemory::getIt();
}

/* ------------------------------------------------------------------------- */

tasking::TaskSyncToken* HostContext::getSyncToken() const
{
    // on Host we will run asynchronous computations as a task

    return new tasking::TaskSyncToken();
}

} /* end namespace hmemo */

} /* end namespace scai */
