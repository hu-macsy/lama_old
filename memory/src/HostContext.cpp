/**
 * @file HostContext.cpp
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
 * @brief HostContext.cpp
 * @author Thomas Brandes
 * @date 11.07.2011
 */

// hpp
#include <memory/HostContext.hpp>
#include <memory/HostMemory.hpp>

// others
#include <common/Exception.hpp>
#include <common/OpenMP.hpp>
#include <common/weak_ptr.hpp>

#include <tasking/TaskSyncToken.hpp>

using  common::shared_ptr;
using  common::weak_ptr;

namespace memory
{

LAMA_LOG_DEF_LOGGER( HostContext::logger, "Context.HostContext" )

HostContext::HostContext() : Context( context::Host )
{
    LAMA_LOG_INFO( logger, "HostContext created" )
}

HostContext::~HostContext()
{
    LAMA_LOG_INFO( logger, "~HostContext" )
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
        LAMA_LOG_WARN( logger, "Context number ignored for HostContext, deviceNr = " << deviceNr )
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

    // same object by pointer can always use same data.

    if ( other.getType() == memtype::HostMemory )
    {
        canUseIt = true;
    }

    return canUseIt;
}

/* ------------------------------------------------------------------------- */

MemoryPtr HostContext::getMemory() const
{
    return HostMemory::getIt();
}

/* ------------------------------------------------------------------------- */

tasking::TaskSyncToken* HostContext::getSyncToken() const
{
    // on Host we will run asynchronous computations as a task

    return new tasking::TaskSyncToken();
}

} // namespace

