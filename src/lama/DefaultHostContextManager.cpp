/**
 * @file DefaultHostContextManager.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief DefaultHostContextManager.cpp
 * @author Thomas Brandes
 * @date 11.07.2011
 * $Id$
 */

// hpp
#include <lama/DefaultHostContextManager.hpp>

// others
#include <lama/DefaultHostContext.hpp>

#include <lama/exception/LAMAAssert.hpp>

// boost
#include <boost/bind.hpp>

using namespace boost;

namespace lama
{

// Explicit definition of static variable required.

DefaultHostContextManager DefaultHostContextManager::theInstance;

boost::weak_ptr<DefaultHostContext> DefaultHostContextManager::contextInstance;

/* ------------------------------------------------------------------------ */

void DefaultHostContextManager::setAsCurrent()
{
    LAMA_LOG_INFO( logger, "DefaultHostContextManager used for " << "ContextFactory::getContext( Context::Host )" )

    theInstance.registerFactory();
}

/* ------------------------------------------------------------------------ */

DefaultHostContextManager::DefaultHostContextManager()
    : ContextManager( Context::Host )
{
    // The singleton object adds itself to the device factory as manager for Host devices.

    registerFactory();
}

/* ------------------------------------------------------------------------ */

DefaultHostContextManager::~DefaultHostContextManager()
{
}

/* ------------------------------------------------------------------------ */

void DefaultHostContextManager::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "DefaultHostContextManager";
}

/* ------------------------------------------------------------------------ */

ContextPtr DefaultHostContextManager::getContext( int deviceNr )
{
    boost::shared_ptr<DefaultHostContext> context;

    if ( deviceNr != LAMA_DEFAULT_DEVICE_NUMBER )
    {
        LAMA_LOG_WARN( logger, "Context number ignored for HostContext, deviceNr = " << deviceNr )
    }

    // use the last contextInstance if it is still valid

    if ( contextInstance.expired() )
    {
        // create a new instance of DefaultHostContext and keep it for further uses

        context = boost::shared_ptr<DefaultHostContext>( new DefaultHostContext() );

        contextInstance = context;
    }
    else
    {
        // the last context instance is still valid, so we return new shared pointer to it

        context = contextInstance.lock();
    }

    return context;
}

} // namespace

