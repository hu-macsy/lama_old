/**
 * @file PGASContextManager.cpp
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
 * @brief Implementation of PGAS host manager class methods.
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */

#include <lama/pgas/PGASContextManager.hpp>
#include <lama/pgas/PGASContext.hpp>

#include <lama/exception/LAMAAssert.hpp>

namespace lama
{

// Explicit definition of static variable required.

PGASContextManager PGASContextManager::theInstance;

boost::weak_ptr<const PGASContext> PGASContextManager::contextInstance;
ContextPtr PGASContextManager::getContext( int )
{
    return getContext();
}

/* ------------------------------------------------------------------------ */

PGASContextManager::PGASContextManager()
    : ContextManager( Context::Host )
{
    // This manager will NOT register itself as default, must be done explicitly
    registerFactory();
}

/* ------------------------------------------------------------------------ */

PGASContextManager::~PGASContextManager()
{
}

/* ------------------------------------------------------------------------ */

ContextPtr PGASContextManager::getContext()
{
    boost::shared_ptr<const PGASContext> context;

    // use the last contextInstance if it is still valid

    if ( contextInstance.expired() )
    {
        // create a new instance of PGASContext and keep it for further uses
        context = boost::shared_ptr<const PGASContext>( new PGASContext() );
        contextInstance = context;
    }
    else
    {
        // the last context instance is still valid, so we return new shared pointer to it
        context = contextInstance.lock();
    }

    return context;
}

void PGASContextManager::setAsCurrent()
{
    theInstance.registerFactory();
}

} // namespace lama

