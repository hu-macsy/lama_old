/**
 * @file ContextFactory.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Definition of a factory class that provides device pointers for different
 *        device types via device manager.
 * @author Jiri Kraus
 * @date 10.07.2011
 * @since 1.0.0
 */

// hpp
#include <lama/ContextFactory.hpp>

// others
#include <lama/ContextManager.hpp>

#include <lama/exception/LAMAAssert.hpp>

using namespace std;
using namespace boost;

namespace lama
{

LAMA_LOG_DEF_LOGGER( ContextFactory::logger, "ContextFactory" )

/** Set NULL manager as defaults for all devices */

const char* ContextFactory::theContextIds[] =
{   "HostContext", "CudaContext", "OpenCLContext"
//                                                , "NewContext"
};

/* --------------------------------------------------------------------------------- */

// Static variables of Context Factory

ContextFactory* ContextFactory::theContextFactory = NULL;

/* --------------------------------------------------------------------------------- */

ContextFactory::ContextFactory()
{
    // initialize the context managers for this factory

    for ( int type = 0; type < int( Context::MaxContext ); ++type )
    {
        mContextManager[ContextType( type )] = NULL;
    }
}

/* --------------------------------------------------------------------------------- */

ContextFactory::~ContextFactory()
{
}

/* --------------------------------------------------------------------------------- */

ContextFactory& ContextFactory::getFactory()
{
    if ( !theContextFactory )
    {
        theContextFactory = new ContextFactory();
    }

    return *theContextFactory;
}

/* --------------------------------------------------------------------------------- */

void ContextFactory::addContextManager( const Context::ContextType type, ContextManager& deviceManager )
{
    LAMA_ASSERT_DEBUG( type < Context::MaxContext, "illegal device" )

    // be careful about logging, method might be called during static initialization

    if ( mContextManager[type] )
    {
        // ContextManager might be replaced, e.g. for Host context (DefaultHost or CUDAHost)

        LAMA_LOG_INFO( logger, "Context manager replaced for context type " << theContextIds[type] )
    }

    mContextManager[type] = &deviceManager;
}

/* --------------------------------------------------------------------------------- */

ContextPtr ContextFactory::getContext( Context::ContextType type, int deviceNr )
{
    LAMA_ASSERT_DEBUG( type < Context::MaxContext, "illegal device type = " << type << ", out of range" )

    LAMA_LOG_DEBUG( logger, "getContext: type = " << theContextIds[type] << ", nr = " << deviceNr )

    const ContextFactory& factory = getFactory();

    ContextManager* manager = factory.mContextManager[type];

    if ( !manager )
    {
        LAMA_THROWEXCEPTION( "Context " << theContextIds[type] << " not supported: no device manager available" )
    }

    ContextPtr context = manager->getContext( deviceNr );

    LAMA_ASSERT_DEBUG( context, "Context manager returned NULL context, should not happen" )

    if ( deviceNr == LAMA_DEFAULT_DEVICE_NUMBER )
    {
        LAMA_LOG_DEBUG( logger, "getContext( type =" << theContextIds[type] << " ) => " << *context )
    }
    else
    {
        LAMA_LOG_DEBUG( logger, "getContext( type = " << theContextIds[type] 
                                 << ", nr = " << deviceNr << ") => " << *context )
    }

    return context;
}

/* --------------------------------------------------------------------------------- */

ContextManager* ContextFactory::getContextManager( Context::ContextType type )
{
    LAMA_ASSERT_DEBUG( type < Context::MaxContext, "illegal device type = " << type << ", out of range" )
    return mContextManager[type];
}

/* --------------------------------------------------------------------------------- */

bool ContextFactory::hasContext( const Context::ContextType type )
{
    const ContextFactory& factory = getFactory();
    return factory.mContextManager[type] != 0;
}

/* ---- release  ------------------------------------------------------------------- */

void ContextFactory::release()
{
    // release/free available context managers on the factory instance

    ContextFactory& factory = getFactory();

    LAMA_LOG_DEBUG( logger, "release context managers" )

    for ( int type = 0; type < int( Context::MaxContext ); ++type )
    {
        factory.mContextManager[ContextType( type )] = NULL;
    }

    LAMA_LOG_INFO( logger, "released all context managers" )
}

}
