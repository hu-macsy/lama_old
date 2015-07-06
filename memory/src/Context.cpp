/**
 * @file Context.cpp
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
 * @brief Context.cpp
 * @author Thomas Brandes
 * @date 11.07.2011
 * @since 1.0.0
 */

// hpp
#include <memory/Context.hpp>

// assert
#include <common/Exception.hpp>

#include <map>

namespace memory
{

LAMA_LOG_DEF_LOGGER( Context::logger, "Context" )

Context::Context( ContextType type )
                : mContextType( type ), mEnabled( false ), mFile( NULL ), mLine( 0 )
{
    LAMA_LOG_DEBUG( logger, "Context( type = " << mContextType << " )" )
}

Context::~Context()
{
    LAMA_LOG_DEBUG( logger, "~Context( type = " << mContextType << " )" )
}

void Context::writeAt( std::ostream& stream ) const
{
    // write identification of this object

    stream << "Context";
}

std::ostream& operator<<( std::ostream& stream, const ContextType& type )
{
    if( type == Context::Host )
    {
        stream << "Host";
    }
    else if( type == Context::CUDA )
    {
        stream << "CUDA";
    }
    else if( type == Context::OpenCL )
    {
        stream << "OpenCL";
    }
    else if( type == Context::MIC )
    {
        stream << "MIC";
    }
    else if( type == Context::UserContext )
    {
        stream << "UserContext";
    }
    else if( type == Context::UserContext1 )
    {
        stream << "UserContext1";
    }
    else
    {
        stream << "Context_" << ( int ) type;
    }

    return stream;
}

/* -----------------------------------------------------------------------------*/

void Context::enable( const char* file, int line ) const
{
    LAMA_LOG_INFO( logger, file << "( line = " << line << ") : enable " << *this )

    if( mEnabled )
    {
        LAMA_LOG_INFO( logger,
                       "Enable " << *this << " at " << file << " ( line = " << line << " )" << ", already enabled at " << mFile << " ( line = " << mLine << " )" )
    }

    mEnabled = true;
    mFile = file;
    mLine = line;
}

/* -----------------------------------------------------------------------------*/

void Context::disable( const char* file, int line ) const
{
    LAMA_LOG_INFO( logger, file << "( line = " << line << ") : disable " << *this )

    if( !mEnabled )
    {
        LAMA_LOG_INFO( logger,
                       "Disable " << *this << " at " << file << " ( line = " << line << " )" << ", context was not enabled before" )
    }

    mEnabled = false;
    mFile = NULL;
    mLine = 0;
}

/* ---------------------------------------------------------------------------------*/

ContextPtr Context::getHostContext() const
{
    // take the default host context of the context factory

    return getContext( Context::Host );
}

/* ---------------------------------------------------------------------------------*/

/** Map container to get for the key the create function. */

typedef std::map<Context::ContextType, Context::CreateFn> CreatorMap;

static CreatorMap& getFactory()
{
    static std::auto_ptr<CreatorMap> factory;

    if( !factory.get() )
    {
        factory = std::auto_ptr<CreatorMap>( new CreatorMap() );
    }

    return *factory;
}

void Context::addCreator( const Context::ContextType type,  CreateFn create )
{
    CreatorMap& factory = getFactory();

    // checks for multiple entries is not really necessary here, so just add entry in map container.

    factory[type] = create;
}

ContextPtr Context::getContext( Context::ContextType type, int deviceNr )
{
    ContextPtr context;

    COMMON_ASSERT( type < Context::MaxContext, "illegal device type = " << type << ", out of range" )

    const CreatorMap& factory = getFactory();

    CreatorMap::const_iterator fn = factory.find( type );

    if ( fn != factory.end() )
    {
        context = fn->second( deviceNr );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Context " << type << " not supported: no creator available" )
    }

    return context;
}

}
