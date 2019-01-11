/**
 * @file Context.cpp
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
 * @brief Implementation of methods for the Context base class and its factory.
 * @author Thomas Brandes
 * @date 10.07.2015
 */

// hpp
#include <scai/hmemo/Context.hpp>

// local library
#include <scai/hmemo/Memory.hpp>

// internal scai libraries
#include <scai/common/macros/throw.hpp>
#include <scai/common/Settings.hpp>

// std
#include <map>
#include <cctype>

namespace scai
{

using common::ContextType;

namespace hmemo
{

/* -----------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( Context::logger, "Context" )

Context::Context( ContextType type ) : 

    mContextType( type ), 
    mUseZeroCopy( false ), 
    mEnabled( false ), 
    mFile( NULL ), 
    mLine( 0 )
{
    SCAI_LOG_DEBUG( logger, "Context( type = " << mContextType << " )" )
}

Context::~Context()
{
    SCAI_LOG_DEBUG( logger, "~Context( type = " << mContextType << " )" )
}

void Context::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "Context";
}

/* -----------------------------------------------------------------------------*/

void Context::enable( const char* file, int line ) const
{
    SCAI_LOG_INFO( logger, file << "( line = " << line << ") : enable " << *this )

    if ( mEnabled )
    {
        SCAI_LOG_INFO( logger,
                       "Enable " << *this << " at " << file << " ( line = " << line << " )" << ", already enabled at " << mFile << " ( line = " << mLine << " )" )
    }

    mEnabled = true;
    mFile = file;
    mLine = line;

    setCurrent();
}

/* -----------------------------------------------------------------------------*/

void Context::disable( const char* file, int line ) const
{
    SCAI_LOG_INFO( logger, file << "( line = " << line << ") : disable " << *this )

    if ( !mEnabled )
    {
        SCAI_LOG_INFO( logger,
                       "Disable " << *this << " at " << file << " ( line = " << line << " )" << ", context was not enabled before" )
    }

    mEnabled = false;
    mFile = NULL;
    mLine = 0;

    unsetCurrent();
}

/* ---------------------------------------------------------------------------------*/

void Context::enableZeroCopy( bool flag ) const
{
    if ( flag == mUseZeroCopy )
    {
        return;
    }

    if ( flag )
    {
        if ( !this->canUseMemory( *this->getHostMemoryPtr() ) )
        {
            SCAI_LOG_WARN( logger, *this << ": does not support zero copy" )
            return;
        }
    }

    mUseZeroCopy = flag;
}

/* ---------------------------------------------------------------------------------*/

MemoryPtr Context::getHostMemoryPtr() const
{
    // take the host memory of the memory factory
    ContextPtr hostContextPtr = Context::getContextPtr( common::ContextType::Host );
    return hostContextPtr->getMemoryPtr();
}

/* ---------------------------------------------------------------------------------*/

MemoryPtr Context::getMemoryPtr() const
{
    if ( mUseZeroCopy )
    {
        return getHostMemoryPtr();
    }
    else
    {
        return getLocalMemoryPtr();
    }
}

/* ---------------------------------------------------------------------------------*/

ContextPtr Context::getContextPtr( ContextType type, int deviceNr )
{
    return create( type, deviceNr );
}

/* ---------------------------------------------------------------------------------*/

ContextPtr Context::getContextPtr()
{
    std::string ctx_string;

    if ( common::Settings::getEnvironment( ctx_string, "SCAI_CONTEXT" ) )
    {
        // ctx_string name not case sensitive, take it upper case
        for ( std::string::iterator p = ctx_string.begin(); ctx_string.end() != p; ++p )
        {
            *p = static_cast<std::string::value_type>( toupper( *p ) );
        }

        if ( ctx_string == "CUDA" )
        {
            return getContextPtr( common::ContextType::CUDA );
        }

        if ( ctx_string == "HOST" )
        {
            return getContextPtr( common::ContextType::Host );
        }

        COMMON_THROWEXCEPTION( "SCAI_CONTEXT=" << ctx_string << ", unknown context type" )
    }

    return getHostPtr();
}

/* ---------------------------------------------------------------------------------*/
/*  thread private : current context                                                */
/* ---------------------------------------------------------------------------------*/

void Context::setCurrent() const
{
    contextStack.push( this );
}

void Context::unsetCurrent() const
{
    if ( contextStack.empty() )
    {
        SCAI_LOG_WARN( logger, "unset this context " << *this << " but not set before" )
        return;
    }

    const Context* current = contextStack.top();

    if ( current != this )
    {
        SCAI_LOG_WARN( logger, "unset this context " << *this << " but was: " << *current )
    }
    else
    {
        SCAI_LOG_INFO( logger, "unset this context " << *this << " was current" )
    }

    contextStack.pop();
}

const Context* Context::getCurrentContext()
{
    if ( contextStack.empty() )
    {
        return NULL;
    }

    return contextStack.top();
}

/** Important: Each thread has its own context stack */

thread_local Context::ContextStack Context::contextStack;

/* ----------------------------------------------------------------------- */

} /* end namespace hmemo */

} /* end namespace scai */

/* ---------------------------------------------------------------------------------*/
