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

namespace scai
{

using common::context::ContextType;

namespace hmemo
{

/* -----------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( Context::logger, "Context" )

Context::Context( ContextType type )
                : mContextType( type ), mUseZeroCopy( false ), mEnabled( false ), mFile( NULL ), mLine( 0 )
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

    if( mEnabled )
    {
        SCAI_LOG_INFO( logger,
                       "Enable " << *this << " at " << file << " ( line = " << line << " )" << ", already enabled at " << mFile << " ( line = " << mLine << " )" )
    }

    mEnabled = true;
    mFile = file;
    mLine = line;
}

/* -----------------------------------------------------------------------------*/

void Context::disable( const char* file, int line ) const
{
    SCAI_LOG_INFO( logger, file << "( line = " << line << ") : disable " << *this )

    if( !mEnabled )
    {
        SCAI_LOG_INFO( logger,
                       "Disable " << *this << " at " << file << " ( line = " << line << " )" << ", context was not enabled before" )
    }

    mEnabled = false;
    mFile = NULL;
    mLine = 0;
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

    ContextPtr hostContextPtr = Context::getContextPtr( common::context::Host );
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
            *p = toupper( *p );
        }

        if ( ctx_string == "CUDA" )
        {
            return getContextPtr( common::context::CUDA );
        }

        if ( ctx_string == "MIC" )
        {
            return getContextPtr( common::context::MIC );
        }

        if ( ctx_string == "HOST" )
        {
            return getContextPtr( common::context::Host );
        }

        COMMON_THROWEXCEPTION( "SCAI_CONTEXT=" << ctx_string << ", unknown communicator type" )
    }

    return getHostPtr();
}

} /* end namespace hmemo */

} /* end namespace scai */

/* ---------------------------------------------------------------------------------*/
