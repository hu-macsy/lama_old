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
#include <scai/common/exception/Exception.hpp>

// std
#include <map>

namespace scai
{

using common::ContextType;

namespace hmemo
{

/* -----------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( Context::logger, "Context" )

Context::Context( ContextType type )
                : mContextType( type ), mEnabled( false ), mFile( NULL ), mLine( 0 )
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

MemoryPtr Context::getHostMemoryPtr() const
{
    // take the host memory of the memory factory

    ContextPtr hostContextPtr = Context::getContextPtr( common::context::Host );
    return hostContextPtr->getMemoryPtr();
}

/* ---------------------------------------------------------------------------------*/

ContextPtr Context::getContextPtr( ContextType type, int deviceNr )
{
    return create( type, deviceNr );
}

} /* end namespace hmemo */

} /* end namespace scai */

/* ---------------------------------------------------------------------------------*/
