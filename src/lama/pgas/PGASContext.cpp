/**
 * @file PGASContext.cpp
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
 * @brief Host context where memory is pinned for fast transfer using PGAS
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */

#include <lama/pgas/PGASContext.hpp>
#include <lama/pgas/PGASInterface.hpp>

#include <lama/ContextAccess.hpp>
#include <lama/TaskSyncToken.hpp>
#include <lama/tracing.hpp>

#include <lama/exception/LAMAAssert.hpp>

#include <boost/bind.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( PGASContext::logger, "Context.PGASContext" )

PGASContext::PGASContext()
{
    LAMA_LOG_DEBUG( logger, "PGASContext created, allows communication via PGAS" )
}

PGASContext::~PGASContext()
{
    LAMA_LOG_DEBUG( logger, "~PGASContext()" )
}

void PGASContext::writeAt( std::ostream& stream ) const
{
    // write identification of this object
    stream << "PGASContext";
}

void* PGASContext::allocate( const size_t size ) const
{
    LAMA_REGION( "PGASContext::allocate" )
    LAMA_LOG_DEBUG( logger, *this << ": allocated " << size << " bytes" )
    PGASInterface::getInstance()->syncronizeAll();
    return PGASInterface::getInstance()->allocate( size );
}

void PGASContext::allocate( ContextData& contextData, const size_t size ) const
{
//    size_t maxsize= PGASInterface::getInstance()->maxToAll(size);
    contextData.pointer = allocate( size );
}

void PGASContext::free( void* pointer, const size_t size ) const
{
    LAMA_REGION( "PGASContext::free" )
    PGASInterface::getInstance()->free( pointer, size );
    LAMA_LOG_DEBUG( logger, *this << ": freed " << size << " bytes, pointer = " << pointer )
}

void PGASContext::free( ContextData& contextData ) const
{
    LAMA_ASSERT_EQUAL_ERROR( contextData.context->getType(), getType() )
    free( contextData.pointer, contextData.size );
}

void PGASContext::memcpy( void* dst, const void* src, const size_t size ) const
{
    ::memcpy( dst, src, size );
}

std::auto_ptr<SyncToken> PGASContext::memcpyAsync( void* dst, const void* src, const size_t size ) const
{
    return std::auto_ptr<SyncToken>( new TaskSyncToken( boost::bind( &::memcpy, dst, src, size ) ) );
}

bool PGASContext::cancpy( const ContextData& dst, const ContextData& src ) const
{
    return dst.context->getType() == getType() && src.context->getType() == getType();
}

void PGASContext::memcpy( ContextData& dst, const ContextData& src, const size_t size ) const
{
    LAMA_ASSERT_ERROR( dst.context->getType() == getType() && src.context->getType() == getType(),
                       "Can not copy from "<< *(src.context) << " to " << *(dst.context) )
    memcpy( dst.pointer, src.pointer, size );
}

std::auto_ptr<SyncToken> PGASContext::memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const
{
    LAMA_ASSERT_ERROR( dst.context->getType() == getType() && src.context->getType() == getType(),
                       "Can not copy from "<< *(src.context) << " to " << *(dst.context) )
    return memcpyAsync( dst.pointer, src.pointer, size );
}

HostContext::HostContextType PGASContext::getHostType() const
{
    return PGAS;
}

} // namespace lama

