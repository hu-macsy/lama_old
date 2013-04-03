/**
 * @file PGASContext.hpp
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
#ifndef LAMA_PGAS_HOST_CONTEXT_HPP_
#define LAMA_PGAS_HOST_CONTEXT_HPP_

#include <boost/weak_ptr.hpp>

#include <lama/HostContext.hpp>

namespace lama
{

/** Alternative context to DefaultHostContext so that memory will be allocated
 *  symmetric on all nodes.
 *
 *  As Openshmem needs symmetric memory to communicate with each other this
 *  context allows communication using Openshmem.
 */

class PGASContext: public lama::HostContext
{
    friend class PGASContextManager; // can only create context

public:

    virtual ~PGASContext();

    virtual void* allocate( const size_t size ) const;

    virtual void allocate( ContextData& contextData, const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void free( ContextData& contextData ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual std::auto_ptr<SyncToken> memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual bool cancpy( const ContextData& dst, const ContextData& src ) const;

    virtual void memcpy( ContextData& dst, const ContextData& src, const size_t size ) const;

    virtual std::auto_ptr<SyncToken> memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const;

    virtual HostContextType getHostType() const;

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    PGASContext( );

    virtual void writeAt( std::ostream& stream ) const;
};

}

#endif // LAMA_PGAS_HOST_CONTEXT_HPP_
