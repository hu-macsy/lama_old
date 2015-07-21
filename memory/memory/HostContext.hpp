/**
 * @file HostContext.hpp
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
 * @brief Definition of context class that stands for execution on the host, 
 *        i.e. convenctional CPU or multiple CPUs.
 * @author Thomas Brandes
 * @date 10.07.2011
 */

#pragma once

// for dll_import
#include <common/config.hpp>
#include <common/Thread.hpp>

// base classes
#include <memory/Context.hpp>
#include <tasking/TaskSyncToken.hpp>

// common
#include <common/shared_ptr.hpp>

namespace memory
{

/** @brief This class implements the default HOST context.
 *
 *  This class is implemented as a singleton, only one default host
 *  context is available.
 *
 *  The default host context allocates/frees data in the usual way.
 */

class COMMON_DLL_IMPORTEXPORT HostContext: 

    public Context, 
    private Context::Register<HostContext>,
    public common::enable_shared_from_this<HostContext>
{

public:

    virtual ~HostContext();

    /** Static method required for create to use in Context::Register */

    static ContextPtr create( int deviceNr );

    /** Static method required for Context::Register */

    static ContextType createValue();

    /** Override Printable::writeAt with version for this class. */

    virtual void writeAt( std::ostream& stream ) const;

    /** Implement pure method Context::canUseMemory */

    virtual bool canUseMemory( const Memory& memory ) const;

    virtual tasking::TaskSyncToken* getSyncToken() const;

    virtual MemoryPtr getMemoryPtr() const;

private:

    HostContext();

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

inline ContextType HostContext::createValue() 
{
    return context::Host;
}

}  // namespace

