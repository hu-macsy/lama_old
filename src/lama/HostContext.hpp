/**
 * @file HostContext.hpp
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
 * @brief Definition of a class that defines a context for HOST computation/allocation.
 * @author Thomas Brandes
 * @date 14.07.2011
 * @since 1.0.0
 */
#ifndef LAMA_HOSTCONTEXT_HPP_
#define LAMA_HOSTCONTEXT_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Context.hpp>
#include <lama/task/TaskSyncToken.hpp>

namespace lama
{

/** @brief This class implements the base class HOST context.
 *
 *  The class HostContext remains abstract as it does not provide
 *  routines to allocate/free memory.
 */

class LAMA_DLL_IMPORTEXPORT HostContext: public Context
{

public:

    /** Due to more efficient transfers between Host and other contexts there might be
     *  different HOST contexts that will be distinguished by their type.
     */

    enum HostContextType
    {
        DefaultHost, CUDAHost, OpenCLHost, PGAS, MaxHostContext
    };

    virtual ~HostContext();

    virtual bool canUseData( const Context& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    virtual HostContextType getHostType() const = 0;

    /** Providing implementation for Context::getSyncToken with covariant return type */

    virtual TaskSyncToken* getSyncToken() const;

protected:

    HostContext();
};

} // namespace LAMA

#endif // LAMA_HOSTCONTEXT_HPP_
