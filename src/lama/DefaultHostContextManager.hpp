/**
 * @file DefaultHostContextManager.hpp
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
 * @brief DefaultHostContextManager.hpp
 * @author Thomas Brandes
 * @date 10.07.2011
 * @since 1.0.0
 */
#ifndef LAMA_DEFAULT_HOST_CONTEXT_MANAGER_HPP_
#define LAMA_DEFAULT_HOST_CONTEXT_MANAGER_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <lama/ContextManager.hpp>

// boost
#include <boost/weak_ptr.hpp>

namespace lama
{

/** This singleton class provides the one and only-one manager
 *  that provides access to a Default Host Context.
 *
 *  The DefaultHostManager will add itself to the ContextFactory so that HostContext
 *  pointers can be obtained by the Context factory.
 */

class COMMON_DLL_IMPORTEXPORT DefaultHostContextManager: public ContextManager
{

public:

    virtual ContextPtr getContext( int deviceNr );

    ~DefaultHostContextManager();

    /** This method sets the DefaultHostContextManager as manager for Context::Host.
     *  Due to singleton, it is a static routine.
     */
    static void setAsCurrent();

    /** Override Printable::writeAt. */

    virtual void writeAt( std::ostream& stream ) const;

private:

    DefaultHostContextManager();

    static boost::weak_ptr<class DefaultHostContext> contextInstance;

    static DefaultHostContextManager theInstance; // singleton instance of this class
};

}

#endif // LAMA_DEFAULT_HOST_CONTEXT_MANAGER_HPP_
