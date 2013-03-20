/**
 * @file PGASContextManager.hpp
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
 * @brief Manager for PGAS host context.
 * @author Michael Drost
 * @date 01.02.2012
 * $Id$
 */
#ifndef LAMA_CUDA_HOST_CONTEXT_MANAGER_HPP_
#define LAMA_CUDA_HOST_CONTEXT_MANAGER_HPP_

#include <boost/weak_ptr.hpp>

#include <lama/config.hpp>
#include <lama/HostContext.hpp>
#include <lama/ContextManager.hpp>

namespace lama
{

class PGASContext;

/** This singleton class provides the one and only-one manager
 *  that provides access to a Host Context.
 *
 *  The HostManager will add itself to the ContextFactory so that OSHContext
 *  pointers can be obtained by the Context factory.
 */

class LAMA_DLL_IMPORTEXPORT PGASContextManager: public ContextManager
{

public:

    virtual ContextPtr getContext();

    virtual ContextPtr getContext( int );

    ~PGASContextManager();

    static void setAsCurrent();

private:

    PGASContextManager();

    static boost::weak_ptr<const class PGASContext> contextInstance;

    static PGASContextManager theInstance; // singleton instance of this class
};

}

#endif // LAMA_CUDA_HOST_CONTEXT_MANAGER_HPP_
