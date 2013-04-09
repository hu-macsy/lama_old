/**
 * @file CUDAContextManager.hpp
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
 * @brief Definition of manager class to manage CUDA contexts.
 * @author Thomas Brandes, Jiri Kraus
 * @date 15.07.2011
 * $Id$
 */
#ifndef LAMA_CUDA_CONTEXT_MANAGER_HPP_
#define LAMA_CUDA_CONTEXT_MANAGER_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/ContextManager.hpp>

// others
#include <lama/Context.hpp>

// boost
#include <boost/weak_ptr.hpp>

#include <string>

namespace lama
{

class CUDAContext;

/** maximal number of CUDA devices that might be available */

#define LAMA_MAX_CUDA_DEVICES 4

class LAMA_DLL_IMPORTEXPORT CUDAContextManager: public ContextManager
{
public:

    ContextPtr getContext( int deviceNr )
    {
        return getInstance( deviceNr );
    }

    ~CUDAContextManager();

    /** Static version of the member method getContext. */

    static ContextPtr getInstance( int deviceNr );

private:

    /** Make constructor private, only single static object will be created. */

    CUDAContextManager();

    /** The manager keeps weak pointers to CUDA contexts so it can return
     *  a shared pointer if the contexts are still in use.
     *  Entry is available for each possible device.
     */

    static boost::weak_ptr<CUDAContext> mCUDAContext[LAMA_MAX_CUDA_DEVICES];

    static CUDAContextManager theInstance; //!< singleton instance of this class
};

}

#endif // LAMA_CUDA_CONTEXT_MANAGER_HPP_
