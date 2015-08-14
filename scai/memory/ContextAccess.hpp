/**
 * @file ContextAccess.hpp
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
 * @brief Definition of an access class to access a context.
 *
 * @author Thomas Brandes
 * @date 14.07.2011
 * @since 1.0.0
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/memory/Access.hpp>

// others
#include <scai/memory/Context.hpp>

// logging
#include <scai/logging.hpp>

namespace memory
{

/** @brief This class is used for accessing a Context to enable computations on it.
 *
 *  This class guarantees that enable of a context is always followed by
 *  a disable of the context in case of exception.
 *
 *  Member variables contain source code location where access to the context
 *  was enabled.
 */
class COMMON_DLL_IMPORTEXPORT ContextAccess: public Access
{
public:

    /** Constructor for an access, enables the context. */

    ContextAccess( ContextPtr context, const char* filename, int line );

    /** Destructor for an access, disables the context. */

    ~ContextAccess();

    /** Mandatory routine for each access class. */

    virtual void release();

    virtual void writeAt( std::ostream& stream ) const;

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    ContextAccess    (); // do not allow the default constructor

    const Context& mContext;//!<  context that is accessed.

    bool mReleased;//!< true, if access has already been released

    const char* mFile;//!<  file name of access

    int mLine;//!< line number of access
};

}

/** Macro that should be used for any context access
 *
 *  \code
 *  {
 *      SCAI_CONTEXT_ACCESS( cudaContext )
 *      ...   // computations on CUDA device
 *      // context is released at end of the scope by destructor of access object.
 *  }
 *  \endcode
 */

#define SCAI_CONTEXT_ACCESS( context ) memory::ContextAccess LAMA_Ctx_Access__( context, __FILE__, __LINE__ );
