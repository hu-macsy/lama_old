/**
 * @file ContextFactory.hpp
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
 * @brief ContextFactory.hpp
 * @author Thomas Brandes
 * @date 10.07.2011
 * $Id$
 */
#ifndef LAMA_CONTEXT_FACTORY_HPP_
#define LAMA_CONTEXT_FACTORY_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/NonCopyable.hpp>

// others
#include <lama/Context.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

namespace lama
{

class ContextManager;

/**
 *   The value here is used as default if no device number is specified.
 *   It should be a value that is not used by any context manager for
 *   a real device.
 */
#define LAMA_DEFAULT_DEVICE_NUMBER -1

/** Singleton class that provides (shared) pointers to a context.
 *  of a given type.
 *
 *  \code
 *  ContextPtr hostContext  = ContextFactory::getContext( Context::Host );
 *  ContextPtr cudaContext  = ContextFactory::getContext( Context::CUDA );
 *  ContextPtr cudaContext1 = ContextFactory::getContext( Context::CUDA, 1 );
 *  \endcode
 *
 *  Each context is returned via a shared pointer so the context is freed
 *  automatically if no more reference is kept.
 *
 *  For each supported context type a context manager must register to the
 *  factory otherwise the type remains unsupported.
 */
class LAMA_DLL_IMPORTEXPORT ContextFactory: NonCopyable
{
public:

    /** This method can be used to get access to a context so that data
     *  can be allocated on it.
     *
     *  The shared pointer guarantees that a context will live at least
     *  as long as arrays are allocated on it.
     *
     *  @param[in] type     is the type of context that is wanted
     *  @param[in] deviceNr is used for multiple devices of the same type (default is LAMA_DEFAULT_DEVICE_NUMBER)
     *
     *  @return             a context of the requested type, if available.
     *
     *  @throws Exception if the context of the requested type is not available
     */
    static ContextPtr getContext( const Context::ContextType type, int deviceNr = LAMA_DEFAULT_DEVICE_NUMBER );

    /** TODO[doxy] Complete Description.
     *
     *  @param[in] type     is the type of context that is wanted
     *  @return             TODO[doxy] Complete Description.
     */
    ContextManager* getContextManager( const Context::ContextType type );

    /** Checks if a context of the passed type is available.
     *
     * @param[in] type  is the type of context that is wanted
     * @return          if a context of the passed type is available
     */
    static bool hasContext( const Context::ContextType type );

    /** Each context type that should be supported must add its context manager so that the
     *  the factory can forward a query for a context type to the corresponding manager.
     *
     *  @param[in] type             is the type of context for which the manager is responsibile
     *  @param[in] contextManager   reference to the manager for all contexts of a given type.
     */
    void addContextManager( const Context::ContextType type, ContextManager& contextManager );

    /** Method that returns a reference to the context factory. */

    static ContextFactory& getFactory();

    /** Release all context managers. */

    static void release();

private:

    ContextFactory();

    virtual ~ContextFactory();

    /** Each possible context will have its own context manager. */

    ContextManager* mContextManager[Context::MaxContext];

    static const char* theContextIds[Context::MaxContext];

    static ContextFactory* theContextFactory;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

}

#endif // LAMA_CONTEXT_FACTORY_HPP_
