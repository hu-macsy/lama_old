/**
 * @file ContextManager.hpp
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
 * @brief ContextManager.hpp
 * @author Thomas Brandes
 * @date 10.07.2011
 * @since 1.0.0
 */
#ifndef LAMA_CONTEXT_MANAGER_HPP_
#define LAMA_CONTEXT_MANAGER_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <common/NonCopyable.hpp>
#include <common/Printable.hpp>

// others
#include <lama/Context.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

namespace lama
{

/** @brief Base class for manager of context objects.
 *
 *  Each derived class of Context should have a derived class counterpart
 *  of derived class of ContextManager.
 *
 *  A context manager registers via the context type to the factory so
 *  that the factory can provide a context for the corresponding type.
 *
 *  The concept of manager and factory provides a mechanism for extendability
 *  with new context classes and types.
 */

class COMMON_DLL_IMPORTEXPORT ContextManager: public Printable, private common::NonCopyable
{
public:

    /** Destructor. */

    virtual ~ContextManager();

    /** This routine must be implemented by each context manager individually.
     *
     */
    virtual ContextPtr getContext( int deviceNr ) = 0;

    /** Override Printable::writeAt. */

    virtual void writeAt( std::ostream& stream ) const;

protected:

    /** Constructor for a context manager.
     *
     *  param[in] type context type for which the context manager will be responsible
     *
     *  The constructor adds itself to the context factory.
     */
    ContextManager( ContextType type );

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    /** This method sets the context manager as responsible for its context type
     *  in the factory.
     */
void    registerFactory();

private:

    ContextType mContextType; //!< context type for which manager is responsible
};

}

#endif // LAMA_CONTEXT_MANAGER_HPP_
