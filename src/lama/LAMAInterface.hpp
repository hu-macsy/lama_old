/**
 * @file LAMAInterface.hpp
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
 * @brief Interface class for context dependent operations to be implemented.
 * @author Thomas Brandes
 * @date 27.04.2011 (revised 02.04.2013)
 * $Id$
 */
#ifndef LAMA_LAMA_INTERFACE_HPP_
#define LAMA_LAMA_INTERFACE_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Printable.hpp>

// interface structures used in LAMAInterface

#include <lama/BLASInterface.hpp>
#include <lama/UtilsInterface.hpp>

namespace lama
{

/**
 * @brief LAMAInterface is a class that provides access to all device dependent
 *        operations used within LAMA.
 *
 * This class specifies all routines that have to be implemented on
 * a new device.
 *
 * Note: Instead of using virtual routines this class uses function pointer variables.
 *       They are stored in one- or two-dimensional arrays indexed by types.
 */
class LAMA_DLL_IMPORTEXPORT LAMAInterface : public Printable
{
public:

    LAMAInterface();

    virtual ~LAMAInterface();

    /** This method writes the name of the interface into the output stream.
     *  Output message helps to identify which interfaces have been used.
     *
     * @see Printable for more details.
     */
    virtual void writeAt( std::ostream& stream ) const;

    // one member variable for each Interface

    BLASInterface       BLAS;        //!< interface table for BLAS routines

    CSRUtilsInterface   CSRUtils;
    DenseUtilsInterface DenseUtils;
    ELLUtilsInterface   ELLUtils;
    JDSUtilsInterface   JDSUtils;
    DIAUtilsInterface   DIAUtils;
    COOUtilsInterface   COOUtils;
    UtilsInterface      Utils;

protected:

    LAMA_LOG_DECL_STATIC_LOGGER(logger)
};

} //namespace lama

/** macros for getting function pointers
 *
 *  @param function variable and name of the function that is searched
 *  @param loc is context where the routine is needed
 *  @param module where to search, e.g. BLAS, CSRUtils, ...
 *  @param structname: structure of the module, e.g. BLAS1, Conversions
 */

#define LAMA_INTERFACE_FN( function, loc, module, structname )                                \
    typename module##Interface::structname::function function =                               \
       loc->getInterface().module.function();                                                 \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_THROWEXCEPTION( "Method " #module "::" #function " not available on " << *loc ); \
    }

/** same as LAMA_INTERFACE_FN but takes routine on Host if not available at loc */

#define LAMA_INTERFACE_FN_DEFAULT( function, loc, module, structname )                        \
    typename module##Interface::structname::function function =                               \
       loc->getInterface().module.function();                                                 \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_UNSUPPORTED( "Method " #module "::" #function " not available on " << *loc );    \
        loc = ContextFactory::getContext( Context::Host );                                    \
        function = loc->getInterface().module.function();                                     \
        if ( function == NULL )                                                               \
        {                                                                                     \
            LAMA_THROWEXCEPTION( "Method " #module "::" #function                             \
                                 " also not available on " << *loc );                         \
        }                                                                                     \
    }

#define LAMA_INTERFACE_FN_T( function, loc, module, structname, T )                           \
    typename module##Interface::structname<T>::function function;                             \
    function = loc->getInterface().module.function<T>();                                      \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_THROWEXCEPTION( "Method " #module "::" #function " not available on " << *loc ); \
    }

#define LAMA_INTERFACE_FN_t( function, loc, module, structname, T )                           \
    module##Interface::structname<T>::function function;                                      \
    function = loc->getInterface().module.function<T>();                                      \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_THROWEXCEPTION( "Method " #module "::" #function " not available on " << *loc ); \
    }

#define LAMA_INTERFACE_FN_DEFAULT_T( function, loc, module, structname, T )                   \
    typename module##Interface::structname<T>::function function =                            \
       loc->getInterface().module.function<T>();                                              \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_UNSUPPORTED( "Method " #module "::" #function " not available on " << *loc );    \
        loc = ContextFactory::getContext( Context::Host );                                    \
        function = loc->getInterface().module.function<T>();                                  \
        if ( function == NULL )                                                               \
        {                                                                                     \
            LAMA_THROWEXCEPTION( "Method " #module "::" #function                             \
                                 " also not available on " << *loc );                         \
        }                                                                                     \
    }

#define LAMA_INTERFACE_FN_TT( function, loc, module, structname, T1, T2 )                     \
    typename module##Interface::structname<T1,T2>::function function;                         \
    function = loc->getInterface().module.function<T1,T2>();                                  \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_THROWEXCEPTION( "Method " #module "::" #function " not available on " << *loc ); \
    }

#endif // LAMA_LAMA_INTERFACE_HPP_
