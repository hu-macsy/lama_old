/**
 * @file LAMAInterfaceRegistry.hpp
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
 * @brief LAMAInterfaceRegistry.hpp
 * @author Thomas Brandes
 * @date 12.04.2013
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <common/config.hpp>

// base classes
#include <common/NonCopyable.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <memory/Context.hpp>

// macros
#include <lama/macros/unique_name.hpp>

#include <map>

namespace lama
{

/** @brief This class provides a registry that gives for a context its corresponding
 *         LAMAInterface.
 *
 *  At program start the interfaces are dynamically filled with function pointers to
 *  implemented routines that are used later for operations on vectors and matrices.
 */

class COMMON_DLL_IMPORTEXPORT LAMAInterfaceRegistry: common::NonCopyable
{

public:

    virtual ~LAMAInterfaceRegistry();

    /** @brief Get a refernce to the registry. */

    static LAMAInterfaceRegistry& getRegistry();

    /** @brief Get const reference to an interface for a certain context type. */

    const LAMAInterface* getInterface( const memory::ContextType location ) const;

    /** @brief Get a modify reference for a LAMAInterface; if not available an
     *         new interface is generated.
     *
     *  @param[in] location context type for which interface is wanted
     *  @return    reference to the corresponding interface
     *
     *  If an interface is not available a default one will be created
     *  (all function pointers are set to NULL).
     */
    LAMAInterface& modifyInterface( const memory::ContextType location );

    /** @brief Query whether an interface for a certain context is available  */

    bool hasInterface( const memory::ContextType location ) const;

private:

    LAMAInterfaceRegistry();

    static LAMAInterfaceRegistry* instance;

    typedef std::map<memory::ContextType,LAMAInterface*> InterfaceMapType;

    InterfaceMapType mInterfaceMap;

    class CGuard
    {
    public:
        CGuard();
        ~CGuard();
    };

    friend class CGuard;
};

} // namespace

