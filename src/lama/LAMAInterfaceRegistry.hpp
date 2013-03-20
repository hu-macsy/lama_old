/**
 * @file LAMAInterfaceRegistry.hpp
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
 * @brief LAMAInterfaceRegistry.hpp
 * @author brandes
 * @date 28.04.2011
 * $Id$
 */

#ifndef LAMA_LAMAINTERREGISTRY_HPP_
#define LAMA_LAMAINTERREGISTRY_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/Context.hpp>

// macros
#include <lama/macros/unique_name.hpp>

#include <map>

namespace lama
{

class LAMA_DLL_IMPORTEXPORT LAMAInterfaceRegistry
{
public:

    virtual ~LAMAInterfaceRegistry();

    static LAMAInterfaceRegistry& getRegistry();

    void addInterface( const ContextType location, LAMAInterface* lamaInterface );

    const LAMAInterface* getInterface( const ContextType location ) const;

    bool hasInterface( const ContextType location ) const;

private:
    LAMAInterfaceRegistry();
    LAMAInterfaceRegistry( const LAMAInterfaceRegistry& other );

    LAMAInterfaceRegistry& operator=( const LAMAInterfaceRegistry& other );

    static LAMAInterfaceRegistry* instance;

    typedef std::map<ContextType,LAMAInterface*> InterfaceMapType;

    InterfaceMapType mInterfaceMap;

    class CGuard
    {
    public:
        CGuard();
        ~CGuard();
    };
    friend class CGuard;
};

template<typename T>
class LAMAInterfaceRegistration
{
public:

    LAMAInterfaceRegistration( const ContextType location )
    {
        LAMAInterfaceRegistry& reg = LAMAInterfaceRegistry::getRegistry();
        reg.addInterface( location, new T() );
    }
};

#define LAMA_LAMAINTERFACE_REGISTRATION( location, type ) \
    static lama::LAMAInterfaceRegistration<type >       \
    LAMA_UNIQUE_NAME( lamaInterfaceRegObj, type )( location )
}

#endif // LAMA_LAMAINTERREGISTRY_HPP_
