/**
 * @file KernelRegistry.hpp
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
 * @brief Registry class for context dependent operations to be implemented.
 * @author Thomas Brandes
 * @date 14.10.2015
 */
#pragma once

#include <scai/kregistry/ContextFunction.hpp>
#include <scai/kregistry/KernelRegistryException.hpp>

#include <scai/logging.hpp>

#include <map>
#include <string>
#include <typeinfo>
#include <cstdlib>
#include <cstring>

#include <iostream>

namespace scai
{

namespace kregistry
{

typedef std::pair<const std::type_info&, const char*> KernelRegistryKey;

// Output operator for kregistry inteface key

static std::ostream& operator<<( std::ostream& stream, const KernelRegistryKey& key )
{
    stream << "KernelRegistryKey( id = " << key.second << ", type = " << key.first.name() << " )";
    return stream;
}

/* --------------------------------------------------------------------------- *
 *  KernelRegistry ( static class )                                           *
 * --------------------------------------------------------------------------- */

class KernelRegistry
{
private:

    /** Method that registers a function pointer for a given context in the registry.  
     *
     *  @param[in] key is key ( pair of name and function type )
     *  @param[in] ctx is the context where the function works
     *  @param[in] replace if TRUE a registered entry might be overwritten
     */

    static void registerContextFunction( const KernelRegistryKey& key, ContextType ctx, VoidFunction fn, bool replace );

public:

    template<typename KernelTrait>
    static void set( typename KernelTrait::FuncType fn, ContextType ctx, bool replace = false )
    {
        KernelRegistryKey key( typeid( typename KernelTrait::FuncType ), KernelTrait::getId() );
        registerContextFunction( key, ctx, ( VoidFunction ) fn, replace );
    }

    template<typename FunctionType>
    static void set( FunctionType fn, const char* name, ContextType ctx, bool replace = false )
    {
        KernelRegistryKey key( typeid( FunctionType ), name );
        registerContextFunction( key, ctx, ( VoidFunction ) fn, replace );
    }

    template<typename FunctionType>
    static void get( FunctionType& fn, const char* name, ContextType ctx )
    {
        KernelRegistryKey key( typeid( FunctionType ), name );

        SCAI_LOG_INFO( logger, "get function pointer for kregistry routine " 
                                << ", key = " << key 
                                << ", context = " << ctx  )

        typename KernelMap::const_iterator it = theKernelMap.find( key );

        if ( it != theKernelMap.end() )
        {
            const _ContextFunction& routine = it->second;
            fn = ( FunctionType ) routine.get( ctx );   // cast required
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "function never registered." )
        }
    }

    /** Get all function pointers for a certain kregistry routine */

    template<typename FunctionType>
    static void get( ContextFunction<FunctionType>& contextFunction, const char* name )
    {
        KernelRegistryKey key( typeid( FunctionType ), name );

        SCAI_LOG_INFO( logger, "get all function pointers for kernel routine by " << key )

        typename KernelMap::const_iterator it = theKernelMap.find( key );

        if ( it != theKernelMap.end() )
        {
            SCAI_LOG_DEBUG( logger, "entry found in kernel registry" )

            contextFunction.assign( it->second );
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "entry not found in kernel registry, set function pointers for each ctxt to NULL" )

            contextFunction.clear();   // just for safety

            SCAI_THROWEXCEPTION( KernelRegistryException, "No context function registered, " << key )
        }
    }

    /** Help routine that prints all registered kernel routines */

    static void printAll()
    {
        KernelMap::const_iterator it;

        std::cout << "KernelRegistry:" << std::endl;
        std::cout << "================" << std::endl;

        for ( it = theKernelMap.begin(); it != theKernelMap.end(); ++it )
        {
            std::cout << "Entry: key = " << it->first;
            std::cout << ", ctx = " << it->second.printIt();
            std::cout << std::endl; 
        }

        std::cout << "================" << std::endl;
    }

protected:

    KernelRegistry();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    class Compare
    {
    public:

        // return x > y

        bool operator()( const KernelRegistryKey& x, const KernelRegistryKey& y );

    };

    typedef std::map<KernelRegistryKey, _ContextFunction, Compare> KernelMap;

    static KernelMap theKernelMap;
};

} /* end namespace kregistry */

} /* end namespace scai */
