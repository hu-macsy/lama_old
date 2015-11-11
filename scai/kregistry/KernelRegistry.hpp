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
 * @brief Registry class for kernel methods implemented on a certain context.
 * @author Thomas Brandes
 * @date 14.10.2015
 */
#pragma once

#include <scai/kregistry/ContextFunction.hpp>
#include <scai/kregistry/exception/KernelRegistryException.hpp>

#include <scai/logging.hpp>

#include <map>
#include <string>
#include <typeinfo>
#include <cstdlib>
#include <cstring>

#include <iostream>

namespace scai
{

/** Namespace for the Kernel Registry and related data structures. */

namespace kregistry
{

/** Defintion of key type used in the map of the kernel registry. 
 *
 *  An entry is identified by the type of the funciton (returned by type_info)
 *  and by its name.
 */

typedef std::pair<const std::type_info&, const char*> KernelRegistryKey;

/** Output operator for kregistry inteface key */

inline std::ostream& operator<<( std::ostream& stream, const KernelRegistryKey& key )
{
    stream << "KernelRegistryKey( id = " << key.second << ", type = " << key.first.name() << " )";
    return stream;
}

/**
 *  @brief Static class where all kernel routines are registered and can be accessed. 
 */

class KernelRegistry
{
private:

    /** Method that registers a function pointer for a given context in the registry.  
     *
     *  @param[in] key is key ( pair of name and function type )
     *  @param[in] ctx is the context for which fn is implemented
     *  @param[in] fn is the pointer to the function that is registered
     *  @param[in] replace if TRUE a registered entry might be overwritten
     */

    static void registerContextFunction( const KernelRegistryKey& key, common::context::ContextType ctx, VoidFunction fn, bool replace );

    /** Method that unregisters a function pointer for a given context in the registry.  
     *
     *  @param[in] key is key ( pair of name and function type )
     *  @param[in] ctx is the context for which fn has been registered
     *  @param[in] fn is the pointer to the function that has been registered
     */

    static void unregisterContextFunction( const KernelRegistryKey& key, common::context::ContextType ctx, VoidFunction fn );

public:

    /** Enumeration type for possible registration flags.
     *
     *  KERNEL_ADD should be used for kernels with lower priority and KERNEL_REPLACE
     *  for kernels with higher priority. This might be used to prefer optimized kernels
     *  when they are available.
     */

    enum KernelRegistryFlag
    {
        KERNEL_ERASE,    //!< Kernel routine will be removed from registry
        KERNEL_ADD,      //!< Kernel routine will be added if no other one is registered
        KERNEL_REPLACE   //!< Kernel routine will be added and will replace an existing one
    };

    /** This method registers a kernel routine for a given context by its name and signature.  
     *
     *  @tparam    FunctionType specifies the signature of the function
     *  @param[in] fn is the function to register
     *  @param[in] name is the name of the function
     *  @param[in] ctx is the context type for which the function is implemented
     *  @param[in] flag specifies kind of registration
     *
     */

    template<typename FunctionType>
    static void set( FunctionType fn, const char* name, const common::context::ContextType ctx, const KernelRegistryFlag flag )
    {
        KernelRegistryKey key( typeid( FunctionType ), name );
        if ( flag == KERNEL_ERASE )
        {
            unregisterContextFunction( key, ctx, ( VoidFunction ) fn );
        }
        else
        {
            bool replace = ( flag == KERNEL_REPLACE );
            registerContextFunction( key, ctx, ( VoidFunction ) fn, replace );
        }
    }

    /** This  method registers a kernel routine for a given context by its kernel trait.
     *
     *  @tparam KernelTrait struct that must contain FuncType type definition and getId method.
     *
     *  @param[in] fn is the function to register
     *  @param[in] ctx is the context type for which the function is implemented
     *  @param[in] flag specifies kind of registration
     *
     *  For the registration, the trait avoid mainly misspelling of names.
     */

    template<typename KernelTrait>
    static void set( typename KernelTrait::FuncType fn, const common::context::ContextType ctx, const KernelRegistryFlag flag )
    {
        set( fn, KernelTrait::getId(), ctx, flag );
    }

    template<typename FunctionType>
    static void get( FunctionType& fn, const char* name, common::context::ContextType ctx )
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

    /** Get for one specific routines the function pointers for all supported contexts. 
     *
     *  @tparam    FunctionType is the function type of queried kernel routine
     *  @param[out] contextFunction function pointer for each context
     *  @param[in] name is the name of the kernel routine
     */

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

    static void printAll();

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
