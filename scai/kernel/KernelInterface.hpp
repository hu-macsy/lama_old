/**
 * @file KernelInterface.hpp
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
 * @brief Interface class for context dependent operations to be implemented.
 * @author Thomas Brandes
 * @date 14.10.2015
 */
#pragma once

#include <scai/kernel/ContextFunction.hpp>

#include <scai/logging.hpp>
#include <scai/common/exception/Exception.hpp>

#include <map>
#include <string>
#include <typeinfo>
#include <cstdlib>
#include <cstring>

#include <iostream>

namespace scai
{

namespace interface
{

typedef std::pair<const std::type_info&, const char*> InterfaceKey;

// Output operator for kernel inteface key

static std::ostream& operator<<( std::ostream& stream, const InterfaceKey& key )
{
    stream << "InterfaceKey( id = " << key.second << ", type = " << key.first.name() << " )";
    return stream;
}

/* --------------------------------------------------------------------------- *
 *  KernelInterface ( static class )                                           *
 * --------------------------------------------------------------------------- */

class KernelInterface
{
private:

    /** Type for unique key in registration */

    static void registerContextFunction( const InterfaceKey& key, ContextType ctx, VoidFunction fn )
    {
        SCAI_LOG_INFO( logger, "register ctx = " << ctx << " with " << key )

        InterfaceMap::iterator it = theInterfaceMap.find( key );

        if ( it == theInterfaceMap.end() )
        {
            SCAI_LOG_DEBUG( logger, "register: no entry yet, will add it" )

            _ContextFunction routine;

            routine.set( ctx, fn );

            theInterfaceMap.insert( std::pair<InterfaceKey, _ContextFunction>( key, routine ) );

            SCAI_LOG_DEBUG( logger, "added" ) 
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "register: entry available, set it for ctx = " << ctx )

            it->second.set( ctx, fn );
        }
    }

    static void getContextFunction( _ContextFunction& contextFunction, const InterfaceKey& key )
    {
        InterfaceMap::const_iterator it = theInterfaceMap.find( key );

        if ( it != theInterfaceMap.end() )
        {
            contextFunction = it->second;
        }
        else
        {
            SCAI_LOG_ERROR( logger, "no context function with name = " << key.second 
                               << ", type = " << key.first.name() << " available" )

            // maybe we should throw an exception

            contextFunction.clear();
        }
    }

public:

    template<typename KernelTrait>
    static void set( typename KernelTrait::FuncType fn, ContextType ctx )
    {
        InterfaceKey key( typeid( typename KernelTrait::FuncType ), KernelTrait::getId() );
        registerContextFunction( key, ctx, ( VoidFunction ) fn );
    }

    template<typename FunctionType>
    static void set( FunctionType fn, const char* name, ContextType ctx )
    {
        InterfaceKey key( typeid( FunctionType ), name );
        registerContextFunction( key, ctx, ( VoidFunction ) fn );
    }

    template<typename FunctionType>
    static void get( FunctionType& fn, const char* name, ContextType ctx )
    {
        SCAI_LOG_INFO( logger, "get function pointer for kernel routine " << name 
                                << ", func type = " << typeid( FunctionType ).name() 
                                << ", context = " << ctx  )

        InterfaceKey key( typeid( FunctionType ), name );

        typename InterfaceMap::const_iterator it = theInterfaceMap.find( key );

        if ( it != theInterfaceMap.end() )
        {
            SCAI_LOG_DEBUG( logger, "function registered" )

            const _ContextFunction& routine = it->second;
            fn = ( FunctionType ) routine.get( ctx );   // cast required

            SCAI_LOG_DEBUG( logger, "function for context = " << fn )
        }
        else
        {
            SCAI_LOG_INFO( logger, "function never registered." )
        }
    }

    /** Get all function pointers for a certain kernel routine */

    template<typename FunctionType>
    static void get( ContextFunction<FunctionType>& contextFunction, const char* name )
    {
        InterfaceKey key( typeid( FunctionType ), name );

        SCAI_LOG_INFO( logger, "get all function pointers for kernel routine by " << key )

        typename InterfaceMap::const_iterator it = theInterfaceMap.find( key );

        if ( it != theInterfaceMap.end() )
        {
            SCAI_LOG_DEBUG( logger, "entry found in interface" )

            contextFunction.assign( it->second );
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "entry not found in interface, set all NULL" )

            contextFunction.clear();   // just for safety

            COMMON_THROWEXCEPTION( "No context function registered, name = " << name 
                                     << ", type = " << typeid( FunctionType ).name() )
        }
    }

    /** Help routine that prints all registered kernel routines */

    static void printAll()
    {
        InterfaceMap::const_iterator it;

        std::cout << "KernelInterface:" << std::endl;
        std::cout << "================" << std::endl;

        for ( it = theInterfaceMap.begin(); it != theInterfaceMap.end(); ++it )
        {
            std::cout << "Entry: key = " << it->first;
            std::cout << ", ctx = " << it->second.printIt();
            std::cout << std::endl; 
        }

        std::cout << "================" << std::endl;
    }

protected:

    KernelInterface();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    class Compare
    {
    public:

        // return x > y

        bool operator()( const InterfaceKey& x, const InterfaceKey& y );

    };

    typedef std::map<InterfaceKey, _ContextFunction, Compare> InterfaceMap;

    static InterfaceMap theInterfaceMap;
};

} /* end namespace interface */

} /* end namespace scai */
