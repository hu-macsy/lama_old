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

#include "ContextType.hpp"

#include <map>
#include <string>
#include <typeinfo>
#include <cstdlib>

#include <iostream>

namespace scai
{

namespace interface
{

/** Type definition of un untyped function pointer */

typedef void ( *VoidFunction )();

/* --------------------------------------------------------------------------- *
 * common base class for ContextFunction                                         *
 * --------------------------------------------------------------------------- */

/** Common base class for template class ContextFunction 
 *
 *  An object of _ContextFunction contains a function pointer for each context
 *  where the function pointer might be NULL for unsupported context
 */

class _ContextFunction
{
public:

    /** Default constructor initializes all function pointers with NULL */

    _ContextFunction();

    _ContextFunction( const _ContextFunction& other );

    void clear();   // sets all function pointers to NULL

    inline VoidFunction get( ContextType ctx ) const
    {
        return mContextFuncArray[ ctx ];
    }

    inline void set( ContextType ctx, VoidFunction fn )
    {
        mContextFuncArray[ ctx ] = fn;
    }

    ContextType validContext( ContextType preferedCtx );

    ContextType validContext( const _ContextFunction& other, ContextType preferedCtx );

protected:

    // array with function pointer for each context

    VoidFunction mContextFuncArray[context::MaxContext];
};

/* --------------------------------------------------------------------------- *
 * template class for ContextFunction                                            *
 * --------------------------------------------------------------------------- */

template<typename FunctionType> 
class ContextFunction : public _ContextFunction
{
public:

    /** Constructor with name of the routine, searches for registered functions. */

    ContextFunction( const std::string& name );

    /** Override default copy constructor */

    ContextFunction( const ContextFunction& other ) : _ContextFunction( other )
    {
    }

    FunctionType get( ContextType ctx ) const
    {
        return ( FunctionType ) mContextFuncArray[ ctx ];
    }

    void set( ContextType ctx, FunctionType fn )
    {
        mContextFuncArray[ ctx ] = ( VoidFunction ) fn;
    }

    FunctionType operator() ( ContextType ctx )
    {
        if ( mContextFuncArray[ ctx ] == NULL )
        {
            // Throw exception
            std::cout << "kernel routine not available context = " << ctx << std::endl;
            std::cout << "STOP" << std::endl;
            exit( - 1 );
        }

        return ( FunctionType ) mContextFuncArray[ ctx ];
    }

    using _ContextFunction::validContext;
};

/**
 * Template class for ContextFunction by using a Trait 
 *
 * @tname <ContextFunctionTrait>  struct that constains signature and name of the context function.
 *
 * \begincode
 *     // Example of ContextFunctionTrait
 *     struct isSorted
 *     {
 *         bool ( *FuncType ) ( const double* array, int n, bool ascending );
 *         static inline const char* getId() { return "isSorted" };
 *     };
 * \endcode
 */

template<typename ContextFunctionTrait> 
class ContextFunctionByTrait : public ContextFunction<typename ContextFunctionTrait::FuncType>
{
public:

    typedef typename ContextFunctionTrait::FuncType ContextFunctionType;

    ContextFunctionByTrait();
};

/* --------------------------------------------------------------------------- *
 *  KernelInterface ( static class )                                           *
 * --------------------------------------------------------------------------- */

class KernelInterface
{

private:

    // Interface as map requires one common function type

    typedef std::pair<const std::type_info&, const std::string> InterfaceKey;

    static void registerContextFunction( const InterfaceKey& key, ContextType ctx, VoidFunction fn )
    {
        std::cout << "registerContextFunction, name = " << key.first.name() << ", " << key.second  << std::endl;

        typename InterfaceMap::iterator it = theInterfaceMap.find( key );

        if ( it == theInterfaceMap.end() )
        {
            std::cout << "registerContextFunction: new routine" << std::endl;

            _ContextFunction routine;

            routine.set( ctx, fn );

            theInterfaceMap.insert( std::make_pair( key, routine ) );

            std::cout << "added" << std::endl;
        }
        else
        {
            std::cout << "registerContextFunction: routine with next context" << std::endl;

            it->second.set( ctx, fn );
        }
    }

    static void getContextFunction( _ContextFunction& contextFunction, const InterfaceKey& key )
    {
        typename InterfaceMap::const_iterator it = theInterfaceMap.find( key );

        if ( it != theInterfaceMap.end() )
        {
            contextFunction = it->second;
        }
        else
        {
            // maybe we should throw an exception

            contextFunction.clear();
        }
    }

public:

    template<typename ContextFunctionTrait>
    static void set( typename ContextFunctionTrait::FuncType fn, ContextType ctx )
    {
        InterfaceKey key( typeid( typename ContextFunctionTrait::FuncType ), ContextFunctionTrait::getId() );
        registerContextFunction( key, ctx, ( VoidFunction ) fn );
    }

    template<typename FunctionType>
    static void set( FunctionType fn, const std::string& name, ContextType ctx )
    {
        InterfaceKey key( typeid( FunctionType ), name );
        registerContextFunction( key, ctx, ( VoidFunction ) fn );
    }

    template<typename FunctionType>
    static void get( FunctionType& fn, const std::string& name, ContextType ctx )
    {
        InterfaceKey key( typeid( FunctionType ), name );

        typename InterfaceMap::const_iterator it = theInterfaceMap.find( key );

        if ( it != theInterfaceMap.end() )
        {
            const _ContextFunction& routine = it->second;
            fn = ( FunctionType ) routine.get( ctx );   // cast required
        }
    }

    /** Get all function pointers for a certain kernel routine */

    template<typename FunctionType>
    static void get( FunctionType mContextFuncArray[], const std::string& name )
    {
        std::cout << "get all function pointers for kernel routine " << name 
                  << ", func type = " << typeid( FunctionType ).name() << std::endl;

        InterfaceKey key( typeid( FunctionType ), name );

        typename InterfaceMap::const_iterator it = theInterfaceMap.find( key );

        if ( it != theInterfaceMap.end() )
        {
            std::cout << "entry found in interface" << std::endl;

            const _ContextFunction& routine = it->second;

            for ( int i = 0; i < context::MaxContext; ++i )
            {
                mContextFuncArray[i] = ( FunctionType ) routine.get( static_cast<ContextType>( i ) );
            }
        }
        else
        {
            std::cout << "entry not found in interface, set all NULL" << std::endl;

            for ( int i = 0; i < context::MaxContext; ++i )
            {
                mContextFuncArray[i] = NULL;
            }

            std::cout << "STOP" << std::endl;
            exit( - 1 );
        }
    }

protected:

    KernelInterface();

private:

    class Compare
    {
    public:
        // return x > y
        bool operator()( const InterfaceKey& x, const InterfaceKey& y )
        {
            const std::string& xstr = x.second;
            const std::string& ystr = y.second;

            if ( xstr > ystr )
            {
                 return true;
            }

            if ( xstr < ystr )
            {
                 return false;
            }
            
            // both have same id, so take typename to distinguish

            return x.first.name() > y.first.name();
        }
    };

    typedef std::map<InterfaceKey, _ContextFunction, Compare> InterfaceMap;

    static InterfaceMap theInterfaceMap;
};

/* --------------------------------------------------------------------------- */

template<typename FunctionType> 
ContextFunction<FunctionType>::ContextFunction( const std::string& name )
{
    std::cout << "ContextFunction<" << typeid( FunctionType ).name() << ">" << std::endl;

    KernelInterface::get( ( FunctionType* ) mContextFuncArray, name );
}

template<typename ContextFunctionTrait> 
ContextFunctionByTrait<ContextFunctionTrait>::ContextFunctionByTrait() 

  : ContextFunction<ContextFunctionType>( ContextFunctionTrait::getId() )

{
    std::cout << "ContextFunctionByTrait<" << typeid( ContextFunctionType ).name()
              << ", id = " << ContextFunctionTrait::getId() << std::endl;
}

} /* end namespace interface */

} /* end namespace scai */
