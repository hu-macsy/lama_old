/**
 * @file ContextFunction.hpp
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
 * @brief Definition of class for typed function pointers, one for each context
 * @author Thomas Brandes
 * @date 14.10.2015
 */
#pragma once

#include <scai/common/ContextType.hpp>

#include <map>
#include <string>
#include <typeinfo>
#include <cstdlib>

#include <iostream>

namespace scai
{

namespace kregistry
{

/** Type definition of un untyped function pointer */

typedef void ( *VoidFunction )();

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

    void assign( const _ContextFunction& other );

    inline VoidFunction get( common::ContextType ctx ) const
    {
        return mContextFuncArray[ ctx ];
    }

    inline void set( common::ContextType ctx, VoidFunction fn )
    {
        mContextFuncArray[ ctx ] = fn;
    }

    std::string printIt() const;

    common::ContextType validContext( common::ContextType preferedCtx );

    common::ContextType validContext( const _ContextFunction& other, common::ContextType preferedCtx );

    common::ContextType validContext( const _ContextFunction& other1, 
                                      const _ContextFunction& other2, 
                                      common::ContextType preferedCtx );

protected:

    // array with function pointer for each context

    VoidFunction mContextFuncArray[common::context::MaxContext];
};

/* --------------------------------------------------------------------------- *
 * template class for ContextFunction                                            *
 * --------------------------------------------------------------------------- */

template<typename FunctionType>
class ContextFunction : public _ContextFunction
{
public:

    /** Only default Constructor at this time */

    ContextFunction()
    {
    }

    /** Override default copy constructor */

    ContextFunction( const ContextFunction& other ) : _ContextFunction( other )
    {
    }

    // provide typed get

    FunctionType get( common::ContextType ctx ) const
    {
        return ( FunctionType ) mContextFuncArray[ ctx ];
    }

    // provide typed set

    void set( common::ContextType ctx, FunctionType fn )
    {
        mContextFuncArray[ ctx ] = ( VoidFunction ) fn;
    }

    using _ContextFunction::validContext;
};

} /* end namespace kregistry */

} /* end namespace scai */
