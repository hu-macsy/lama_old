/**
 * @file ContextFunction.cpp
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
 * @brief Implementation of routines for ContextFunction
 * @author Thomas Brandes
 * @date 13.10.2015
 */

#include <scai/kregistry/ContextFunction.hpp>

#include <sstream>

namespace scai
{

using namespace common;

namespace kregistry
{

_ContextFunction:: _ContextFunction()
{
    clear();
}

void _ContextFunction::clear()
{
    for ( int i = 0; i < context::MaxContext; ++i )
    {
        mContextFuncArray[i] =  NULL;
    }
}

_ContextFunction::_ContextFunction( const _ContextFunction& other )
{
    assign( other );
}

void _ContextFunction::assign( const _ContextFunction& other )
{
    for ( int i = 0; i < context::MaxContext; ++i )
    {
        mContextFuncArray[i] =  other.mContextFuncArray[i];
    }
}

context::ContextType _ContextFunction::validContext( context::ContextType preferedCtx )
{
    if ( mContextFuncArray[preferedCtx] != NULL )
    {
        return preferedCtx;
    }

    for ( int i = 0; i < context::MaxContext; ++i )
    {
        if ( mContextFuncArray[i] != NULL )
        {
            return static_cast<context::ContextType>( i );
        }
    }

    // throw exception

    return static_cast<context::ContextType>( context::MaxContext );
}

context::ContextType _ContextFunction::validContext( const _ContextFunction& other, context::ContextType preferedCtx )
{
    if ( mContextFuncArray[preferedCtx] != NULL && other.mContextFuncArray[preferedCtx] != NULL )
    {
        // std::cout << "both valid at context = " << preferedCtx << std::endl;

        return preferedCtx;
    }

    for ( int i = 0; i < context::MaxContext; ++i )
    {
        if ( mContextFuncArray[i] != NULL && other.mContextFuncArray[i] != NULL )
        {
            return static_cast<context::ContextType>( i );
        }
    }

    // throw exception

    return static_cast<ContextType>( context::MaxContext );
}

_ContextFunction::ContextType _ContextFunction::validContext( 
    const _ContextFunction& other1, 
    const _ContextFunction& other2,
    ContextType preferedCtx )
{
    if ( mContextFuncArray[preferedCtx] != NULL && other1.mContextFuncArray[preferedCtx] != NULL 
                                                && other2.mContextFuncArray[preferedCtx] != NULL )
    {
        // std::cout << "both valid at context = " << preferedCtx << std::endl;

        return preferedCtx;
    }

    for ( int i = 0; i < context::MaxContext; ++i )
    {
        if ( mContextFuncArray[i] != NULL && other1.mContextFuncArray[i] != NULL && other2.mContextFuncArray[i] )
        {
            return static_cast<ContextType>( i );
        }
    }

    // throw exception

    return static_cast<ContextType>( context::MaxContext );
}

std::string _ContextFunction::printIt() const
{
    std::ostringstream msg;

    for ( int i = 0; i < context::MaxContext; ++i )
    {
        if ( mContextFuncArray[i] != NULL )
        {
            msg << "+" << static_cast<ContextType>( i ) << "+";
        }
    }

    return msg.str();
}

} /* end namespace kregistry */

} /* end namespace scai */

