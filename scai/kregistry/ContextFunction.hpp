/**
 * @file ContextFunction.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Definition of class for typed function pointers, one for each context
 * @author Thomas Brandes
 * @date 14.10.2015
 */
#pragma once

#include <scai/common/ContextType.hpp>
#include <scai/common/SCAITypes.hpp>

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

class _ContextFunction : public common::context
{
public:

    /** Default constructor initializes all function pointers with NULL */

    _ContextFunction();

    _ContextFunction( const _ContextFunction& other );

    void clear();   // sets all function pointers to NULL

    void assign( const _ContextFunction& other );

    inline VoidFunction get( ContextType ctx ) const
    {
        return mContextFuncArray[ ctx ];
    }

    inline void set( ContextType ctx, VoidFunction fn )
    {
        mContextFuncArray[ ctx ] = fn;
    }

    inline bool isEmpty() const
    {
    	bool status = true;
    	for(IndexType i = 0; i < common::context::MaxContext; ++i)
    	{
    		if ( mContextFuncArray[i] != NULL )
    		{
    			status = false;
    			break;
    		}
    	}

    	return status;
    }

    std::string printIt() const;

    common::context::ContextType validContext( ContextType preferedCtx ) const;

    ContextType validContext( const _ContextFunction& other, ContextType preferedCtx ) const;

    ContextType validContext( const _ContextFunction& other1, 
                                      const _ContextFunction& other2, 
                                      ContextType preferedCtx ) const;

protected:

    // array with function pointer for each context

    VoidFunction mContextFuncArray[MaxContext];
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

    FunctionType get( ContextType ctx ) const
    {
        return ( FunctionType ) mContextFuncArray[ ctx ];
    }

    // provide typed set

    void set( ContextType ctx, FunctionType fn )
    {
        mContextFuncArray[ ctx ] = ( VoidFunction ) fn;
    }

    using _ContextFunction::validContext;
};

} /* end namespace kregistry */

} /* end namespace scai */
