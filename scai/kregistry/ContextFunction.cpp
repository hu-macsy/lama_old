/**
 * @file ContextFunction.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
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

using common::ContextType;

namespace kregistry
{

_ContextFunction:: _ContextFunction()
{
    clear();
}

void _ContextFunction::clear()
{
    for ( int i = 0; i < static_cast<int>( ContextType::MaxContext ); ++i )
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
    for ( int i = 0; i < static_cast<int>( ContextType::MaxContext ); ++i )
    {
        mContextFuncArray[i] =  other.mContextFuncArray[i];
    }
}

ContextType _ContextFunction::validContext( ContextType preferedCtx ) const
{
    if ( mContextFuncArray[static_cast<int>( preferedCtx )] != NULL )
    {
        return preferedCtx;
    }

    for ( int i = 0; i < static_cast<int>( ContextType::MaxContext ); ++i )
    {
        if ( mContextFuncArray[i] != NULL )
        {
            return ContextType( i );
        }
    }

    // throw exception
    return ContextType::MaxContext;
}

ContextType _ContextFunction::validContext( const _ContextFunction& other, ContextType preferedCtx ) const
{
    int p = static_cast<int>( preferedCtx );

    if ( mContextFuncArray[p] != NULL && other.mContextFuncArray[p] != NULL )
    {
        // std::cout << "both valid at context = " << preferedCtx << std::endl;
        return preferedCtx;
    }

    for ( int i = 0; i < static_cast<int>( ContextType::MaxContext ); ++i )
    {
        if ( mContextFuncArray[i] != NULL && other.mContextFuncArray[i] != NULL )
        {
            return ContextType( i );
        }
    }

    // throw exception
    return ContextType::MaxContext;
}

ContextType _ContextFunction::validContext(
    const _ContextFunction& other1,
    const _ContextFunction& other2,
    ContextType preferedCtx ) const
{
    int p = static_cast<int>( preferedCtx );

    if ( mContextFuncArray[p] != NULL && other1.mContextFuncArray[p] != NULL && other2.mContextFuncArray[p] != NULL )
    {
        return preferedCtx;
    }

    for ( int i = 0; i < static_cast<int>( ContextType::MaxContext ); ++i )
    {
        if ( mContextFuncArray[i] != NULL && other1.mContextFuncArray[i] != NULL && other2.mContextFuncArray[i] )
        {
            return ContextType( i );
        }
    }

    // throw exception
    return ContextType::MaxContext;
}

std::string _ContextFunction::printIt() const
{
    std::ostringstream msg;

    for ( int i = 0; i < static_cast<int>( ContextType::MaxContext ); ++i )
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

