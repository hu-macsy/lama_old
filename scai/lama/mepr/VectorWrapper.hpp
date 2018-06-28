/**
 * @file lama/mepr/VectorWrapper.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Wrapper for templated calls in derived vector classes
 * @author Eric Schricker
 * @date 10.03.2016
 */

#pragma once

#include <scai/common/mepr/TypeList.hpp>

namespace scai
{

namespace lama
{

namespace mepr
{

/*
 * This class provides meta-programming to resolve (virtual) methods that
 * take an untyped _Vector object and translate it to a corresponding call of 
 * a typed (template) method.
 *
 * \code
 *    vector.assign( const _Vector& x )
 *    vector.assignImpl<ValueType>( const Vector<ValueType>& x )   
 * \endcode
 */
template<typename Derived, typename TList>
struct VectorWrapper;

/*
 * Termination
 */
template<typename Derived>
struct VectorWrapper<Derived, common::mepr::NullType>
{
    static void assignImpl(
        Derived* obj,
        const _Vector& otherVector )
    {
        COMMON_THROWEXCEPTION( "assign to " << obj << ": unsupported type = " << otherVector.getValueType() );
    }

    static void setDenseValuesImpl(
        Derived* obj,
        const hmemo::_HArray& values )
    {
        COMMON_THROWEXCEPTION( "setDenseValues " << obj << ": unsupported type = " << values.getValueType() );
    }
};

/*
 * Step n
 */
template<typename Derived, typename H, typename T>
struct VectorWrapper<Derived, common::mepr::TypeList<H, T> >
{
    static void assignImpl( Derived* obj, const _Vector& otherVector )
    {
        if ( otherVector.getValueType() == common::getScalarType<H>() )
        {
            obj->assignImpl( reinterpret_cast<const Vector<H>& >( otherVector ) );
        }
        else
        {
            VectorWrapper<Derived, T>::assignImpl( obj, otherVector );
        }
    }

    static void setDenseValuesImpl( Derived* obj, const hmemo::_HArray& values )
    {
        if ( values.getValueType() == common::getScalarType<H>() )
        {
            obj->setDenseValuesImpl( reinterpret_cast<const hmemo::HArray<H>& >( values ) );
        }
        else
        {
            VectorWrapper<Derived, T>::setDenseValuesImpl( obj, values );
        }
    }
};

} /* end namespace mepr */

} /* end namespace lama */

} /* end namespace scai */
