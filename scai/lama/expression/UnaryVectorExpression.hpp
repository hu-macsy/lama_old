/**
 * @file UnaryVectorExpression.hpp
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
 * @brief Functions to build symbolic expressions unaryFunction( Vector<ValueType> )
 * @author Thomas Brandes
 * @date 21.11.2017
 */
#pragma once

#include <scai/common/UnaryOp.hpp>

namespace scai
{

namespace lama
{

template<typename ValueType>
class Vector;

template<typename ValueType>
class UnaryVectorExpression
{
public:
    common::UnaryOp mOp;
    const Vector<ValueType>& mV;

    UnaryVectorExpression( common::UnaryOp op, const Vector<ValueType>& v ) : mOp( op ), mV( v ) 
    {}
};

/** Symbolic expression: apply conj element-wise to the elements. */

template<typename ValueType>
UnaryVectorExpression<ValueType> conj( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::CONJ, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> abs( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::ABS, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> exp( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::EXP, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> log( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::LOG, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> floor( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::FLOOR, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> ceil( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::CEIL, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> sqrt( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::SQRT, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> sin( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::SIN, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> cos( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::COS, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> tan( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::TAN, v );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> atan( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( common::UnaryOp::ATAN, v );
}

} /* end namespace lama */

} /* end namespace scai */

