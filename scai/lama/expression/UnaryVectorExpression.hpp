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

/**
 * @brief The template class UnaryExpression represents a unary expression.
 *
 * @param T     the type of the operand of the expression
 */
template<typename T>
class UnaryExpression
{
public:
    typedef T ArgType;
    typedef const UnaryExpression ExpressionMemberType;

private:

    const common::UnaryOp mOp;

    typename ArgType::ExpressionMemberType mArg;

public:

    /**
     * @brief This constructor creates an UnaryExpression for the given type.
     *
     * @param arg is the operand of the expression
     * @param op  specifies the unary operation that is used in the expression
     */
    UnaryExpression( const ArgType& arg, const common::UnaryOp op ) :

        mOp( op ),
        mArg( arg )

    {
    }

    /**
     * @brief The destructor destroys this UnaryExpression.
     */
    virtual ~UnaryExpression()
    {
    }

    /**
     * @brief getExpressionType returns the expression type of this Expression.
     *
     * @return the type of this Expression.
     */
    inline common::UnaryOp getOp() const
    {
        return mOp;
    }

    /**
     * @brief getArg() returns a reference to the operand of this Expression.
     *
     * @return the operand of this Expression.
     */
    inline const ArgType& getArg() const
    {
        return mArg;
    }
};

template<typename ValueType>
class Vector;


/* ============================================================================ */
/*    Vector expressions                                                        */
/* ============================================================================ */

/** Symbolic expression op ( Vector ) */

template<typename ValueType>
using UnaryVectorExpression = UnaryExpression<Vector<ValueType> >;

/** Symbolic expression: apply conj element-wise to the elements. */

template<typename ValueType>
UnaryVectorExpression<ValueType> conj( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::CONJ );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> abs( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::ABS );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> exp( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::EXP );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> log( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::LOG );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> floor( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::FLOOR );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> ceil( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::CEIL );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> sqrt( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::SQRT );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> sin( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::SIN );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> cos( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::COS );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> tan( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::TAN );
}

template<typename ValueType>
UnaryVectorExpression<ValueType> atan( const Vector<ValueType>& v )
{
    return UnaryVectorExpression<ValueType>( v, common::UnaryOp::ATAN );
}


} /* end namespace lama */

} /* end namespace scai */

