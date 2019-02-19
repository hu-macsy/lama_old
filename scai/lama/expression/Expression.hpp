/**
 * @file Expression.hpp
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
 * @brief Definition of template class Expression used for symbolic expressions.
 * @author Thomas Brandes
 * @date 28.03.2011
 */
#pragma once

#include <scai/common/BinaryOp.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

/**
 * @brief The namespace lama holds everything of the LAMA Library.
 */
namespace lama
{

/**
 * @brief The template class Expression represents a binary expression for two arguments.
 *
 * The template class Expression represents a mathematical expression with two
 * operands. The supported operators are defined with BinaryOp. The
 * expression holds references to the two operands.
 *
 * @param T1    the type of the first operand of the expression
 * @param T2    the type of the second operand of the expression
 * @param type  the type of the operation
 *
 * An object of class expression is comparable to an abstract syntax tree
 * for a binary expression. Operations on this object are just the constructor,
 * the getter for the expression type and the two getters for the operands.
 */
template<typename T1, typename T2, common::BinaryOp type>
class Expression
{
public:
    typedef T1 Arg1Type;
    typedef T2 Arg2Type;
    typedef Expression ExpressionType;
    typedef const ExpressionType ExpressionMemberType;

private:

    const common::BinaryOp mExpressionType;

    typename Arg1Type::ExpressionMemberType mArg1;
    typename Arg2Type::ExpressionMemberType mArg2;

public:

    /**
     * @brief This constructor creates a Expression for the given types.
     *
     * @param arg1 the first operand of the expression
     * @param arg2 the second operand of the expression
     */
    Expression( const Arg1Type& arg1, const Arg2Type& arg2 ) : 

        mExpressionType( type ), 
        mArg1( arg1 ), 
        mArg2( arg2 )

    {
    }

    /**
     * @brief The destructor destroys this Expression.
     */
    virtual ~Expression()
    {
    }

    /**
     * @brief getExpressionType returns the expression type of this Expression.
     *
     * @return the type of this Expression.
     */
    inline common::BinaryOp getExpressionType() const
    {
        return mExpressionType;
    }

    /**
     * @brief getArg1() returns a reference to the first operand of this Expression.
     *
     * @return the first operand of this Expression.
     */
    inline const Arg1Type& getArg1() const
    {
        return mArg1;
    }

    /**
     * @brief getArg2() returns a reference to the second operand of this Expression.
     *
     * @return the second operand of this Expression.
     */
    inline const Arg2Type& getArg2() const
    {
        return mArg2;
    }
};

namespace intern
{
    class Scalar;
}

template<typename ValueType> class Matrix;
template<typename ValueType> class Vector;

/* ============================================================================ */
/*    Vector expressions                                                        */
/* ============================================================================ */

/** Symbolic expression 'Scalar * Vector */

template<typename ValueType>
using Expression_SV = Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MULT>;

/** Symbolic expression 'Vector1 * Vector2 element-wise */

template<typename ValueType>
using Expression_VV = Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::MULT>;

/** Symbolic expression 'Scalar * Vector * Vector', element-wise */

template<typename ValueType>
using Expression_SVV = Expression<intern::Scalar, Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::MULT>, common::BinaryOp::MULT>;

/** Symbolic expression 'Scalar * Vector + Scalar' */

template<typename ValueType>
using Expression_SV_S = Expression<Expression_SV<ValueType>, intern::Scalar, common::BinaryOp::ADD>;

/** Symbolic expression 'Scalar * Vector + Scalar * Vector' */

template<typename ValueType>
using Expression_SV_SV = Expression<Expression_SV<ValueType>, Expression_SV<ValueType>, common::BinaryOp::ADD>;

/* ============================================================================ */
/*    OpMatrix =   matrix, transpose( matrix ), conj( matrix ), ...             */
/* ============================================================================ */

template <typename ValueType>
class OpMatrix
{
public:

    typedef const OpMatrix<ValueType> ExpressionMemberType;

    OpMatrix( const Matrix<ValueType>& matrix, const common::MatrixOp op ) :
        mMatrix( matrix ),
        mOp( op )
    {
    }

    OpMatrix<ValueType>( const OpMatrix<ValueType>& ) = default;

    OpMatrix<ValueType>& operator=( const OpMatrix<ValueType>& ) = delete;

    const Matrix<ValueType>& getMatrix() const 
    {
        return mMatrix;
    }

    common::MatrixOp getOp() const 
    {
        return mOp;
    }

private:

    const Matrix<ValueType>& mMatrix;
    const common::MatrixOp mOp;
};

/* ============================================================================ */
/*    Matrix expressions                                                        */
/* ============================================================================ */

/** Symbolic expression 'Scalar * Matrix' */

template<typename ValueType>
using Expression_SM = Expression<intern::Scalar, OpMatrix<ValueType>, common::BinaryOp::MULT>;

/** Symbolic expression 'Scalar * Matrix * Matrix' */

template<typename ValueType>
using Expression_SMM = Expression<Expression_SM<ValueType>, OpMatrix<ValueType>, common::BinaryOp::MULT>;

/** Symbolic expression 'Scalar * Matrix + Scalar * Matrix' */

template<typename ValueType>
using Expression_SM_SM = Expression<Expression_SM<ValueType>, Expression_SM<ValueType>, common::BinaryOp::ADD>;

/** Symbolic expression 'Scalar * Matrix * Matrix + Scalar * Matrix' */

template<typename ValueType>
using Expression_SMM_SM = Expression<Expression_SMM<ValueType>, Expression_SM<ValueType>, common::BinaryOp::ADD>;

/* ============================================================================ */
/*    Matrix - Vector expressions                                               */
/* ============================================================================ */

/** Symbolic expression '(op)Matrix * Vector' */

template<typename ValueType>
using Expression_MV = Expression<OpMatrix<ValueType>, Vector<ValueType>, common::BinaryOp::MULT>;

/** Symbolic expression 'Scalar * (op)Matrix * Vector' */

template<typename ValueType>
using Expression_SMV = Expression<intern::Scalar, Expression_MV<ValueType>, common::BinaryOp::MULT>;

/** Symbolic expression 'Scalar * (op)Matrix * Vector + Scalar * Vector' */

template<typename ValueType>
using Expression_SMV_SV = Expression<Expression_SMV<ValueType>, Expression_SV<ValueType>, common::BinaryOp::ADD>;

} /* end namespace lama */

} /* end namespace scai */
