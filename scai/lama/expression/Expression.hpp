/**
 * @file Expression.hpp
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
 * @brief Definition of template class Expression used for symbolic expressions.
 * @author brandes
 * @date 28.03.2011
 */
#pragma once

namespace scai
{

/**
 * @brief The namespace lama holds everything of the LAMA Library.
 */
namespace lama
{

/**
 * @brief ExpressionTypes expresses the type of an Expression.
 */
enum ExpressionTypes
{
    Plus, Minus, Times, Divide
};

/**
 * @brief The template class Expression represents a mathematical expression.
 *
 * The template class Expression represents a mathematical expression with two
 * operands. The supported operators are defined with ExpressionTypes. The
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
template<typename T1, typename T2, ExpressionTypes type>
class Expression
{
public:
    typedef T1 Arg1Type;
    typedef T2 Arg2Type;
    typedef Expression ExpressionType;
    typedef const ExpressionType ExpressionMemberType;

private:
    const ExpressionTypes mExpressionType;
    const Arg1Type& mArg1;
    const Arg2Type& mArg2;
    // typename Arg1Type::ExpressionMemberType mArg1;
    // typename Arg2Type::ExpressionMemberType mArg2;
public:

    /**
     * @brief This constructor creates a Expression for the given types.
     *
     * @param arg1 the first operand of the expression
     * @param arg2 the second operand of the expression
     */
    Expression( const Arg1Type& arg1, const Arg2Type& arg2 )
        : mExpressionType( type ), mArg1( arg1 ), mArg2( arg2 )
    {
    }
    ;

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
    inline ExpressionTypes getExpressionType() const
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

template<typename ValueType> class Matrix;
template<typename ValueType> class Vector;

// Vector Expressions

/** Symbolic expression 'Scalar * Vector */

template<typename ValueType>
using Expression_SV = Expression<ValueType, Vector<ValueType>, Times>;

/** Symbolic expression 'Vector1 * Vector2 compWise */

template<typename ValueType>
using Expression_VV = Expression<Vector<ValueType>, Vector<ValueType>, Times>;

/** Symbolic expression 'Scalar * Vector * Vector' */

template<typename ValueType>
using Expression_SVV = Expression<ValueType, Expression<Vector<ValueType>, Vector<ValueType>, Times>, Times>;

/** Symbolic expression 'Scalar * Vector + Scalar' */

template<typename ValueType>
using Expression_SV_S = Expression<Expression_SV<ValueType>, ValueType, Plus>;

/** Symbolic expression 'Scalar * Vector + Scalar * Vector' */

template<typename ValueType>
using Expression_SV_SV = Expression<Expression_SV<ValueType>, Expression_SV<ValueType>, Plus>;

/* ============================================================================ */
/*    Matrix expressions                                                        */
/* ============================================================================ */

/** Symbolic expression 'Scalar * Matrix' */

template<typename ValueType>
using Expression_SM = Expression<ValueType, Matrix<ValueType>, Times>;

/* ============================================================================ */
/*    Matrix - Vector expressions                                               */
/* ============================================================================ */

/** Symbolic expression 'Matrix * Vector' */

template<typename ValueType>
using Expression_MV = Expression<Matrix<ValueType>, Vector<ValueType>, Times>;

/** Symbolic expression 'Vector * Matrix' */

template<typename ValueType>
using Expression_VM = Expression<Vector<ValueType>, Matrix<ValueType>, Times>;

/** Symbolic expression 'Scalar * Matrix * Vector' */

template<typename ValueType>
using Expression_SMV = Expression<ValueType, Expression<Matrix<ValueType>, Vector<ValueType>, Times>, Times>;

/** Symbolic expression 'Scalar * Vector * Matrix' */

template<typename ValueType>
using Expression_SVM = Expression<ValueType, Expression<Vector<ValueType>, Matrix<ValueType>, Times>, Times>;

/** Symbolic expression 'Scalar * Matrix * Vector + Scalar * Vector' */

template<typename ValueType>
using Expression_SMV_SV = Expression<Expression_SMV<ValueType>, Expression_SV<ValueType>, Plus>;

/** Symbolic expression 'Scalar * Vector * Matrix + Scalar * Vector' */

template<typename ValueType>
using Expression_SVM_SV = Expression<Expression_SVM<ValueType>, Expression_SV<ValueType>, Plus>;

// Matrix-Matrix Expressions

/** Symbolic expression 'Scalar * Matrix * Matrix' */

template<typename ValueType>
using Expression_SMM = Expression<Expression_SM<ValueType>, Matrix<ValueType>, Times>;

/** Symbolic expression 'Scalar * Matrix + Scalar * Matrix' */

template<typename ValueType>
using Expression_SM_SM = Expression<Expression_SM<ValueType>, Expression_SM<ValueType>, Plus>;

/** Symbolic expression 'Scalar * Matrix * Matrix + Scalar * Matrix' */

template<typename ValueType>
using Expression_SMM_SM = Expression<Expression_SMM<ValueType>, Expression_SM<ValueType>, Plus>;

} /* end namespace lama */

} /* end namespace scai */
