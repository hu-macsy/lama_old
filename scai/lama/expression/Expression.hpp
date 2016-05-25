/**
 * @file Expression.hpp
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
template<typename T1,typename T2,ExpressionTypes type>
class Expression
{
public:
    typedef T1 Arg1Type;
    typedef T2 Arg2Type;
    typedef Expression ExpressionType;
    typedef const ExpressionType ExpressionMemberType;

private:
    const ExpressionTypes mExpressionType;
    typename Arg1Type::ExpressionMemberType mArg1;
    typename Arg2Type::ExpressionMemberType mArg2;
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

class Scalar;
class Vector;
class Matrix;

/** Symbolic expression 'Scalar * Vector' */

typedef Expression<Scalar,Vector,Times> Expression_SV;

/** Symbolic expression 'Scalar * Vector + Scalar * Vector' */

typedef Expression<Expression_SV,Expression_SV,Plus> Expression_SV_SV;

typedef Expression<Matrix,Vector,Times> Expression_MV;

typedef Expression<Vector,Matrix,Times> Expression_VM;

typedef Expression<Scalar,Expression<Matrix,Vector,Times>,Times> Expression_SMV;

typedef Expression<Scalar,Expression<Vector,Matrix,Times>,Times> Expression_SVM;

typedef Expression<Expression_SMV,Expression_SV,Plus> Expression_SMV_SV;

typedef Expression<Expression_SVM,Expression_SV,Plus> Expression_SVM_SV;

typedef Expression<Scalar,Matrix,Times> Expression_SM;

typedef Expression<Expression_SM,Matrix,Times> Expression_SMM;

typedef Expression<Expression_SM,Expression_SM,Plus> Expression_SM_SM;

typedef Expression<Expression_SMM,Expression_SM,Plus> Expression_SMM_SM;

} /* end namespace lama */

} /* end namespace scai */
