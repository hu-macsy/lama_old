/**
 * @file Expression.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Definition of template class Expression used for symbolic expressions.
 * @author brandes
 * @date 28.03.2011
 * @since 1.0.0
 */
#ifndef LAMA_EXPRESSION_HPP_
#define LAMA_EXPRESSION_HPP_

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

typedef Expression<Scalar, Vector, Times> Expression_SV;

/** Symbolic expression 'Scalar * Vector + Scalar * Vector' */

typedef Expression<Expression_SV, Expression_SV, Plus> Expression_SV_SV;

typedef Expression<Matrix, Vector, Times> Expression_MV;

typedef Expression<Scalar, Expression<Matrix, Vector, Times>, Times> Expression_SMV;

typedef Expression<Expression_SMV, Expression_SV, Plus> Expression_SMV_SV;

typedef Expression<Scalar, Matrix, Times> Expression_SM;

typedef Expression<Expression_SM, Matrix, Times> Expression_SMM;

typedef Expression<Expression_SM, Expression_SM, Plus> Expression_SM_SM;

typedef Expression<Expression_SMM, Expression_SM, Plus> Expression_SMM_SM;

} //namespace lama

#endif // LAMA_EXPRESSION_HPP_
