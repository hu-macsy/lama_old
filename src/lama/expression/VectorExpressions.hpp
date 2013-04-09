/**
 * @file VectorExpressions.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief VectorExpressions.hpp
 * @author Jiri Kraus
 * @date 01.06.2011
 * $Id$
 */
#ifndef LAMA_VECTOREXPRESSIONS_HPP_
#define LAMA_VECTOREXPRESSIONS_HPP_

#include <lama/Scalar.hpp>
#include <lama/Vector.hpp>
#include <lama/expression/Expression.hpp>

namespace lama
{

/**
 * @brief The times operator creates an expression that represents Scalar times Vector
 *
 * @param[in] alpha     The scalar.
 * @param[in] x         The vector.
 * @return              The expression representing this difference.
 */

inline Expression<Scalar,Vector,Times> operator*( const Scalar& alpha, const Vector& x )
{
    return Expression<Scalar,Vector,Times>( alpha, x );
}

/**
 * @brief The times operator creates an expression that represents Vector times Scalar
 *
 * @param[in] x         The vector.
 * @param[in] alpha     The scalar.
 * @return              The expression representing this difference.
 */

inline Expression<Scalar,Vector,Times> operator*( const Vector& x, const Scalar& alpha )
{
    return Expression<Scalar,Vector,Times>( alpha, x );
}

/**
 * @brief The plus operator creates an expression that represents the sum
 *        of two vectors.
 *
 * @param[in] x     The first vector.
 * @param[in] y     The second vector.
 * @return          The expression representing this sum.
 */
inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Vector& x,
    const Vector& y )
{
    const Scalar alpha( 1.0 );
    const Scalar beta( 1.0 );

    Expression<Scalar,Vector,Times> exp1( alpha, x );
    Expression<Scalar,Vector,Times> exp2( beta, y );
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( exp1, exp2 );
}

/**
 * @brief The plus operator creates an expression that represents sum of Vector and Vector times Scalar
 *
 * @param[in] x         The vector.
 * @param[in] exp2      The Vector times Scalar.
 * @return              The expression representing this difference.
 */

inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Vector& x,
    const Expression<Scalar,Vector,Times>& exp2 )
{
    Expression<Scalar,Vector,Times> exp1( 1.0, x );
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( exp1, exp2 );
}

/**
 * @brief The plus operator creates an expression that represents sum of Vector times Scalar and Vector times Scalar
 *
 * @param[in] exp1      The vector times Scalar
 * @param[in] exp2      The vector times Scalar
 * @return              The expression representing this difference.
 */

inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Scalar,Vector,Times>& exp1,
    const Expression<Scalar,Vector,Times>& exp2 )
{
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( exp1, exp2 );
}

/**
 * @brief The plus operator creates an expression that represents sum of Vector times Scalar plus Vector
 *
 * @param[in] x         The vector.
 * @param[in] exp2      Expression of Scalar times Vector.
 * @return              The expression representing this difference.
 */

inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Scalar,Vector,Times>& exp2,
    const Vector& x )
{
    Expression<Scalar,Vector,Times> exp1( 1.0, x );
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( exp1, exp2 );
}

/**
 * @brief The plus operator creates an expression that represents sum of Scalar times Vector plus Vector
 *
 * @param[in] x         The vector.
 * @param[in] exp2      Expression of Scalar times Vector.
 * @return              The expression representing this difference.
 */

inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Vector,Scalar,Times>& exp2,
    const Vector& x )
{
    Expression<Scalar,Vector,Times> exp( exp2.getArg2(), exp2.getArg1() );
    Expression<Scalar,Vector,Times> exp1( 1.0, x );
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( exp, exp1 );
}

/* ------------------------------------------------------------------------- */

/**
 * @brief The minus operator creates an expression that represents the difference of Vector and Vector times Scalar
 *
 * @param[in] x         The vector.
 * @param[in] exp2      Expression of Scalar times Vector.
 * @return              The expression representing this difference.
 */

inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Vector& x,
    const Expression<Scalar,Vector,Times>& exp2 )
{
    Expression<Scalar,Vector,Times> exp1( 1.0, x );
    Expression<Scalar,Vector,Times> tmpExp( -exp2.getArg1(), exp2.getArg2() );
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( exp1, tmpExp );
}

/**
 * @brief The minus operator creates an expression that represents the difference
 *        of two scalar vector products.
 *
 * @param[in] exp1     The first expression.
 * @param[in] exp2     The second expression.
 * @return             The expression representing this difference.
 */

inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Expression<Scalar,Vector,Times>& exp1,
    const Expression<Scalar,Vector,Times>& exp2 )
{
    Expression<Scalar,Vector,Times> tmpExp( -exp2.getArg1(), exp2.getArg2() );
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( exp1, tmpExp );
}

/**
 * @brief The minus operator creates an expression that represents the difference
 *        of two vectors.
 *
 * @param[in] x     The first vector.
 * @param[in] y     The second vector.
 * @return          The expression representing this difference.
 */
inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Vector& x,
    const Vector& y )
{
    const Scalar alpha( 1.0 );
    const Scalar beta( -1.0 );

    Expression<Scalar,Vector,Times> exp1( alpha, x );
    Expression<Scalar,Vector,Times> exp2( beta, y );
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( exp1, exp2 );
}

/**
 * @brief The minus operator creates an expression that represents the difference of Vector times Scalar and Vector
 *
 * @param[in] x         The vector.
 * @param[in] exp2      Expression of Scalar times Vector.
 * @return              The expression representing this difference.
 *
 * \code
 *       exp2             x                     return
 *   ( a * vector1 ) - vector1  ->  ( a * vector1 ) +  ( (-1.0) * vectorX )
 * \endcode
 */

inline Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Expression<Scalar,Vector,Times>& exp2,
    const Vector& x )
{
    Expression<Scalar,Vector,Times> tmpExp( exp2.getArg1(), exp2.getArg2() );
    Expression<Scalar,Vector,Times> exp1( -1.0, x );
    return Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>( tmpExp, exp1 );
}

}

#endif // LAMA_VECTOREXPRESSIONS_HPP_
