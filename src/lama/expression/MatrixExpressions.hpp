/**
 * @file MatrixVectorExpressions.hpp
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
 * @brief MatrixVectorExpressions.hpp
 * @author brandes
 * @date 28.03.2011
 * $Id$
 */
#ifndef LAMA_MATRIXEXPRESSIONS_HPP_
#define LAMA_MATRIXEXPRESSIONS_HPP_

#include <lama/matrix/Matrix.hpp>

#include <lama/Vector.hpp>
#include <lama/Scalar.hpp>

#include <lama/expression/Expression.hpp>

namespace lama
{

/**
 * @brief The add operator to add two matrixes
 *
 * @param[in] matrixA The first matrix.
 * @param[in] matrixB The second matrix.
 * @return            The sum of both matrices.
 */
inline Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus> operator+(
    const Matrix& matrixA,
    const Matrix& matrixB )
{
    const Scalar alpha( 1.0 );
    const Scalar beta( 1.0 );
    Expression<Scalar,Matrix,Times> exp1( alpha, matrixA );
    Expression<Scalar,Matrix,Times> exp2( beta, matrixB );
    return Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus>( exp1, exp2 );
}

/**
 * @brief The minus operator to sub two matrixes
 *
 * @param[in] matrixA The first matrix.
 * @param[in] matrixB The second matrix.
 * @return            The difference of both matrices.
 */
inline Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus> operator-(
    const Matrix& matrixA,
    const Matrix& matrixB )
{
    const Scalar alpha( 1.0 );
    const Scalar beta( -1.0 );
    Expression<Scalar,Matrix,Times> exp1( alpha, matrixA );
    Expression<Scalar,Matrix,Times> exp2( beta, matrixB );
    return Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus>( exp1, exp2 );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix.
 *
 * @param[in] scalar  The input scalar.
 * @param[in] matrix  The input matrix.
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Matrix,Times> operator*( const Scalar& scalar, const Matrix& matrix )
{
    return Expression<Scalar,Matrix,Times>( scalar, matrix );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix.
 *
 * @param[in] matrix  The input matrix.
 * @param[in] scalar  The input scalar.
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Matrix,Times> operator*( const Matrix& matrix, const Scalar& scalar )
{
    return Expression<Scalar,Matrix,Times>( scalar, matrix );
}

//A*B
/**
 * @brief This times operator creates an expression that represents the product
 *        of a Matrix and a Matrix.
 *
 * @param[in] m1      The first input matrix.
 * @param[in] m2      The second input matrix.
 * @return            The expression representing this product.
 */
inline Expression<Matrix,Matrix,Times> operator*( const Matrix& m1, const Matrix& m2 )
{
    return Expression<Matrix,Matrix,Times>( m1, m2 );
}

//alpha times A times B

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix.
 *
 * @param[in] m1      The first input matrix.
 * @param[in] exp     The expression scalar times matrix.
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> operator*(
    const Matrix& m1,
    const Expression<Scalar,Matrix,Times>& exp )
{
    return Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>(
               exp.getArg1(), Expression<Matrix,Matrix,Times>( m1, exp.getArg2() ) );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix.
 *
 * @param[in] exp     The expression scalar times matrix.
 * @param[in] m1      The first input matrix.
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> operator*(
    const Expression<Scalar,Matrix,Times>& exp,
    const Matrix& m1 )
{
    return Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>(
               exp.getArg1(), Expression<Matrix,Matrix,Times>( exp.getArg2(), m1 ) );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix.
 *
 * @param[in] s1      The Scalar
 * @param[in] exp     The expression scalar times matrix.
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> operator*(
    const Scalar& s1,
    const Expression<Matrix,Matrix,Times>& exp )
{
    return Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>( s1, exp );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix.
 *
 * @param[in] exp     The expression scalar times matrix.
 * @param[in] s1      The Scalar
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> operator*(
    const Expression<Matrix,Matrix,Times>& exp,
    const Scalar& s1 )
{
    return Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>( s1, exp );
}

//alpha*A*B + beta*C

/**
 * @brief This plus operator creates an expression that represents the sum
 *        of Scalar times Matrix times Matrix plus Scalar times Matrix.
 *
 * @param[in] exp1    The expression Scalar times Matrix times Matrix.
 * @param[in] exp2    The expression Scalar times Matrix.
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> operator+(
    const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& exp1,
    const Expression<Scalar,Matrix,Times>& exp2 )
{
    return Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus>(
               exp1, exp2 );
}

/**
 * @brief This plus operator creates an expression that represents the sum
 *        of Scalar times Matrix times Matrix plus Scalar times Matrix.
 *
 * @param[in] exp2    The expression Scalar times Matrix.
 * @param[in] exp1    The expression Scalar times Matrix times Matrix.
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> operator+(
    const Expression<Scalar,Matrix,Times>& exp2,
    const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& exp1 )
{
    return Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus>(
               exp1, exp2 );
}

/**
 * @brief This times operator creates an expression that represents the sum
 *        of Scalar*(Matrix*Matrix) + Scalar*Matrix.
 *
 * @param[in] exp1    The expression Scalar*(Matrix*Matrix).
 * @param[in] exp2    The expression Scalar*Matrix.
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> operator-(
    const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& exp1,
    const Expression<Scalar,Matrix,Times>& exp2 )
{
    return Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus>(
               exp1, Expression<Scalar,Matrix,Times>( exp2.getArg1() * -1.0, exp2.getArg2() ) );
}

/**
 * @brief This times operator creates an expression that represents the sum
 *        of Scalar*(Matrix*Matrix) + Scalar*Matrix.
 *
 * @param[in] exp2    The expression Scalar*Matrix.
 * @param[in] exp1    The expression Scalar*(Matrix*Matrix).
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> operator-(
    const Expression<Scalar,Matrix,Times>& exp2,
    const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& exp1 )
{
    Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> e( exp1.getArg1() * -1.0, exp1.getArg2() );
    return Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus>(
               e, exp2 );
}

}

#endif // LAMA_MATRIXEXPRESSIONS_HPP_
