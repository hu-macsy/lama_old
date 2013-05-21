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
 * @author Thomas Brandes
 * @date 28.03.2011
 * $Id$
 */
#ifndef LAMA_MATRIXVECTOREXPRESSIONS_HPP_
#define LAMA_MATRIXVECTOREXPRESSIONS_HPP_

#include <lama/matrix/Matrix.hpp>

#include <lama/Vector.hpp>
#include <lama/Scalar.hpp>

#include <lama/expression/Expression.hpp>

namespace lama
{

/* ------------------------------------------------------------------------- */

/**
 * @brief This times operator creates an expression that represents the product
 *        of a Matrix and a Vector.
 *
 * The times operator creates an Expression that represents the product of the
 * a matrix and a vector. To give an example this Expression is
 * used by a assignment operator Vector to assign the result of this expression to a Vector
 * without generating any unnecessary temporaries.
 *
 * @param[in] matrix  The input matrix.
 * @param[in] vector  The input vector.
 * @return            The expression representing this product.
 */
inline Expression<Matrix,Vector,Times> operator*( const Matrix& matrix, const Vector& vector )
{
    return Expression<Matrix,Vector,Times>( matrix, vector );
}

// alpha*(A*x)
/**
 * @brief This plus operator creates an expression that represents a * A * x, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] scalar  The Scalar.
 * @param[in] exp     The Expression A*x
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Expression<Matrix,Vector,Times>,Times> operator*(
    const Scalar& scalar,
    const Expression<Matrix,Vector,Times>& exp )
{
    return Expression<Scalar,Expression<Matrix,Vector,Times>,Times>( scalar, exp );
}

/**
 * @brief This times operator creates an expression that represents a * A * x, where
 *        x is vector, A is a matrix and a is a scalar
 *
 * @param[in] exp     The Expression A*x
 * @param[in] scalar  The Scalar.
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Expression<Matrix,Vector,Times>,Times> operator*(
    const Expression<Matrix,Vector,Times>& exp,
    const Scalar& scalar )
{
    return Expression<Scalar,Expression<Matrix,Vector,Times>,Times>( scalar, exp );
}
/**
 * @brief This times operator creates an expression that represents the product
 *        of a Matrix, Vector and Scalar.
 *
 * @param[in] matrix  The matrix.
 * @param[in] exp     The expression a*x
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Expression<Matrix,Vector,Times>,Times> operator*(
    const Matrix& matrix,
    const Expression<Scalar,Vector,Times>& exp )
{
    return Expression<Scalar,Expression<Matrix,Vector,Times>,Times>(
               exp.getArg1(), Expression<Matrix,Vector,Times>( matrix, exp.getArg2() ) );
}

/**
 * @brief This plus operator creates an expression that represents A * B * x, where
 *        x is vector, A and B are matrices.
 *
 * @param[in] exp     The expression a*A
 * @param[in] vector  The Vector.
 * @return            The expression representing this product.
 */
inline Expression<Scalar,Expression<Matrix,Vector,Times>,Times> operator*(
    const Expression<Scalar,Matrix,Times>& exp,
    const Vector& vector )
{
    return Expression<Scalar,Expression<Matrix,Vector,Times>,Times>(
               exp.getArg1(), Expression<Matrix,Vector,Times>( exp.getArg2(), vector ) );
}

/* ------------------------------------------------------------------------- */

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression a * A * x.
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& exp,
    const Vector& vector )
{
    Expression<Scalar,Vector,Times> exp1( Scalar( -1.0 ), vector );

    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp, exp1 );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x * a.
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Expression<Matrix,Expression<Vector,Scalar,Times>,Times>& exp,
    const Vector& vector )
{
    Expression<Scalar,Vector,Times> exp1( Scalar( -1.0 ), vector );

    Expression<Matrix,Vector,Times> e( exp.getArg1(), exp.getArg2().getArg1() );
    Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp2( exp.getArg2().getArg2(), e );
    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp2, exp1 );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp1    The expression a*A*x
 * @param[in] exp2    The expression b*y
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& exp1,
    const Expression<Scalar,Vector,Times>& exp2 )
{
    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp1, Expression<Scalar,Vector,Times>( exp2.getArg1() * -1.0, exp2.getArg2() ) );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp2    The expression a*x
 * @param[in] exp1    The expression b*A*y
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Expression<Scalar,Vector,Times>& exp2,
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& exp1 )
{

    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               Expression<Scalar,Expression<Matrix,Vector,Times>,Times>( exp1.getArg1() * -1.0, exp1.getArg2() ),
               exp2 );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Vector& vector,
    const Expression<Matrix,Vector,Times>& exp )
{
    const Expression<Scalar,Vector,Times> exp1( Scalar( 1.0 ), vector );
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp2( Scalar( -1.0 ), exp );

    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp2, exp1 );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Expression<Matrix,Vector,Times>& exp,
    const Vector& vector )
{
    const Expression<Scalar,Vector,Times> exp1( Scalar( -1.0 ), vector );
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp2( Scalar( 1.0 ), exp );

    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp2, exp1 );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this sum.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator-(
    const Expression<Matrix,Expression<Scalar,Vector,Times>,Times>& exp,
    const Vector& vector )
{
    const Expression<Scalar,Vector,Times> exp1( Scalar( -1.0 ), vector );
    Expression<Matrix,Vector,Times> e( exp.getArg1(), exp.getArg2().getArg2() );
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp2( exp.getArg2().getArg1(), e );

    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp2, exp1 );
}

/* ------------------------------------------------------------------------- */

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this addition.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Vector& vector,
    const Expression<Matrix,Vector,Times>& exp )
{
    const Expression<Scalar,Vector,Times> exp1( Scalar( 1.0 ), vector );
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp2( Scalar( 1.0 ), exp );
    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp2, exp1 );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression a * A * y.
 * @return            The expression representing this addition.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& exp,
    const Vector& vector )
{
    Expression<Scalar,Vector,Times> exp2( 1.0, vector );
    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp, exp2 );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * a * x.
 * @return            The expression representing this addition.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Matrix,Expression<Scalar,Vector,Times>,Times>& exp,
    const Vector& vector )
{
    Expression<Matrix,Vector,Times> exp1( exp.getArg1(), exp.getArg2().getArg2() );
    Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp2( exp.getArg2().getArg1(), exp1 );
    Expression<Scalar,Vector,Times> exp3( 1.0, vector );

    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp2, exp3 );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars.
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x * a.
 * @return            The expression representing this addition.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Matrix,Expression<Vector,Scalar,Times>,Times>& exp,
    const Vector& vector )
{
    Expression<Matrix,Vector,Times> exp1( exp.getArg1(), exp.getArg2().getArg1() );
    Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp2( exp.getArg2().getArg2(), exp1 );
    Expression<Scalar,Vector,Times> exp3( 1.0, vector );

    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp2, exp3 );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars.
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this addition.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Matrix,Vector,Times>& exp,
    const Vector& vector )
{
    Expression<Matrix,Vector,Times> exp1( exp.getArg1(), exp.getArg2() );
    Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp2( 1.0, exp1 );
    Expression<Scalar,Vector,Times> exp3( 1.0, vector );
    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp2, exp3 );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp1    The expression a*A*x
 * @param[in] exp2    The expression b*y
 * @return            The expression representing this product.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& exp1,
    const Expression<Scalar,Vector,Times>& exp2 )
{
    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp1, exp2 );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp2    The expression a*x
 * @param[in] exp1    The expression b*A*y
 * @return            The expression representing this product.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Scalar,Vector,Times>& exp2,
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& exp1 )
{
    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp1, exp2 );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp2    The expression a*x
 * @param[in] exp1    The expression b*A*y
 * @return            The expression representing this product.
 */
inline Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> operator+(
    const Expression<Matrix,Vector,Times>& exp2,
    const Expression<Scalar,Vector,Times>& exp1 )
{
    Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp( 1.0, exp2 );
    return Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>(
               exp, exp1 );
}

}

#endif // LAMA_MATRIXVECTOREXPRESSIONS_HPP_
