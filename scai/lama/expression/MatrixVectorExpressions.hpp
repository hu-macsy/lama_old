/**
 * @file MatrixVectorExpressions.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Operators to form symbolic expressions Scalar * Matrix * Vector + Scalar * Vector
 * @author Thomas Brandes
 * @date 28.03.2011
 * @since 1.0.0
 */
#pragma once

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/Vector.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/common/Constants.hpp>

#include <scai/lama/expression/Expression.hpp>

namespace scai
{

using common::Constants;

namespace lama
{

/* ------------------------------------------------------------------------- */
/*  Expressions return 'Scalar * Matrix * Vector'                            */
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
inline Expression_SMV operator*( const Matrix& matrix, const Vector& vector )
{
    return Expression_SMV( Scalar( Constants<IndexType>::one ), Expression_MV( matrix, vector ) );
}

inline Expression_SVM operator*( const Vector& vector, const Matrix& matrix )
{
    return Expression_SVM( Scalar( Constants<IndexType>::one ), Expression_VM( vector, matrix ) );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] scalar  The Scalar.
 * @param[in] exp     The Expression A*x
 * @return            The expression representing this product.
 */
inline Expression_SMV operator*( const Scalar& scalar, const Expression_MV& exp )
{
    return Expression_SMV( scalar, exp );
}

/**
 * @brief This plus operator creates an expression that represents a * x * A, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] scalar  The Scalar.
 * @param[in] exp     The Expression x*A
 * @return            The expression representing this product.
 */
inline Expression_SVM operator*( const Scalar& scalar, const Expression_VM& exp )
{
    return Expression_SVM( scalar, exp );
}

inline Expression_SMV operator*( const Scalar& scalar, const Expression_SMV& exp )
{
    const Expression_MV& mv = exp.getArg2();

    return Expression_SMV( scalar * exp.getArg1(), mv );
}

inline Expression_SVM operator*( const Scalar& scalar, const Expression_SVM& exp )
{
    const Expression_VM& vm = exp.getArg2();

    return Expression_SVM( scalar * exp.getArg1(), vm );
}

inline Expression_SMV operator*( const Expression_SMV& exp, const Scalar& scalar )
{
    return Expression_SMV( exp.getArg1() * scalar, exp.getArg2() );
}

inline Expression_SVM operator*( const Expression_SVM& exp, const Scalar& scalar )
{
    return Expression_SVM( exp.getArg1() * scalar, exp.getArg2() );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of a Matrix, Vector and Scalar.
 *
 * @param[in] matrix  The matrix.
 * @param[in] exp     The expression a*x
 * @return            The expression representing this product.
 */
inline Expression_SMV operator*( const Matrix& matrix, const Expression_SV& exp )
{
    return Expression_SMV( exp.getArg1(), Expression_MV( matrix, exp.getArg2() ) );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of a Vector, Scalar and M .
 *
 * @param[in] exp     The expression a*x
 * @param[in] matrix  The matrix.
 * @return            The expression representing this product.
 */
inline Expression_SVM operator*( const Expression_SV& exp, const Matrix& matrix )
{
    return Expression_SVM( exp.getArg1(), Expression_VM( exp.getArg2(), matrix ) );
}

/**
 * @brief This plus operator creates an expression that represents A * B * x, where
 *        x is vector, A and B are matrices.
 *
 * @param[in] exp     The expression a*A
 * @param[in] vector  The Vector.
 * @return            The expression representing this product.
 */
inline Expression_SMV operator*( const Expression_SM& exp, const Vector& vector )
{
    return Expression_SMV( exp.getArg1(), Expression_MV( exp.getArg2(), vector ) );
}

inline Expression_SVM operator*( const Vector& vector, const Expression_SM& exp )
{
    return Expression_SVM( exp.getArg1(), Expression_VM( vector, exp.getArg2() ) );
}

/* ------------------------------------------------------------------------- */
/*  Expressions 'Scalar * Matrix * Vector' +/- vector                        */
/* ------------------------------------------------------------------------- */

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression a * A * x.
 * @return            The expression representing this sum.
 */
inline Expression_SMV_SV operator-( const Expression_SMV& exp, const Vector& vector )
{
    return Expression_SMV_SV( exp, Expression_SV( Scalar( Constants<IndexType>::minusone ), vector ) );
}

/**
 * @brief This minus operator creates an expression that represents a * x * A + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression a * x * A.
 * @return            The expression representing this sum.
 */
inline Expression_SVM_SV operator-( const Expression_SVM& exp, const Vector& vector )
{
    return Expression_SVM_SV( exp, Expression_SV( Scalar( Constants<IndexType>::minusone ), vector ) );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this sum.
 */
inline Expression_SMV_SV operator-( const Vector& vector, const Expression_SMV& exp )
{
    Expression_SMV minusExp( -exp.getArg1(), exp.getArg2() );

    return Expression_SMV_SV( minusExp, Expression_SV( Scalar( Constants<IndexType>::one ), vector ) );
}

/**
 * @brief This minus operator creates an expression that represents a * x * A + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression x * A .
 * @return            The expression representing this sum.
 */
inline Expression_SVM_SV operator-( const Vector& vector, const Expression_SVM& exp )
{
    Expression_SVM minusExp( -exp.getArg1(), exp.getArg2() );

    return Expression_SVM_SV( minusExp, Expression_SV( Scalar( Constants<IndexType>::one ), vector ) );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this addition.
 */
inline Expression_SMV_SV operator+( const Vector& vector, const Expression_SMV& exp )
{
    return Expression_SMV_SV( exp, Expression_SV( Scalar( Constants<IndexType>::one ), vector ) );
}

/**
 * @brief This plus operator creates an expression that represents a * x * A + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression x * A.
 * @return            The expression representing this addition.
 */
inline Expression_SVM_SV operator+( const Vector& vector, const Expression_SVM& exp )
{
    return Expression_SVM_SV( exp, Expression_SV( Scalar( Constants<IndexType>::one ), vector ) );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression a * A * y.
 * @return            The expression representing this addition.
 */
inline Expression_SMV_SV operator+( const Expression_SMV& exp, const Vector& vector )
{
    return Expression_SMV_SV( exp, Expression_SV( Scalar( Constants<IndexType>::one ), vector ) );
}

/**
 * @brief This plus operator creates an expression that represents a * x * A + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression a * y * A.
 * @return            The expression representing this addition.
 */
inline Expression_SVM_SV operator+( const Expression_SVM& exp, const Vector& vector )
{
    return Expression_SVM_SV( exp, Expression_SV( Scalar( Constants<IndexType>::one ), vector ) );
}

/* ------------------------------------------------------------------------- */
/*  Expressions 'Scalar * Matrix * Vector' +/- 'Scalar * Vector'             */
/* ------------------------------------------------------------------------- */

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp1    The expression a*A*x
 * @param[in] exp2    The expression b*y
 * @return            The expression representing this sum.
 */
inline Expression_SMV_SV operator-( const Expression_SMV& exp1, const Expression_SV& exp2 )
{
    Expression_SV minusExp2( -exp2.getArg1(), exp2.getArg2() );

    return Expression_SMV_SV( exp1, minusExp2 );
}

/**
 * @brief This minus operator creates an expression that represents a * x * A + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp1    The expression a*x*A
 * @param[in] exp2    The expression b*y
 * @return            The expression representing this sum.
 */
inline Expression_SVM_SV operator-( const Expression_SVM& exp1, const Expression_SV& exp2 )
{
    Expression_SV minusExp2( -exp2.getArg1(), exp2.getArg2() );

    return Expression_SVM_SV( exp1, minusExp2 );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp2    The expression a*x
 * @param[in] exp1    The expression b*A*y
 * @return            The expression representing this sum.
 */
inline Expression_SMV_SV operator-( Expression_SV& exp1, const Expression_SMV& exp2 )
{
    return Expression_SMV_SV( Expression_SMV( -exp2.getArg1(), exp2.getArg2() ), exp1 );
}

inline Expression_SMV_SV operator+( const Expression_SV& exp1, const Expression_SMV& exp2 )
{
    return Expression_SMV_SV( exp2, exp1 );
}

inline Expression_SMV_SV operator+( const Expression_SMV& exp1, const Expression_SV& exp2 )
{
    return Expression_SMV_SV( exp1, exp2 );
}

/**
 * @brief This minus operator creates an expression that represents a * x * A + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp2    The expression a*x
 * @param[in] exp1    The expression b*y*A
 * @return            The expression representing this sum.
 */
inline Expression_SVM_SV operator-( Expression_SV& exp1, const Expression_SVM& exp2 )
{
    return Expression_SVM_SV( Expression_SVM( -exp2.getArg1(), exp2.getArg2() ), exp1 );
}

inline Expression_SVM_SV operator+( const Expression_SV& exp1, const Expression_SVM& exp2 )
{
    return Expression_SVM_SV( exp2, exp1 );
}

inline Expression_SVM_SV operator+( const Expression_SVM& exp1, const Expression_SV& exp2 )
{
    return Expression_SVM_SV( exp1, exp2 );
}

} /* end namespace lama */

} /* end namespace scai */

