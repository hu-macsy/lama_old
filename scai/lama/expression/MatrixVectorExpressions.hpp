/**
 * @file MatrixVectorExpressions.hpp
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
 * @brief Operators to form symbolic expressions scalar * matrix * vector + scalar * vector
 * @author Thomas Brandes
 * @date 28.03.2011
 */
#pragma once

#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/Expression.hpp>

namespace scai
{

namespace lama
{

/* ------------------------------------------------------------------------- */
/*  Expressions return 'Scalar * matrix * Vector'                            */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline OpMatrix<ValueType> transpose( const Matrix<ValueType>& matrix )
{
    return OpMatrix<ValueType>( matrix, common::MatrixOp::TRANSPOSE );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of a matrix and a Vector.
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
template<typename ValueType>
inline Expression_SMV<ValueType> operator*( const Matrix<ValueType>& matrix, const Vector<ValueType>& vector )
{
    return Expression_SMV<ValueType>( intern::Scalar( 1 ), 
                                      Expression_MV<ValueType>( OpMatrix<ValueType>( matrix, common::MatrixOp::NORMAL ), 
                                                                vector ) );
}

template<typename ValueType>
inline Expression_SMV<ValueType> operator*( const OpMatrix<ValueType>& opMatrix, const Vector<ValueType>& vector )
{
    return Expression_SMV<ValueType>( intern::Scalar( 1 ), 
                                      Expression_MV<ValueType>( opMatrix, vector ) );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of a matrix, Vector and scalar.
 *
 * @param[in] matrix  The matrix.
 * @param[in] exp     The expression a*x
 * @return            The expression representing this product.
 */
template<typename ValueType>
inline Expression_SMV<ValueType> operator*( const Matrix<ValueType>& matrix, const Expression_SV<ValueType>& exp )
{
    return Expression_SMV<ValueType>( 
               exp.getArg1(), 
               Expression_MV<ValueType>( 
                   OpMatrix<ValueType>( matrix, common::MatrixOp::NORMAL ), 
                   exp.getArg2() ) );
}

template<typename ValueType>
inline Expression_SMV<ValueType> operator*( const OpMatrix<ValueType>& opMatrix, const Expression_SV<ValueType>& exp )
{
    return Expression_SMV<ValueType>( 
               exp.getArg1(), 
               Expression_MV<ValueType>( opMatrix, exp.getArg2() ) );
}

/**
 * @brief This plus operator creates an expression that represents A * B * x, where
 *        x is vector, A and B are matrices.
 *
 * @param[in] exp     The expression a*A
 * @param[in] vector  The vector.
 * @return            The expression representing this product.
 */
template<typename ValueType>
inline Expression_SMV<ValueType> operator*( const Expression_SM<ValueType>& exp, const Vector<ValueType>& vector )
{
    return Expression_SMV<ValueType>( exp.getArg1(), 
                                      Expression_MV<ValueType>( exp.getArg2(), vector ) );
}

/* ------------------------------------------------------------------------- */
/*  Expressions 'scalar * matrix * Vector' +/- vector                        */
/* ------------------------------------------------------------------------- */

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression a * A * x.
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator-( const Expression_SMV<ValueType>& exp, const Vector<ValueType>& vector )
{
    return Expression_SMV_SV<ValueType>( exp, Expression_SV<ValueType>( intern::Scalar( -1 ), vector ) );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator-( const Vector<ValueType>& vector, const Expression_SMV<ValueType>& exp )
{
    Expression_SMV<ValueType> minusExp( -exp.getArg1(), exp.getArg2() );
    return Expression_SMV_SV<ValueType>( minusExp, Expression_SV<ValueType>( intern::Scalar( 1 ), vector ) );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression A * x.
 * @return            The expression representing this addition.
 */
template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator+( const Vector<ValueType>& vector, const Expression_SMV<ValueType>& exp )
{
    return Expression_SMV_SV<ValueType>( exp, Expression_SV<ValueType>( intern::Scalar( 1 ), vector ) );
}

/**
 * @brief This plus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] vector  The input vector y.
 * @param[in] exp     The expression a * A * y.
 * @return            The expression representing this addition.
 */
template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator+( const Expression_SMV<ValueType>& exp, const Vector<ValueType>& vector )
{
    return Expression_SMV_SV<ValueType>( exp, Expression_SV<ValueType>( intern::Scalar( 1 ), vector ) );
}

/* ------------------------------------------------------------------------- */
/*  Expressions 'Scalar * matrix * Vector' +/- 'Scalar * Vector'             */
/* ------------------------------------------------------------------------- */

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp1    The expression a*A*x
 * @param[in] exp2    The expression b*y
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator-( const Expression_SMV<ValueType>& exp1, const Expression_SV<ValueType>& exp2 )
{
    Expression_SV<ValueType> minusExp2( -exp2.getArg1(), exp2.getArg2() );
    return Expression_SMV_SV<ValueType>( exp1, minusExp2 );
}

/**
 * @brief This minus operator creates an expression that represents a * A * x + b * y, where
 *        x and y are vectors, A is a matrix and a and b are scalars
 *
 * @param[in] exp2    The expression a*x
 * @param[in] exp1    The expression b*A*y
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator-( Expression_SV<ValueType>& exp1, const Expression_SMV<ValueType>& exp2 )
{
    return Expression_SMV_SV<ValueType>( Expression_SMV<ValueType>( -exp2.getArg1(), exp2.getArg2() ), exp1 );
}

template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator+( const Expression_SV<ValueType>& exp1, const Expression_SMV<ValueType>& exp2 )
{
    return Expression_SMV_SV<ValueType>( exp2, exp1 );
}

template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator+( const Expression_SMV<ValueType>& exp1, const Expression_SV<ValueType>& exp2 )
{
    return Expression_SMV_SV<ValueType>( exp1, exp2 );
}

/* ------------------------------------------------------------------------- */
/*   scaling of matrix-vector operations                                     */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline Expression_SMV<ValueType> operator*( const Expression_SMV<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SMV<ValueType>( exp.getArg1() * alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator*( const Expression_SMV_SV<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SMV_SV<ValueType>( exp.getArg1() * alpha, exp.getArg2() * alpha );
}

template<typename ValueType>
inline Expression_SMV<ValueType> operator*( const intern::Scalar alpha, const Expression_SMV<ValueType>& exp )
{
    return Expression_SMV<ValueType>( alpha * exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator*( const intern::Scalar alpha, const Expression_SMV_SV<ValueType>& exp )
{
    return Expression_SMV_SV<ValueType>( alpha * exp.getArg1(), alpha * exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMV<ValueType> operator/( const Expression_SMV<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SMV<ValueType>( exp.getArg1() / alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator/( const Expression_SMV_SV<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SMV_SV<ValueType>( exp.getArg1() / alpha, exp.getArg2() / alpha );
}

/* ------------------------------------------------------------------------- */
/*   unary operator -                                                        */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline Expression_SMV<ValueType> operator-( const Expression_SMV<ValueType>& exp )
{
    return Expression_SMV<ValueType>( -exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMV_SV<ValueType> operator-( const Expression_SMV_SV<ValueType>& exp )
{
    return Expression_SMV_SV<ValueType>( -exp.getArg1(), -exp.getArg2() );
}

} /* end namespace lama */

} /* end namespace scai */

