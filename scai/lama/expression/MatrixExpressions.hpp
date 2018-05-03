/**
 * @file MatrixExpressions.hpp
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
 * @brief Operators to build symbolic matrix expressions for alpha * matrix + beta * matrix.
 * @author brandes
 * @date 28.03.2011
 */
#pragma once

#include <scai/lama/expression/Expression.hpp>

namespace scai
{

namespace lama
{

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   Matrix +/- Matrix                               */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief Expression object that stands for the sum of two matrices
 *
 * @param[in] matrixA The first matrix.
 * @param[in] matrixB The second matrix.
 * @return            symbolic expression 1.0 * matrixA + 1.0 * matrixB
 */
template<typename ValueType>
inline Expression_SM_SM<ValueType> operator+( const Matrix<ValueType>& matrixA, const Matrix<ValueType>& matrixB )
{
    auto op = common::MatrixOp::NORMAL;

    return Expression_SM_SM<ValueType>( Expression_SM<ValueType>( intern::Scalar( 1 ), OpMatrix<ValueType>( matrixA, op ) ), 
                                        Expression_SM<ValueType>( intern::Scalar( 1 ), OpMatrix<ValueType>( matrixB, op ) ) );
}

/**
 * @brief Make symbolic expression for the difference of two matrices
 *
 * @param[in] matrixA The first matrix.
 * @param[in] matrixB The second matrix.
 * @return            symbolic expression 1.0 * matrixA - 1.0 * matrixB
 */
template<typename ValueType>
inline Expression_SM_SM<ValueType> operator-( const Matrix<ValueType>& matrixA, const Matrix<ValueType>& matrixB )
{
    auto op = common::MatrixOp::NORMAL;

    return Expression_SM_SM<ValueType>( Expression_SM<ValueType>( intern::Scalar( 1 ), OpMatrix<ValueType>( matrixA, op ) ), 
                                        Expression_SM<ValueType>( intern::Scalar( -1 ), OpMatrix<ValueType>( matrixB, op ) ) );
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   Matrix * Scalar                                 */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix.
 *
 * @param[in] scalar  The input scalar.
 * @param[in] matrix  The input matrix.
 * @return            The expression representing this product.
 */
template<typename ValueType>
inline Expression_SM<ValueType> operator*( const intern::Scalar& scalar, const Matrix<ValueType>& matrix )
{
    return Expression_SM<ValueType>( scalar, OpMatrix<ValueType>( matrix, common::MatrixOp::NORMAL ) );
}

template<typename ValueType>
inline Expression_SM<ValueType> operator*( const intern::Scalar& scalar, const OpMatrix<ValueType>& opMatrix )
{
    return Expression_SM<ValueType>( scalar, opMatrix );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix.
 *
 * @param[in] matrix  The input matrix.
 * @param[in] scalar  The input scalar.
 * @return            The expression representing this product.
 */
template<typename ValueType>
inline Expression_SM<ValueType> operator*( const Matrix<ValueType>& matrix, const intern::Scalar& scalar )
{
    return Expression_SM<ValueType>( scalar, OpMatrix<ValueType>( matrix, common::MatrixOp::NORMAL ) );
}

template<typename ValueType>
inline Expression_SM<ValueType> operator*( const OpMatrix<ValueType>& opMatrix, const intern::Scalar& scalar )
{
    return Expression_SM<ValueType>( scalar, opMatrix );
}

/**
 * @brief Create a symbolic expression 'Scalar * Matrix' if the - operator is applied to a matrix.
 */
template<typename ValueType>
inline Expression_SM<ValueType> operator-( const Matrix<ValueType>& matrix )
{
    return Expression_SM<ValueType>( intern::Scalar( -1 ), OpMatrix<ValueType>( matrix, common::MatrixOp::NORMAL ) );
}

template<typename ValueType>
inline Expression_SM<ValueType> operator-( const OpMatrix<ValueType>& opMatrix )
{
    return Expression_SM<ValueType>( intern::Scalar( -1 ), opMatrix );
}

/**
 * @brief Create a symbolic expression 'Scalar * Matrix' for the division matrix / alpha
 *
 * @param[in] matrix   The matrix.
 * @param[in] alpha    The scalar.
 * @return             symbolic expression [1.0/alpha] *  matrixA
 */

template<typename ValueType>
inline Expression_SM<ValueType> operator/( const Matrix<ValueType>& matrix, const intern::Scalar& alpha )
{
    // build 1.0/ alpha as new scalar for a symbolic expression Scalar * Matrix
    return Expression_SM<ValueType>( intern::Scalar( 1 ) / alpha, OpMatrix<ValueType>( matrix, common::MatrixOp::NORMAL ) );
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   '[Scalar *] Matrix '  * Matrix                  */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief This times operator creates an expression that represents the product
 *        of a Matrix and a Matrix.
 *
 * @param[in] m1      The first input matrix.
 * @param[in] m2      The second input matrix.
 * @return            Symbolic expression representing the product 'm1 * m2'
 */
template<typename ValueType>
inline Expression_SMM<ValueType> operator*( const Matrix<ValueType>& m1, const Matrix<ValueType>& m2 )
{
    auto op = common::MatrixOp::NORMAL;
    return Expression_SMM<ValueType>( Expression_SM<ValueType>( intern::Scalar( 1 ), OpMatrix<ValueType>( m1, op ) ),  OpMatrix<ValueType>( m2, op ) );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix, alpha * matrixA * matrixB
 *
 * @param[in] matrix  The first input matrix.
 * @param[in] exp     The expression scalar times matrix.
 * @return            The expression representing this product.
 */
template<typename ValueType>
inline Expression_SMM<ValueType> operator*( const Matrix<ValueType>& matrix, const Expression_SM<ValueType>& exp )
{
    auto op = common::MatrixOp::NORMAL;
    return Expression_SMM<ValueType>( Expression_SM<ValueType>( exp.getArg1(), OpMatrix<ValueType>( matrix, op ) ), exp.getArg2() );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix.
 *
 * @param[in] exp     The expression scalar times matrix.
 * @param[in] matrix  The first input matrix.
 * @return            The expression representing this product.
 */
template<typename ValueType>
inline Expression_SMM<ValueType> operator*( const Expression_SM<ValueType>& exp, const Matrix<ValueType>& matrix )
{
    return Expression_SMM<ValueType>( exp, OpMatrix<ValueType>( matrix, common::MatrixOp::NORMAL ) );
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   Scalar * Matrix * Matrix +/- Scalar * Matrix    */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief Build symbolic sum of 'Scalar * Matrix * Matrix' + 'Scalar * Matrix'
 *
 * @param[in] exp1    The expression Scalar times Matrix times Matrix.
 * @param[in] exp2    The expression Scalar times Matrix.
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator+( const Expression_SMM<ValueType>& exp1, const Expression_SM<ValueType>& exp2 )
{
    return Expression_SMM_SM<ValueType>( exp1, exp2 );
}

/**
 * @brief This plus operator creates an expression that represents the sum
 *        of Scalar times Matrix times Matrix plus Scalar times Matrix.
 *
 * @param[in] exp2    The expression Scalar times Matrix.
 * @param[in] exp1    The expression Scalar times Matrix times Matrix.
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator+( const Expression_SM<ValueType>& exp2, const Expression_SMM<ValueType>& exp1 )
{
    return Expression_SMM_SM<ValueType>( exp1, exp2 );
}

/**
 * @brief This times operator creates an expression that represents the sum
 *        of Scalar*(Matrix*Matrix) + Scalar*Matrix.
 *
 * @param[in] exp1    The expression Scalar*(Matrix*Matrix).
 * @param[in] exp2    The expression Scalar*Matrix.
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator-( const Expression_SMM<ValueType>& exp1, const Expression_SM<ValueType>& exp2 )
{
    Expression_SM<ValueType> minusExp2( -exp2.getArg1(), exp2.getArg2() );
    return Expression_SMM_SM<ValueType>( exp1, minusExp2 );
}

/**
 * @brief This times operator creates an expression that represents the sum
 *        of Scalar*(Matrix*Matrix) + Scalar*Matrix.
 *
 * @param[in] exp2    The expression Scalar*Matrix.
 * @param[in] exp1    The expression Scalar*(Matrix*Matrix).
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator-( const Expression_SM<ValueType>& exp2, const Expression_SMM<ValueType>& exp1 )
{
    Expression_SM<ValueType> exp1SM = exp1.getArg1();
    Expression_SMM<ValueType> minusExp1( Expression_SM<ValueType>( -exp1SM.getArg1(), exp1SM.getArg2() ), exp1.getArg2() );
    return Expression_SMM_SM<ValueType>( minusExp1, exp2 );
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   Scalar * Matrix * Matrix +/- Matrix             */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief Build symbolic sum of 'Scalar * Matrix * Matrix' + 'Matrix'
 *
 * @param[in] exp     The expression Scalar times Matrix times Matrix.
 * @param[in] matrix  Matrix summand
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator+( const Expression_SMM<ValueType>& exp, const Matrix<ValueType>& matrix )
{
    auto op = common::MatrixOp::NORMAL;   // identity on matrix

    return Expression_SMM_SM<ValueType>( exp, Expression_SM<ValueType>( intern::Scalar( 1 ), OpMatrix<ValueType>( matrix, op ) ) );
}

/**
 * @brief Build symbolic sum of 'Matrix' + 'scalar * Matrix * Matrix'
 *
 * @param[in] matrix  Matrix summand
 * @param[in] exp     The expression scalar times Matrix times Matrix.
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator+( const Matrix<ValueType>& matrix, const Expression_SMM<ValueType>& exp )
{
    auto op = common::MatrixOp::NORMAL;   // identity on matrix

    return Expression_SMM_SM<ValueType>( exp, Expression_SM<ValueType>( intern::Scalar( 1 ), OpMatrix<ValueType>( matrix, op ) ) );
}

/**
 * @brief Build symbolic difference of 'Scalar * Matrix * Matrix' - 'Matrix'
 *
 * @param[in] exp     The expression Scalar times Matrix times Matrix.
 * @param[in] matrix  Matrix as subtrahend
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator-( const Expression_SMM<ValueType>& exp, const Matrix<ValueType>& matrix )
{
    auto op = common::MatrixOp::NORMAL;   // identity on matrix

    return Expression_SMM_SM<ValueType>( exp, Expression_SM<ValueType>( intern::Scalar( -1 ), OpMatrix<ValueType>( matrix, op ) ) );
}

/**
 * @brief Build symbolic difference of 'scalar * Matrix * Matrix' - 'Matrix'
 *
 * @param[in] matrix  Matrix as minuend
 * @param[in] exp     The expression Scalar times Matrix times Matrix.
 * @return            The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator-( const Matrix<ValueType>& matrix, const Expression_SMM<ValueType>& exp )
{
    auto op = common::MatrixOp::NORMAL;   // identity on matrix

    // Build temporary expression for -exp
    Expression_SM<ValueType> expSM = exp.getArg1();
    Expression_SMM<ValueType> minusExp( Expression_SM<ValueType>( -expSM.getArg1(), expSM.getArg2() ), exp.getArg2() );
    return Expression_SMM_SM<ValueType>( minusExp, Expression_SM<ValueType>( intern::Scalar( 1 ), OpMatrix<ValueType>( matrix, op ) ) );
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   Scalar * Matrix +/- Scalar * Matrix             */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        sum of two symbolic expressions 'Scalar * Matrix' + 'Scalar * Matrix'
 *
 * @param[in] exp1      The vector times Scalar
 * @param[in] exp2      The vector times Scalar
 * @return              Symbolic expression for the sum of the two expressions
 */
template<typename ValueType>
inline Expression_SM_SM<ValueType> operator+( const Expression_SM<ValueType>& exp1, const Expression_SM<ValueType>& exp2 )
{
    return Expression_SM_SM<ValueType>( exp1, exp2 );
}

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        difference of two symbolic expressions 'Scalar * Matrix' - 'scalar * Matrix'
 *
 * @param[in] exp1      symbolic expression alpha * matrixA
 * @param[in] exp2      symbolic expression beta * matrixB
 * @return              Symbolic expression for the difference of the two expressions
 */
template<typename ValueType>
inline Expression_SM_SM<ValueType> operator-( const Expression_SM<ValueType>& exp1, const Expression_SM<ValueType>& exp2 )
{
    ValueType minusBeta = -exp2.getArg1();
    Expression_SM<ValueType> minusExp2( -exp2.getArg1(), exp2.getArg2() );
    return Expression_SM_SM<ValueType>( exp1, minusExp2 );
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   scalar * Matrix +/- Matrix                      */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief Make a symbolic expression 'scalar * Matrix + scalar * Matrix' for the
 *        sum of matrix + 'scalar * Matrix'
 *
 * @param[in] matrix    first summand
 * @param[in] exp       symbolic expression 'scalar * Matrix'
 * @return              Symbolic expression for the sum of the two expressions
 */
template<typename ValueType>
inline Expression_SM_SM<ValueType> operator+( const Matrix<ValueType>& matrix, const Expression_SM<ValueType>& exp )
{
    return Expression_SM_SM<ValueType>( Expression_SM<ValueType>( intern::Scalar( 1 ), matrix ), exp );
}

/**
 * @brief Make a symbolic expression 'scalar * Matrix + scalar * Matrix' for the
 *        difference of matrix - 'scalar * Matrix'
 *
 * @param[in] matrix    first summand
 * @param[in] exp       symbolic expression 'scalar * Matrix'
 * @return              Symbolic expression for the difference of the two expressions
 */

template<typename ValueType>
inline Expression_SM_SM<ValueType> operator-( const Matrix<ValueType>& matrix, const Expression_SM<ValueType>& exp )
{
    Expression_SM<ValueType> minusExp( -exp.getArg1(), exp.getArg2() );
    return Expression_SM_SM<ValueType>( Expression_SM<ValueType>( intern::Scalar( 1 ), matrix ), minusExp );
}

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        sum of 'Scalar * Matrix' + matrix
 *
 * @param[in] exp       symbolic expression 'Scalar * Matrix'
 * @param[in] matrix    second summand
 * @return              Symbolic expression for the sum of the two expressions
 */
template<typename ValueType>
inline Expression_SM_SM<ValueType> operator+( const Expression_SM<ValueType>& exp, const Matrix<ValueType>& matrix )
{
    return Expression_SM_SM<ValueType>( exp, Expression_SM<ValueType>( intern::Scalar( 1 ), matrix ) );
}

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        difference 'Scalar * Matrix' - matrix
 *
 * @param[in] exp       symbolic expression 'Scalar * Matrix' as minuend^
 * @param[in] matrix    subtrahend
 * @return              Symbolic expression for the difference of the two expressions
 */
template<typename ValueType>
inline Expression_SM_SM<ValueType> operator-( const Expression_SM<ValueType>& exp, const Matrix<ValueType>& matrix )
{
    return Expression_SM_SM<ValueType>( exp, Expression_SM<ValueType>( intern::Scalar( -1 ), matrix ) );
}

/* ------------------------------------------------------------------------- */
/*   Scaling of matrix expressions                                           */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline Expression_SM<ValueType> operator*( const intern::Scalar alpha, const Expression_SM<ValueType>& exp )
{
    return Expression_SM<ValueType>( alpha * exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMM<ValueType> operator*( const intern::Scalar alpha, const Expression_SMM<ValueType>& exp )
{
    return Expression_SMM<ValueType>( alpha * exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SM_SM<ValueType> operator*( const intern::Scalar alpha, const Expression_SM_SM<ValueType>& exp )
{
    return Expression_SM_SM<ValueType>( alpha * exp.getArg1(), alpha * exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator*( const intern::Scalar alpha, const Expression_SMM_SM<ValueType>& exp )
{
    return Expression_SMM_SM<ValueType>( alpha * exp.getArg1(), alpha * exp.getArg2() );
}

template<typename ValueType>
inline Expression_SM<ValueType> operator*( const Expression_SM<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SM<ValueType>( exp.getArg1() * alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMM<ValueType> operator*( const Expression_SMM<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SMM<ValueType>( exp.getArg1() * alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_SM_SM<ValueType> operator*( const Expression_SM_SM<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SM_SM<ValueType>( exp.getArg1() * alpha, exp.getArg2() * alpha );
}

template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator*( const Expression_SMM_SM<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SMM_SM<ValueType>( exp.getArg1() * alpha, exp.getArg2() * alpha );
}

template<typename ValueType>
inline Expression_SM<ValueType> operator/( const Expression_SM<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SM<ValueType>( exp.getArg1() / alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMM<ValueType> operator/( const Expression_SMM<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SMM<ValueType>( exp.getArg1() / alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_SM_SM<ValueType> operator/( const Expression_SM_SM<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SM_SM<ValueType>( exp.getArg1() / alpha, exp.getArg2() / alpha );
}

template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator/( const Expression_SMM_SM<ValueType>& exp, const intern::Scalar alpha )
{
    return Expression_SMM_SM<ValueType>( exp.getArg1() / alpha, exp.getArg2() / alpha );
}

/* ------------------------------------------------------------------------- */
/*   unary operator -  same as scaling with -1                               */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline Expression_SM<ValueType> operator-( const Expression_SM<ValueType>& exp )
{
    return Expression_SM<ValueType>( -exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMM<ValueType> operator-( const Expression_SMM<ValueType>& exp )
{
    return Expression_SMM<ValueType>( -exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SM_SM<ValueType> operator-( const Expression_SM_SM<ValueType>& exp )
{
    return Expression_SM_SM<ValueType>( -exp.getArg1(), -exp.getArg2() );
}

template<typename ValueType>
inline Expression_SMM_SM<ValueType> operator-( const Expression_SMM_SM<ValueType>& exp )
{
    return Expression_SMM_SM<ValueType>( -exp.getArg1(), -exp.getArg2() );
}

} /* end namespace lama */

} /* end namespace scai */
