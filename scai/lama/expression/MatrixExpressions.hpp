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

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/Scalar.hpp>

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
inline Expression_SM_SM operator+( const Matrix& matrixA, const Matrix& matrixB )
{
    return Expression_SM_SM( Expression_SM( Scalar( 1.0 ), matrixA ), Expression_SM( Scalar( 1.0 ), matrixB ) );
}

/**
 * @brief Make symbolic expression for the difference of two matrices
 *
 * @param[in] matrixA The first matrix.
 * @param[in] matrixB The second matrix.
 * @return            symbolic expression 1.0 * matrixA - 1.0 * matrixB
 */
inline Expression_SM_SM operator-( const Matrix& matrixA, const Matrix& matrixB )
{
    return Expression_SM_SM( Expression_SM( Scalar( 1.0 ), matrixA ), Expression_SM( Scalar( -1.0 ), matrixB ) );
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
inline Expression_SM operator*( const Scalar& scalar, const Matrix& matrix )
{
    return Expression_SM( scalar, matrix );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix.
 *
 * @param[in] matrix  The input matrix.
 * @param[in] scalar  The input scalar.
 * @return            The expression representing this product.
 */
inline Expression_SM operator*( const Matrix& matrix, const Scalar& scalar )
{
    return Expression_SM( scalar, matrix );
}

/**
 * @brief Create a symbolic expression 'Scalar * Matrix' for the division matrix / alpha
 *
 * @param[in] matrix   The matrix.
 * @param[in] alpha    The scalar.
 * @return             symbolic expression [1.0/alpha] *  matrixA
 */

inline Expression_SM operator/( const Matrix& matrix, const Scalar& alpha )
{
    // build 1.0/ alpha as new scalar for a symbolic expression Scalar * Matrix
    return Expression_SM( Scalar( 1.0 ) / alpha, matrix );
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
inline Expression_SMM operator*( const Matrix& m1, const Matrix& m2 )
{
    return Expression_SMM( Expression_SM( Scalar( 1.0 ), m1 ), m2 );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix, alpha * matrixA * matrixB
 *
 * @param[in] matrix  The first input matrix.
 * @param[in] exp     The expression scalar times matrix.
 * @return            The expression representing this product.
 */
inline Expression_SMM operator*( const Matrix& matrix, const Expression_SM& exp )
{
    return Expression_SMM( Expression_SM( exp.getArg1(), matrix ), exp.getArg2() );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix.
 *
 * @param[in] exp     The expression scalar times matrix.
 * @param[in] matrix  The first input matrix.
 * @return            The expression representing this product.
 */
inline Expression_SMM operator*( const Expression_SM& exp, const Matrix& matrix )
{
    return Expression_SMM( exp, matrix );
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   'Scalar * Matrix * Matrix'  * Scalar            */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix.
 *
 * @param[in] scalar  The Scalar
 * @param[in] exp     The expression 'Scalar * Matrix * Matrix'
 * @return            The expression representing this product.
 */
inline Expression_SMM operator*( const Scalar& scalar, const Expression_SMM& exp )
{
    const Expression_SM sm = exp.getArg1();
    return Expression_SMM( Expression_SM( scalar * sm.getArg1(), sm.getArg2() ), exp.getArg2() );
}

/**
 * @brief This times operator creates an expression that represents the product
 *        of Scalar times Matrix times Matrix.
 *
 * @param[in] exp     The expression scalar times matrix.
 * @param[in] scalar  The Scalar
 * @return            The expression representing this product.
 */
inline Expression_SMM operator*( const Expression_SMM& exp, const Scalar& scalar )
{
    const Expression_SM sm = exp.getArg1();
    return Expression_SMM( Expression_SM( sm.getArg1() * scalar, sm.getArg2() ), exp.getArg2() );
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
inline Expression_SMM_SM operator+( const Expression_SMM& exp1, const Expression_SM& exp2 )
{
    return Expression_SMM_SM( exp1, exp2 );
}

/**
 * @brief This plus operator creates an expression that represents the sum
 *        of Scalar times Matrix times Matrix plus Scalar times Matrix.
 *
 * @param[in] exp2    The expression Scalar times Matrix.
 * @param[in] exp1    The expression Scalar times Matrix times Matrix.
 * @return            The expression representing this sum.
 */
inline Expression_SMM_SM operator+( const Expression_SM& exp2, const Expression_SMM& exp1 )
{
    return Expression_SMM_SM( exp1, exp2 );
}

/**
 * @brief This times operator creates an expression that represents the sum
 *        of Scalar*(Matrix*Matrix) + Scalar*Matrix.
 *
 * @param[in] exp1    The expression Scalar*(Matrix*Matrix).
 * @param[in] exp2    The expression Scalar*Matrix.
 * @return            The expression representing this sum.
 */
inline Expression_SMM_SM operator-( const Expression_SMM& exp1, const Expression_SM& exp2 )
{
    Expression_SM minusExp2( -exp2.getArg1(), exp2.getArg2() );
    return Expression_SMM_SM( exp1, minusExp2 );
}

/**
 * @brief This times operator creates an expression that represents the sum
 *        of Scalar*(Matrix*Matrix) + Scalar*Matrix.
 *
 * @param[in] exp2    The expression Scalar*Matrix.
 * @param[in] exp1    The expression Scalar*(Matrix*Matrix).
 * @return            The expression representing this sum.
 */
inline Expression_SMM_SM operator-( const Expression_SM& exp2, const Expression_SMM& exp1 )
{
    Expression_SM exp1SM = exp1.getArg1();
    Expression_SMM minusExp1( Expression_SM( -exp1SM.getArg1(), exp1SM.getArg2() ), exp1.getArg2() );
    return Expression_SMM_SM( minusExp1, exp2 );
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
inline Expression_SMM_SM operator+( const Expression_SMM& exp, const Matrix& matrix )
{
    return Expression_SMM_SM( exp, Expression_SM( Scalar( 1.0 ), matrix ) );
}

/**
 * @brief Build symbolic sum of 'Matrix' + 'Scalar * Matrix * Matrix'
 *
 * @param[in] matrix  Matrix summand
 * @param[in] exp     The expression Scalar times Matrix times Matrix.
 * @return            The expression representing this sum.
 */
inline Expression_SMM_SM operator+( const Matrix& matrix, const Expression_SMM& exp )
{
    return Expression_SMM_SM( exp, Expression_SM( Scalar( 1.0 ), matrix ) );
}

/**
 * @brief Build symbolic difference of 'Scalar * Matrix * Matrix' - 'Matrix'
 *
 * @param[in] exp     The expression Scalar times Matrix times Matrix.
 * @param[in] matrix  Matrix as subtrahend
 * @return            The expression representing this sum.
 */
inline Expression_SMM_SM operator-( const Expression_SMM& exp, const Matrix& matrix )
{
    return Expression_SMM_SM( exp, Expression_SM( Scalar( -1.0 ), matrix ) );
}

/**
 * @brief Build symbolic difference of 'Scalar * Matrix * Matrix' - 'Matrix'
 *
 * @param[in] matrix  Matrix as minuend
 * @param[in] exp     The expression Scalar times Matrix times Matrix.
 * @return            The expression representing this sum.
 */
inline Expression_SMM_SM operator-( const Matrix& matrix, const Expression_SMM& exp )
{
    // Build temporary expression for -exp
    Expression_SM expSM = exp.getArg1();
    Expression_SMM minusExp( Expression_SM( -expSM.getArg1(), expSM.getArg2() ), exp.getArg2() );
    return Expression_SMM_SM( minusExp, Expression_SM( Scalar( 1.0 ), matrix ) );
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

inline Expression_SM_SM operator+( const Expression_SM& exp1, const Expression_SM& exp2 )
{
    return Expression_SM_SM( exp1, exp2 );
}

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        difference of two symbolic expressions 'Scalar * Matrix' - 'Scalar * Matrix'
 *
 * @param[in] exp1      symbolic expression alpha * matrixA
 * @param[in] exp2      symbolic expression beta * matrixB
 * @return              Symbolic expression for the difference of the two expressions
 */

inline Expression_SM_SM operator-( const Expression_SM& exp1, const Expression_SM& exp2 )
{
    Scalar minusBeta = -exp2.getArg1();
    Expression_SM minusExp2( -exp2.getArg1(), exp2.getArg2() );
    return Expression_SM_SM( exp1, minusExp2 );
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*  expressions:   Scalar * Matrix +/- Matrix                      */
/*                                                                 */
/* --------------------------------------------------------------- */

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        sum of matrix + 'Scalar * Matrix'
 *
 * @param[in] matrix    first summand
 * @param[in] exp       symbolic expression 'Scalar * Matrix'
 * @return              Symbolic expression for the sum of the two expressions
 */

inline Expression_SM_SM operator+( const Matrix& matrix, const Expression_SM& exp )
{
    return Expression_SM_SM( Expression_SM( Scalar( 1.0 ), matrix ), exp );
}

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        difference of matrix - 'Scalar * Matrix'
 *
 * @param[in] matrix    first summand
 * @param[in] exp       symbolic expression 'Scalar * Matrix'
 * @return              Symbolic expression for the difference of the two expressions
 */

inline Expression_SM_SM operator-( const Matrix& matrix, const Expression_SM& exp )
{
    Expression_SM minusExp( -exp.getArg1(), exp.getArg2() );
    return Expression_SM_SM( Expression_SM( Scalar( 1.0 ), matrix ), minusExp );
}

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        sum of 'Scalar * Matrix' + matrix
 *
 * @param[in] exp       symbolic expression 'Scalar * Matrix'
 * @param[in] matrix    second summand
 * @return              Symbolic expression for the sum of the two expressions
 */

inline Expression_SM_SM operator+( const Expression_SM& exp, const Matrix& matrix )
{
    return Expression_SM_SM( exp, Expression_SM( Scalar( 1.0 ), matrix ) );
}

/**
 * @brief Make a symbolic expression 'Scalar * Matrix + Scalar * Matrix' for the
 *        difference 'Scalar * Matrix' - matrix
 *
 * @param[in] exp       symbolic expression 'Scalar * Matrix' as minuend^
 * @param[in] matrix    subtrahend
 * @return              Symbolic expression for the difference of the two expressions
 */

inline Expression_SM_SM operator-( const Expression_SM& exp, const Matrix& matrix )
{
    return Expression_SM_SM( exp, Expression_SM( Scalar( -1.0 ), matrix ) );
}

} /* end namespace lama */

} /* end namespace scai */
