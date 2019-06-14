/**
 * @file VectorExpressions.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Operators to build symbolic expressions scalar1 * vector1 + scalar2 * vector2.
 * @author Jiri Kraus
 * @date 01.06.2011
 */
#pragma once

#include <scai/lama/Scalar.hpp>
#include <scai/lama/expression/Expression.hpp>

namespace scai
{

namespace lama
{

/**
 * @brief Create a symbolic expression for the product alpha * vectorX
 *
 * @param[in] alpha     The scalar.
 * @param[in] vectorX   The vector.
 * @return              Symbolic expression alpha * vectorX
 */

template<typename ValueType>
inline Expression_SV<ValueType> operator*( const intern::Scalar& alpha, const Vector<ValueType>& vectorX )
{
    return Expression_SV<ValueType>( alpha, vectorX );
}

/**
 * @brief Create a symbolic expression for the product vectorX * alpha
 *
 * @param[in] vectorX   The vector.
 * @param[in] alpha     The scalar.
 * @return              Symbolic expression alpha * vectorX
 *
 * Note: due to normalization the arguments are switched in the symbolic expression
 */

template<typename ValueType>
inline Expression_SV<ValueType> operator*( const Vector<ValueType>& vectorX, const intern::Scalar& alpha )
{
    return Expression_SV<ValueType>( alpha, vectorX );
}

/**
 * @brief Create a symbolic expression for the division vector / alpha
 *
 * @param[in] vector   The vector.
 * @param[in] alpha    The scalar.
 * @return             Symbolic expression [1.0/alpha] * x
 */
template<typename ValueType>
inline Expression_SV<ValueType> operator/( const Vector<ValueType>& vector, const intern::Scalar& alpha )
{
    // build 1.0/ alpha as new scalar for a symbolic expression Scalar * Vector
    return Expression_SV<ValueType>( intern::Scalar( 1 ) / alpha, vector );
}

/**
 * @brief The plus operator creates an expression that represents the sum
 *        of two vectors.
 *
 * @param[in] x     The first vector.
 * @param[in] y     The second vector.
 * @return          The expression representing this sum.
 */
template<typename ValueType>
inline Expression_VV<ValueType> operator*( const Vector<ValueType>& x, const Vector<ValueType>& y )
{
    return Expression_VV<ValueType>( x, y );
}

/* ------------------------------------------------------------------------- */
/*   operator* to generate Expression_SVV                                    */
/* ------------------------------------------------------------------------- */

/**
 * @brief The times operator creates an expression that represents scaling a elementwise vector times
 *        vector expression.
 *
 * @param[in] alpha factor used for multiplication
 * @param[in] exp   The vector times vector expression.
 * @return          The expression representing this SVV.
 */
template<typename ValueType>
inline Expression_SVV<ValueType> operator*( const intern::Scalar& alpha, const Expression_VV<ValueType>& exp )
{
    return Expression_SVV<ValueType>( alpha, exp );
}
 
/** @brief Multiplication vector * vector with a scalar on the rhs */

template<typename ValueType>
inline Expression_SVV<ValueType> operator*( const Expression_VV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_SVV<ValueType>( alpha, exp );
}

template<typename ValueType>
inline Expression_SVV<ValueType> operator/( const Expression_VV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_SVV<ValueType>( 1 / intern::Scalar( alpha ), exp );
}

/**
 * @brief The times operator creates an expression that represents scaling a elementwise vector times
 *        vector expression.
 *
 * @param[in] exp   an existing expression scalar * vector
 * @param[in] v     the second multiplicator
 * @return          The expression representing this SVV.
 */
template<typename ValueType>
inline Expression_SVV<ValueType> operator*( const Expression_SV<ValueType>& exp, const Vector<ValueType>& v )
{
    return Expression_SVV<ValueType>( exp.getArg1(), Expression_VV<ValueType>( exp.getArg2(), v ) );
}

/**
 * @brief The times operator creates an expression that represents scaling a elementwise vector times
 *        vector expression.
 *
 * @param[in] v     the first multiplicator
 * @param[in] exp   an existing expression scalar * vector
 * @return          The expression representing this SVV.
 */
template<typename ValueType>
inline Expression_SVV<ValueType> operator*( const Vector<ValueType>& v, const Expression_SV<ValueType>& exp )
{
    return Expression_SVV<ValueType>( exp.getArg1(), Expression_VV<ValueType>( v, exp.getArg2() ) );
}

template<typename ValueType>
inline Expression_SVV<ValueType> operator*( const Expression_SV<ValueType>& exp1, const Expression_SV<ValueType>& exp2 )
{
    return Expression_SVV<ValueType>( exp1.getArg1() * exp2.getArg1(), 
                                      Expression_VV<ValueType>( exp1.getArg2(), exp2.getArg2() ) );
}

/* ------------------------------------------------------------------------- */
/*   operator+ to generate Expression_S_SV                                   */
/* ------------------------------------------------------------------------- */

/**
 * @brief The plus operator creates an expression that represents the sum
 *        a vector adding a scalar to each vector element.
 *
 * @param[in] alpha is the scalar value to be added.
 * @param[in] x     The vector.
 * @return          The expression representing this sum.
 */

template<typename ValueType>
inline Expression_S_SV<ValueType> operator+( const intern::Scalar& alpha, const Vector<ValueType>& x )
{
    return Expression_S_SV<ValueType>( alpha, Expression_SV<ValueType>( intern::Scalar( 1 ), x ) );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator-( const intern::Scalar& alpha, const Vector<ValueType>& x )
{
    return Expression_S_SV<ValueType>( alpha, Expression_SV<ValueType>( intern::Scalar( -1 ), x ) );
}

/**
 * @brief The plus operator creates an expression that represents the sum
 *        a vector adding a scalar to each vector element.
 *
 * @param[in] x     The vector.
 * @param[in] alpha The scalar.
 * @return          The expression representing this sum.
 */

template<typename ValueType>
inline Expression_S_SV<ValueType> operator+( const Vector<ValueType>& x, const intern::Scalar& alpha )
{
    // return Expression_S_V( alpha, vectorX );
    return Expression_S_SV<ValueType>( alpha, Expression_SV<ValueType>( intern::Scalar( 1 ), x ) );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator-( const Vector<ValueType>& x, const intern::Scalar& alpha )
{
    // return Expression_S_V( alpha, vectorX );
    return Expression_S_SV<ValueType>( -alpha, Expression_SV<ValueType>( intern::Scalar( 1 ), x ) );
}

/* ------------------------------------------------------------------------- */
/*   operator+, operator- to generate Expression_SV_V                        */
/* ------------------------------------------------------------------------- */

/**
 * @brief The plus operator creates an expression that represents the sum
 *        a vector adding a scalar to each vector element.
 *
 * @param[in] exp   an existing expression scalar * vector
 * @param[in] beta  The scalar.
 * @return          The expression representing this sum.
 */

template<typename ValueType>
inline Expression_S_SV<ValueType> operator+( const Expression_SV<ValueType>& exp, const intern::Scalar& beta )
{
    return Expression_S_SV<ValueType>( beta, exp );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator+( const intern::Scalar& beta, const Expression_SV<ValueType>& exp )
{
    return Expression_S_SV<ValueType>( beta, exp );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator-( const Expression_SV<ValueType>& exp, const intern::Scalar& beta )
{
    return Expression_S_SV<ValueType>( -beta, exp );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator-( const intern::Scalar& beta, const Expression_SV<ValueType>& exp )
{
    return Expression_S_SV<ValueType>( beta, -exp );
}

/* ------------------------------------------------------------------------- */
/*   S_SV +- scalar, scalar +- S_SV remains S_SV                             */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline Expression_S_SV<ValueType> operator+( const intern::Scalar& beta, const Expression_S_SV<ValueType>& exp )
{
    return Expression_S_SV<ValueType>( beta + exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator+( const Expression_S_SV<ValueType>& exp, const intern::Scalar& beta )
{
    return Expression_S_SV<ValueType>( beta + exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator-( const intern::Scalar& beta, const Expression_S_SV<ValueType>& exp )
{
    return Expression_S_SV<ValueType>( beta - exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator-( const Expression_S_SV<ValueType>& exp, const intern::Scalar& beta )
{
    return Expression_S_SV<ValueType>( exp.getArg1(), exp.getArg2() - beta );
}

/* ------------------------------------------------------------------------- */
/*   operator+ to generate Expression_SV_SV                                  */
/* ------------------------------------------------------------------------- */

/**
 * @brief The plus operator creates an expression that represents the sum
 *        of two vectors.
 *
 * @param[in] x     The first vector.
 * @param[in] y     The second vector.
 * @return          The expression representing this sum.
 */
template<typename ValueType>
inline Expression_SV_SV<ValueType> operator+( const Vector<ValueType>& x, const Vector<ValueType>& y )
{
    return Expression_SV_SV<ValueType>( Expression_SV<ValueType>( intern::Scalar( 1 ), x ), 
                                        Expression_SV<ValueType>( intern::Scalar( 1 ), y ) );
}

/**
 * @brief The plus operator creates an expression that represents sum of Vector and Vector times Scalar
 *
 * @param[in] vector    The vector.
 * @param[in] exp       The Vector times Scalar.
 * @return              The expression representing this difference.
 */
template<typename ValueType>
inline Expression_SV_SV<ValueType> operator+( const Vector<ValueType>& vector, const Expression_SV<ValueType>& exp )
{
    return Expression_SV_SV<ValueType>( Expression_SV<ValueType>( intern::Scalar( 1 ), vector ), exp );
}

/**
 * @brief The plus operator creates an expression that represents sum of Vector times Scalar plus Vector
 *
 * @param[in] vector    The vector.
 * @param[in] exp       Expression of Scalar times Vector.
 * @return              The expression representing this difference.
 */
template<typename ValueType>
inline Expression_SV_SV<ValueType> operator+( const Expression_SV<ValueType>& exp, const Vector<ValueType>& vector )
{
    return Expression_SV_SV<ValueType>( exp, Expression_SV<ValueType>( intern::Scalar( 1 ), vector ) );
}

/**
 * @brief The plus operator creates an expression that represents sum of vector times scalar plus vector times scalar
 *
 * @param[in] exp1      The vector times Scalar
 * @param[in] exp2      The vector times Scalar
 * @return              The expression representing this difference.
 */
template<typename ValueType>
inline Expression_SV_SV<ValueType> operator+( const Expression_SV<ValueType>& exp1, const Expression_SV<ValueType>& exp2 )
{
    return Expression_SV_SV<ValueType>( exp1, exp2 );
}

/* ------------------------------------------------------------------------- */
/*   operator- to generate Expression_SV_SV                                  */
/* ------------------------------------------------------------------------- */

/**
 * @brief The minus operator creates an expression that represents the difference
 *        of two vectors.
 *
 * @param[in] x     The first vector as minuend
 * @param[in] y     The second vector is subtrahend
 * @return          Normalized symbolic expression ' 1 * x + (-1) * y'
 */
template<typename ValueType>
inline Expression_SV_SV<ValueType> operator-( const Vector<ValueType>& x, const Vector<ValueType>& y )
{
    return Expression_SV_SV<ValueType>( Expression_SV<ValueType>( intern::Scalar( 1 ), x ), 
                                        Expression_SV<ValueType>( intern::Scalar( -1 ), y ) );
}

/**
 * @brief Build symbolic expression for 'Scalar * Vector' + Vector
 *
 * @param[in] exp       Expression of Scalar times Vector.
 * @param[in] vector    The vector.
 * @return              The expression representing this difference.
 *
 */
template<typename ValueType>
inline Expression_SV_SV<ValueType> operator-( const Expression_SV<ValueType>& exp, const Vector<ValueType>& vector )
{
    return Expression_SV_SV<ValueType>( exp, Expression_SV<ValueType>( intern::Scalar( -1 ), vector ) );
}

/**
 * @brief Build symbolic expression for Vector - 'Scalar * Vector'
 *
 * @param[in] vector    The vector.
 * @param[in] exp       Expression of Scalar times Vector.
 * @return              The expression representing this difference.
 */
template<typename ValueType>
inline Expression_SV_SV<ValueType> operator-( const Vector<ValueType>& vector, const Expression_SV<ValueType>& exp )
{
    return Expression_SV_SV<ValueType>( Expression_SV<ValueType>( intern::Scalar( 1 ), vector ),
                                        Expression_SV<ValueType>( -exp.getArg1(), exp.getArg2() ) );
}

/**
 * @brief The minus operator creates an expression that represents the difference
 *        of two scalar vector products.
 *
 * @param[in] exp1     the minuend as symbolic 'scalar * vector'
 * @param[in] exp2     The subtranhend as symbolic 'scalar * vector'
 * @return             symoblic expression for the difference normed as 'scalar * vector + scalar * vector'
 */
template<typename ValueType>
inline Expression_SV_SV<ValueType> operator-( const Expression_SV<ValueType>& exp1, const Expression_SV<ValueType>& exp2 )
{
    Expression_SV<ValueType> minusExp2( -exp2.getArg1(), exp2.getArg2() );
    return Expression_SV_SV<ValueType>( exp1, minusExp2 );
}

/* ------------------------------------------------------------------------- */
/*   operator* to scale supported expressions                                */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline Expression_SV<ValueType> operator*( const intern::Scalar& alpha, const Expression_SV<ValueType>& exp )
{
    return Expression_SV<ValueType>( alpha * exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SV<ValueType> operator*( const Expression_SV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_SV<ValueType>( exp.getArg1() * alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_SVV<ValueType> operator*( const intern::Scalar& alpha, const Expression_SVV<ValueType>& exp )
{
    return Expression_SVV<ValueType>( alpha * exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SVV<ValueType> operator*( const Expression_SVV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_SVV<ValueType>( exp.getArg1() * alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator*( const intern::Scalar& alpha, const Expression_S_SV<ValueType>& exp )
{
    return Expression_S_SV<ValueType>( alpha * exp.getArg1(), alpha * exp.getArg2() );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator*( const Expression_S_SV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_S_SV<ValueType>( exp.getArg1() * alpha, exp.getArg2() * alpha );
}

template<typename ValueType>
inline Expression_SV_SV<ValueType> operator*( const intern::Scalar& alpha, const Expression_SV_SV<ValueType>& exp )
{
    return Expression_SV_SV<ValueType>( alpha * exp.getArg1(), alpha * exp.getArg2() );
}

template<typename ValueType>
inline Expression_SV_SV<ValueType> operator*( const Expression_SV_SV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_SV_SV<ValueType>( exp.getArg1() * alpha, exp.getArg2() * alpha );
}

/* ------------------------------------------------------------------------- */
/*   operator/ ( exp, s ) is same as operator*( exp, 1/s )                   */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline Expression_SV<ValueType> operator/( const Expression_SV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_SV<ValueType>( exp.getArg1() / alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_SVV<ValueType> operator/( const Expression_SVV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_SVV<ValueType>( exp.getArg1() / alpha, exp.getArg2() );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator/( const Expression_S_SV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_S_SV<ValueType>( exp.getArg1() / alpha, exp.getArg2() / alpha );
}

template<typename ValueType>
inline Expression_SV_SV<ValueType> operator/( const Expression_SV_SV<ValueType>& exp, const intern::Scalar& alpha )
{
    return Expression_SV_SV<ValueType>( exp.getArg1() / alpha, exp.getArg2() / alpha );
}

/* ------------------------------------------------------------------------- */
/*   unary operator -, similiar to scale with -1                             */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
inline Expression_SV<ValueType> operator-( const Vector<ValueType>& v )
{
    return Expression_SV<ValueType>( intern::Scalar( -1 ), v );
}

template<typename ValueType>
inline Expression_SV<ValueType> operator-( const Expression_SV<ValueType>& exp )
{
    return Expression_SV<ValueType>( -exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_SVV<ValueType> operator-( const Expression_SVV<ValueType>& exp )
{
    return Expression_SVV<ValueType>( -exp.getArg1(), exp.getArg2() );
}

template<typename ValueType>
inline Expression_S_SV<ValueType> operator-( const Expression_S_SV<ValueType>& exp )
{
    return Expression_S_SV<ValueType>( -exp.getArg1(), -exp.getArg2() );
}

template<typename ValueType>
inline Expression_SV_SV<ValueType> operator-( const Expression_SV_SV<ValueType>& exp )
{
    return Expression_SV_SV<ValueType>( -exp.getArg1(), -exp.getArg2() );
}

/* ------------------------------------------------------------------------- */
/*   Failing operator's to give useful error messages                        */
/* ------------------------------------------------------------------------- */

// Note: using of delete more convenient but does not give a useful error message

template<typename ValueType, typename X>
inline Expression_SV_SV<ValueType> operator+( const Expression_SV_SV<ValueType>& exp, const X& any ) 
{
    static_assert( sizeof( X ) == -1, "too complex expression, introduce temporary" );
    return exp;
}

template<typename ValueType, typename X>
inline Expression_SV_SV<ValueType> operator+( const X& any, const Expression_SV_SV<ValueType>& exp )
{
    static_assert( sizeof( X ) == -1, "too complex expression, introduce temporary" );
    return exp;
}

/* ========================================================================= */
/*   other binary operators                                                  */
/* ========================================================================= */

/*  min( alpha, x ), min( x, alpha ), min( x, y ) */

template<typename ValueType>
inline Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MIN> min( const intern::Scalar& alpha, const Vector<ValueType>& x )
{
    return Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MIN>( alpha, x );
}

template<typename ValueType>
inline Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MIN> min( const Vector<ValueType>& x, const intern::Scalar& alpha )
{
    return Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MIN>( alpha, x );
}

template<typename ValueType>
inline Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::MIN> 
    min( const Vector<ValueType>& x, const Vector<ValueType>& y )
{
    return Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::MIN>( x, y );
}

/*  max( alpha, x ), max( x, alpha ), max( x, y ) */

template<typename ValueType>
inline Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MAX> max( const intern::Scalar& alpha, const Vector<ValueType>& x )
{
    return Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MAX>( alpha, x );
}

template<typename ValueType>
inline Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MAX> max( const Vector<ValueType>& x, const intern::Scalar& alpha )
{
    return Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::MAX>( alpha, x );
}

template<typename ValueType>
inline Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::MAX> 
    max( const Vector<ValueType>& x, const Vector<ValueType>& y )
{
    return Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::MAX>( x, y );
}

/*  pow( alpha, x ), pow( x, alpha ), pow( x, y ) */

template<typename ValueType>
inline Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::POW> pow( const intern::Scalar& alpha, const Vector<ValueType>& x )
{
    return Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::POW>( alpha, x );
}

template<typename ValueType>
inline Expression<Vector<ValueType>, intern::Scalar, common::BinaryOp::POW> pow( const Vector<ValueType>& x, const intern::Scalar& alpha )
{
    return Expression<Vector<ValueType>, intern::Scalar, common::BinaryOp::POW>( x, alpha );
}

template<typename ValueType>
inline Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::POW> 
    pow( const Vector<ValueType>& x, const Vector<ValueType>& y )
{
    return Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::POW>( x, y );
}

/*  alpha / x, x / y,  but NOT here: x / alpha is handled as x * ( 1/alpha ) */

template<typename ValueType>
inline Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::DIVIDE> operator/( const intern::Scalar& alpha, const Vector<ValueType>& x )
{
    return Expression<intern::Scalar, Vector<ValueType>, common::BinaryOp::DIVIDE>( alpha, x );
}

template<typename ValueType>
inline Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::DIVIDE> operator/( 
    const Vector<ValueType>& x, const Vector<ValueType>& y )
{
    return Expression<Vector<ValueType>, Vector<ValueType>, common::BinaryOp::DIVIDE>( x, y );
}

} /* end namespace lama */

} /* end namespace scai */

