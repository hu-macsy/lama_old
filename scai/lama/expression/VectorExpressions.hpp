/**
 * @file VectorExpressions.hpp
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
 * @brief Operators to build symbolic expressions scalar1 * vector1 + scalar2 * vector2.
 * @author Jiri Kraus
 * @date 01.06.2011
 */
#pragma once

#include <scai/lama/Scalar.hpp>
#include <scai/lama/_Vector.hpp>
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
inline Expression_SV<ValueType> operator*( const Scalar& alpha, const Vector<ValueType>& vectorX )
{
    return Expression_SV<ValueType>( alpha.getValue<ValueType>(), vectorX );
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
inline Expression_SV<ValueType> operator*( const Vector<ValueType>& vectorX, const Scalar& alpha )
{
    return Expression_SV<ValueType>( alpha.getValue<ValueType>(), vectorX );
}

/**
 * @brief Create a symbolic expression for the division vector / alpha
 *
 * @param[in] vector   The vector.
 * @param[in] alpha    The scalar.
 * @return             Symbolic expression [1.0/alpha] * x
 */
template<typename ValueType>
inline Expression_SV<ValueType> operator/( const Vector<ValueType>& vector, const Scalar& alpha )
{
    // build 1.0/ alpha as new scalar for a symbolic expression Scalar * Vector
    return Expression_SV<ValueType>( ValueType( 1 ) / alpha.getValue<ValueType>(), vector );
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

/**
 * @brief The times operator creates an expression that represents scaling a elementwise vector times
 *        vector expression.
 *
 * @param[in] alpha factor used for multiplication
 * @param[in] exp   The vector times vector expression.
 * @return          The expression representing this SVV.
 */
template<typename ValueType>
inline Expression_SVV<ValueType> operator*( const Scalar& alpha, const Expression_VV<ValueType>& exp )
{
    return Expression_SVV<ValueType>( alpha.getValue<ValueType>(), exp );
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
inline Expression_SVV<ValueType> operator*( const Expression_SV<ValueType> exp, const Vector<ValueType>& v )
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
inline Expression_SVV<ValueType> operator*( const Vector<ValueType>& v, const Expression_SV<ValueType> exp )
{
    return Expression_SVV<ValueType>( exp.getArg1(), Expression_VV<ValueType>( v, exp.getArg2() ) );
}

/* ------------------------------------------------------------------------- */
/*   operator+ to generate Expression_S_V                                    */
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
inline Expression_SV_S<ValueType> operator+( const Scalar& alpha, const Vector<ValueType>& x )
{
    return Expression_SV_S<ValueType>( Expression_SV<ValueType>( ValueType( 1 ), x ), alpha.getValue<ValueType>() );
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
inline Expression_SV_S<ValueType> operator+( const Vector<ValueType>& x, const Scalar& alpha )
{
    // return Expression_S_V( alpha, vectorX );
    return Expression_SV_S<ValueType>( Expression_SV<ValueType>( ValueType( 1 ), x ), alpha.getValue<ValueType>() );
}

/* ------------------------------------------------------------------------- */
/*   operator+ to generate Expression_SV_V                                   */
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
inline Expression_SV_S<ValueType> operator+( const Expression_SV<ValueType>& exp, const Scalar& beta )
{
    return Expression_SV_S<ValueType>( exp, beta.getValue<ValueType>() );
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
    return Expression_SV_SV<ValueType>( Expression_SV<ValueType>( ValueType( 1 ), x ), 
                                        Expression_SV<ValueType>( ValueType( 1 ), y ) );
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
    return Expression_SV_SV<ValueType>( Expression_SV<ValueType>( ValueType( 1 ), vector ), exp );
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
    return Expression_SV_SV<ValueType>( exp, Expression_SV<ValueType>( ValueType( 1 ), vector ) );
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
    return Expression_SV_SV<ValueType>( Expression_SV<ValueType>( ValueType( 1 ), x ), 
                                        Expression_SV<ValueType>( ValueType( -1 ), y ) );
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
    return Expression_SV_SV<ValueType>( exp, Expression_SV<ValueType>( ValueType( -1 ), vector ) );
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
    return Expression_SV_SV<ValueType>( Expression_SV<ValueType>( ValueType( 1 ), vector ),
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

} /* end namespace lama */

} /* end namespace scai */

