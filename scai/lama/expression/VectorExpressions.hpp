/**
 * @file VectorExpressions.hpp
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
 * @brief Operators to build symbolic expressions scalar1 * vector1 + scalar2 * vector2.
 * @author Jiri Kraus
 * @date 01.06.2011
 * @since 1.0.0
 */
#pragma once 

#include <scai/lama/Scalar.hpp>
#include <scai/lama/Vector.hpp>
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

inline Expression_SV operator*( const Scalar& alpha, const Vector& vectorX )
{
    return Expression_SV( alpha, vectorX );
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

inline Expression_SV operator*( const Vector& vectorX, const Scalar& alpha )
{
    return Expression_SV( alpha, vectorX );
}

/**
 * @brief Create a symbolic expression for the division vector / alpha
 *
 * @param[in] vector   The vector.
 * @param[in] alpha    The scalar.
 * @return             Symbolic expression [1.0/alpha] * x
 */

inline Expression_SV operator/( const Vector& vector, const Scalar& alpha )
{
    // build 1.0/ alpha as new scalar for a symbolic expression Scalar * Vector

    return Expression_SV( Scalar( 1.0 ) / alpha, vector );
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
inline Expression_SV_SV operator+( const Vector& x, const Vector& y )
{
    return Expression_SV_SV( Expression_SV( Scalar( 1.0 ), x ), Expression_SV( Scalar( 1.0 ), y ) );
}

/**
 * @brief The plus operator creates an expression that represents sum of Vector and Vector times Scalar
 *
 * @param[in] vector    The vector.
 * @param[in] exp       The Vector times Scalar.
 * @return              The expression representing this difference.
 */

inline Expression_SV_SV operator+( const Vector& vector, const Expression_SV& exp )
{
    return Expression_SV_SV( Expression_SV( Scalar( 1.0 ), vector ), exp );
}

/**
 * @brief The plus operator creates an expression that represents sum of Vector times Scalar plus Vector
 *
 * @param[in] vector    The vector.
 * @param[in] exp       Expression of Scalar times Vector.
 * @return              The expression representing this difference.
 */

inline Expression_SV_SV operator+( const Expression_SV& exp, const Vector& vector )
{
    return Expression_SV_SV( exp, Expression_SV( Scalar( 1.0 ), vector ) );
}

/**
 * @brief The plus operator creates an expression that represents sum of Vector times Scalar and Vector times Scalar
 *
 * @param[in] exp1      The vector times Scalar
 * @param[in] exp2      The vector times Scalar
 * @return              The expression representing this difference.
 */

inline Expression_SV_SV operator+( const Expression_SV& exp1, const Expression_SV& exp2 )
{
    return Expression_SV_SV( exp1, exp2 );
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

inline Expression_SV_SV operator-( const Vector& x, const Vector& y )
{
    return Expression_SV_SV( Expression_SV( Scalar( 1.0 ), x ), Expression_SV( Scalar( -1.0 ), y ) );
}

/**
 * @brief Build symbolic expression for 'Scalar * Vector' + Vector
 *
 * @param[in] exp       Expression of Scalar times Vector.
 * @param[in] vector    The vector.
 * @return              The expression representing this difference.
 *
 */

inline Expression_SV_SV operator-( const Expression_SV& exp, const Vector& vector )
{
    return Expression_SV_SV( exp, Expression_SV( -1.0, vector ) );
}

/**
 * @brief Build symbolic expression for Vector - 'Scalar * Vector'
 *
 * @param[in] vector    The vector.
 * @param[in] exp       Expression of Scalar times Vector.
 * @return              The expression representing this difference.
 */

inline Expression_SV_SV operator-( const Vector& vector, const Expression_SV& exp )
{
    Expression_SV minusExp( -exp.getArg1(), exp.getArg2() );

    return Expression_SV_SV( Expression_SV( 1.0, vector ), minusExp );
}

/**
 * @brief The minus operator creates an expression that represents the difference
 *        of two scalar vector products.
 *
 * @param[in] exp1     the minuend as symbolic 'scalar * vector'
 * @param[in] exp2     The subtranhend as symbolic 'scalar * vector'
 * @return             symoblic expression for the difference normed as 'scalar * vector + scalar * vector'
 */

inline Expression_SV_SV operator-( const Expression_SV& exp1, const Expression_SV& exp2 )
{
    Expression_SV minusExp2( -exp2.getArg1(), exp2.getArg2() );
    return Expression_SV_SV( exp1, minusExp2 );
}

} /* end namespace lama */

} /* end namespace scai */

