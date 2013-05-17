/**
 * @file Vector.cpp
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
 * @brief Vector.cpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * $Id$
 */

// hpp
#include <lama/Vector.hpp>

// others
#include <lama/DenseVector.hpp>

#include <lama/distribution/NoDistribution.hpp>

#include <lama/matrix/Matrix.hpp>

// tracing
#include <lama/tracing.hpp>

using namespace boost;

namespace lama
{

LAMA_LOG_DEF_LOGGER( Vector::logger, "Vector" )

Vector* Vector::createVector( const Scalar::ScalarType valueType, DistributionPtr distribution )
{
    switch ( valueType )
    {
    case Scalar::FLOAT:
        return new DenseVector<float>( distribution );
    case Scalar::DOUBLE:
        return new DenseVector<double>( distribution );
//        case Scalar::LONG_DOUBLE:
//            return std::auto_ptr<Vector>( new DenseVector<long double>( distribution) );
//            break;
//        case Scalar::COMPLEX:
//            return std::auto_ptr<Vector>( new DenseVector<std::complex<float> >( distribution) );
//            break;
//        case Scalar::DOUBLE_COMPLEX:
//            return std::auto_ptr<Vector>( new DenseVector<std::complex<double> >( distribution) );
//            break;
//        case Scalar::LONG_DOUBLE_COMPLEX:
//            return std::auto_ptr<Vector>( new DenseVector<std::complex<long double> >( distribution) );
//            break;
    default:
        LAMA_THROWEXCEPTION( "createVector does not support " << valueType )
    }
}

Vector::Vector( const IndexType size, ContextPtr context )
    : Distributed( shared_ptr<Distribution>( new NoDistribution( size ) ) ), mContext( context )
{
    if ( !mContext )
    {
        mContext = ContextFactory::getContext( Context::Host );
    }

    LAMA_LOG_INFO( logger, "Vector(" << size << "), replicated, on " << *mContext )
}

Vector::Vector( DistributionPtr distribution, ContextPtr context )
    : Distributed( distribution ), mContext( context )
{
    if ( !mContext )
    {
        mContext = ContextFactory::getContext( Context::Host );
    }

    LAMA_LOG_INFO( logger,
                   "Vector(" << distribution->getGlobalSize() << ") with " << getDistribution() << " constructed" )
}

Vector::Vector( const Vector& other )
    : Distributed( other ), mContext( other.getContext() )
{
    LAMA_ASSERT_ERROR( mContext, "NULL context not allowed" )
    LAMA_LOG_INFO( logger, "Vector(" << other.getDistribution().getGlobalSize() << "), distributed, copied" )
}

Vector::~Vector()
{
    LAMA_LOG_INFO( logger, "~Vector(" << getDistribution().getGlobalSize() << ")" )
}

Vector& Vector::operator=( const Expression<Matrix,Vector,Times>& expression )
{
    LAMA_LOG_DEBUG( logger, "this = matrix * vector1 -> this = 1.0 * matrix * vector1 + 0.0 * this" )

    // expression = A * x, generalized to A * x * 1.0 + 0.0 * this
    // but be careful: this might not be resized correctly, so we do it here

    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times> exp1( Scalar( 1.0 ), expression );

    const Expression<Scalar,Vector,Times> exp2( Scalar( 0.0 ), *this );

    Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> tempExpression(
        exp1, exp2 );

    // due to alias of result/vector2 resize already here

    resize( expression.getArg1().getDistributionPtr() );

    return *this = tempExpression;
}

Vector& Vector::operator=(
    const Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus>& expression )
{
    LAMA_LOG_DEBUG( logger, "this = a * vector1 + b * vector2, check vector1.size() == vector2.size()" )

    if ( expression.getArg1().getArg2().size() != expression.getArg2().getArg2().size() )
    {
        LAMA_THROWEXCEPTION(
            "size of input vector 1 " << expression.getArg1().getArg2().size() << " mismatches size of input vector 2 " << expression.getArg2().getArg2().size() )
    }

    assign( expression );
    return *this;
}

Vector& Vector::operator=( const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& expression )
{
    LAMA_LOG_INFO( logger, "this = alpha * matrix * vectorX -> this = alpha * matrix * vectorX + 0.0 * this" )

    const Scalar& beta = 0.0;

    Expression<Scalar,Vector,Times> exp2( beta, *this );

    Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> tmpExp(
        expression, exp2 );

    const Vector& vectorX = expression.getArg2().getArg2();

    if ( &vectorX != this )
    {
        // so this is not aliased to the vector on the rhs
        // as this will be used on rhs we do allocate it here
        // distribution is given by the row distribution of the matrix

        const Matrix& matrix = expression.getArg2().getArg1();

        DistributionPtr dist = matrix.getDistributionPtr();

        resize( dist );

        // values remain uninitialized as we assume that 0.0 * this (undefined) will
        // never be executed as an operation
    }

    return operator=( tmpExp );
}

Vector& Vector::operator=(
    const Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus>& expression )
{
    LAMA_LOG_INFO( logger,
                   "Vector::operator=(const Expression<Expression<Scalar, Expression<Matrix, Vector, Times>, Times>," << "Expression<Scalar, Vector, Times>, Plus>& expression)" )
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& exp1 = expression.getArg1();
    const Expression<Scalar,Vector,Times>& exp2 = expression.getArg2();
    const Scalar& alpha = exp1.getArg1();
    const Expression<Matrix,Vector,Times>& matrixTimesVectorExp = exp1.getArg2();
    const Scalar& beta = exp2.getArg1();
    const Vector& vectorY = exp2.getArg2();

    const Matrix& matrix = matrixTimesVectorExp.getArg1();
    const Vector& vectorX = matrixTimesVectorExp.getArg2();

    Vector* resultPtr = this;

    boost::shared_ptr<Vector> tmpResult;

    if ( &vectorX == this )
    {
        LAMA_LOG_DEBUG( logger, "Temporary for X required" )
        tmpResult = boost::shared_ptr<Vector>( this->create( getDistributionPtr() ) );
        resultPtr = tmpResult.get();
    }

    LAMA_LOG_DEBUG( logger, "call matrixTimesVector with matrix = " << matrix )

    matrix.matrixTimesVector( *resultPtr, alpha, vectorX, beta, vectorY );

    if ( resultPtr != this )
    {
        swap( *tmpResult );
    }

    return *this;
}

Vector& Vector::operator=( const Expression<Scalar,Vector,Times>& expression )
{
    LAMA_LOG_DEBUG( logger, "a * vector1 -> a * vector1 + 0.0 * vector1" )

    Expression<Scalar,Vector,Times> exp1( 0.0, expression.getArg2() );
    Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> tmpExp( expression, exp1 );

    // calling operator=( tmpExp ) would imply unnecessary checks, so call assign directly

    assign( tmpExp );
    return *this;
}

Vector& Vector::operator=( const Expression<Vector,Vector,Plus>& expression )
{
    LAMA_LOG_DEBUG( logger, "vector1 + vector2 -> 1.0 * vector1 + 1.0 * vector2" )

    Expression<Scalar,Vector,Times> exp1( 1.0, expression.getArg1() );
    Expression<Scalar,Vector,Times> exp2( 1.0, expression.getArg2() );
    Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> tmpExp( exp1, exp2 );

    // do not call assign( tmpExp ) here as it would skip size checks

    return operator=( tmpExp );
}

Vector& Vector::operator=( const Vector& other )
{
    Distributed::operator=( other );
    assign( other );
    return *this;
}

Vector& Vector::operator=( const Scalar value )
{
    assign( value );
    return *this;
}

Vector& Vector::operator*=( const Scalar value )
{
    Expression<Scalar, Vector, Times> exp1( value, *this );
    return operator=( exp1 );
}

Vector& Vector::operator/=( const Scalar value )
{
    Expression<Scalar, Vector, Times> exp1( Scalar( 1.0 ) / value , *this );
    return operator=( exp1 );
}

Vector& Vector::operator+=( const Vector& other )
{
    Expression<Vector, Vector, Plus> exp1( other, *this );
    return operator=( exp1 );
}

Vector& Vector::operator+=( const Expression<Scalar, Vector, Times>& expression )
{
    Expression<Scalar, Vector, Times> exp1( 1.0, *this );
    Expression<Expression<Scalar, Vector, Times>,
               Expression<Scalar, Vector, Times>,
               Plus> exp3( exp1, expression );
    return operator=( exp3 );
}

Vector& Vector::operator-=( const Expression<Scalar, Vector, Times>& expression )
{
    Expression<Scalar, Vector, Times> exp1( 1.0, *this );
    Expression<Scalar, Vector, Times> exp2( - expression.getArg1(), expression.getArg2() );
    Expression<Expression<Scalar, Vector, Times>,
               Expression<Scalar, Vector, Times>,
               Plus> exp3( exp1, exp2 );

    return operator=( exp3 );
}

Vector& Vector::operator+=( const Expression<Matrix, Vector, Times>& expression )
{
    // build exp3 = 1.0 * this + 1.0 * A * x

    Expression<Scalar, Vector, Times> exp1( 1.0, *this );
    Expression<Scalar, Expression<Matrix, Vector, Times>, Times> exp2 ( 1.0, expression );

    Expression<Expression<Scalar, Expression<Matrix, Vector, Times>, Times>,
               Expression<Scalar, Vector, Times>,
               Plus> exp3( exp2, exp1 );

    return operator=( exp3 );
}

Vector& Vector::operator-=( const Expression<Matrix, Vector, Times>& expression )
{
    // build exp3 = 1.0 * this - 1.0 * A * x

    Expression<Scalar, Vector, Times> exp1( 1.0, *this );
    Expression<Scalar, Expression<Matrix, Vector, Times>, Times> exp2 ( -1.0, expression );

    Expression<Expression<Scalar, Expression<Matrix, Vector, Times>, Times>,
               Expression<Scalar, Vector, Times>,
               Plus> exp3( exp2, exp1 );

    return operator=( exp3 );
}

Vector& Vector::operator+=( const Expression<Scalar, Expression<Matrix, Vector, Times>, Times>& expression )
{
    // build exp3 = 1.0 * this + alpha * A * x

    Expression<Scalar, Vector, Times> exp1( 1.0, *this );
    Expression<Expression<Scalar, Expression<Matrix, Vector, Times>, Times>,
               Expression<Scalar, Vector, Times>,
               Plus> exp3( expression, exp1 );

    return operator=( exp3 );
}

Vector& Vector::operator-=( const Expression<Scalar, Expression<Matrix, Vector, Times>, Times>& expression )
{
    Expression<Scalar, Vector, Times> exp1( 1.0, *this );
    Expression<Scalar, Expression<Matrix, Vector, Times>, Times> exp2( - expression.getArg1(), expression.getArg2() );
    Expression<Expression<Scalar, Expression<Matrix, Vector, Times>, Times>,
               Expression<Scalar, Vector, Times>,
               Plus> exp3( exp2, exp1 );

    return operator=( exp3 );
}

Vector& Vector::operator-=( const Vector& other )
{
    Expression<Scalar,Vector,Times> exp1( 1.0, *this );
    Expression<Scalar,Vector,Times> exp2( -1.0, other );

    Expression<Expression<Scalar,Vector,Times>,Expression<Scalar,Vector,Times>,Plus> exp3( exp1, exp2 );

    return operator=( exp3 );
}

const Scalar Vector::operator()( const IndexType i ) const
{
    return getValue( i );
}

Scalar Vector::operator*( const Vector& other ) const
{
    LAMA_REGION( "Vector.dotP" )
    return dotProduct( other );
}

void Vector::swapVector( Vector& other )
{
    // swaps only on this base class, not whole vectors

    mContext.swap( other.mContext );
    Distributed::swap( other );
}

void Vector::writeAt( std::ostream& stream ) const
{
    stream << "Vector(" << getDistributionPtr()->getGlobalSize() << ")";
}

void Vector::setContext( ContextPtr context )
{
    if ( mContext->getType() != context->getType() )
    {
        LAMA_LOG_DEBUG( logger, *this << ": new context = " << *context << ", old context = " << *mContext )
    }

    mContext = context;
}

void Vector::prefetch() const
{
    prefetch( mContext );
}

void Vector::resize( DistributionPtr distributionPtr )
{
    setDistributionPtr( distributionPtr );
    resizeImpl();
}

}

