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
#include <scai/lama/Vector.hpp>

// local library
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/distribution/NoDistribution.hpp>

#include <scai/lama/matrix/Matrix.hpp>

// tracing
#include <scai/tracing.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <map>
#include <ostream>

using namespace scai::common;
using namespace scai::hmemo;

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( Vector::logger, "Vector" )

/* ---------------------------------------------------------------------------------------*/
/*    Factory to create a vector                                                          */
/* ---------------------------------------------------------------------------------------*/

Vector* Vector::getVector( const VectorKind kind, common::scalar::ScalarType type )
{
	using ::operator<<;

    VectorCreateKeyType key( kind, type );

    SCAI_LOG_INFO( logger, "getVector uses Factory::create " << key )

    // get it from the factory by building a pair as key the creator fn

    return create( key );
}

Vector* Vector::createVector( const common::scalar::ScalarType valueType, DistributionPtr distribution )
{
    Vector* v = getVector( DENSE, valueType );

    v->resize( distribution );
    return v;
}

/* ---------------------------------------------------------------------------------------*/
/*    Constructor / Destructor                                                            */
/* ---------------------------------------------------------------------------------------*/

Vector::Vector( const IndexType size, ContextPtr context )
                : Distributed( shared_ptr<Distribution>( new NoDistribution( size ) ) ), mContext( context )
{
    if( !mContext )
    {
        mContext = Context::getHostPtr();
    }

    SCAI_LOG_INFO( logger, "Vector(" << size << "), replicated, on " << *mContext )
}

Vector::Vector( DistributionPtr distribution, ContextPtr context )
                : Distributed( distribution ), mContext( context )
{
    if( !mContext )
    {
        mContext = Context::getHostPtr();
    }

    SCAI_LOG_INFO( logger,
                   "Vector(" << distribution->getGlobalSize() << ") with " << getDistribution() << " constructed" )
}

Vector::Vector( const Vector& other )
                : Distributed( other ), mContext( other.getContextPtr() )
{
    SCAI_ASSERT_ERROR( mContext, "NULL context not allowed" )
    SCAI_LOG_INFO( logger, "Vector(" << other.getDistribution().getGlobalSize() << "), distributed, copied" )
}

Vector::~Vector()
{
    SCAI_LOG_INFO( logger, "~Vector(" << getDistribution().getGlobalSize() << ")" )
}

/* ---------------------------------------------------------------------------------------*/
/*    Assignment operator                                                                 */
/* ---------------------------------------------------------------------------------------*/

Vector& Vector::operator=( const Expression_MV& expression )
{
    SCAI_LOG_DEBUG( logger, "this = matrix * vector1 -> this = 1.0 * matrix * vector1 + 0.0 * this" )

    // expression = A * x, generalized to A * x * 1.0 + 0.0 * this
    // but be careful: this might not be resized correctly, so we do it here

    const Expression_SMV exp1( Scalar( 1.0 ), expression );

    const Expression_SV exp2( Scalar( 0.0 ), *this );

    const Expression_SMV_SV tempExpression( exp1, exp2 );

    // due to alias of result/vector2 resize already here

    resize( expression.getArg1().getDistributionPtr() );

    return *this = tempExpression;
}

Vector& Vector::operator=( const Expression_VM& expression )
{
    SCAI_LOG_DEBUG( logger, "this = matrix * vector1 -> this = 1.0 * vector1 * matrix + 0.0 * this" )

    // expression = A * x, generalized to A * x * 1.0 + 0.0 * this
    // but be careful: this might not be resized correctly, so we do it here

    const Expression_SVM exp1( Scalar( 1.0 ), expression );

    const Expression_SV exp2( Scalar( 0.0 ), *this );

    const Expression_SVM_SV tempExpression( exp1, exp2 );

    // due to alias of result/vector2 resize already here

    resize( expression.getArg1().getDistributionPtr() );

    return *this = tempExpression;
}

Vector& Vector::operator=( const Expression_SV_SV& expression )
{
    SCAI_LOG_DEBUG( logger, "this = a * vector1 + b * vector2, check vector1.size() == vector2.size()" )

    const Vector& x = expression.getArg1().getArg2();
    const Vector& y = expression.getArg2().getArg2();

    SCAI_ASSERT_EQUAL_ERROR( x.size(), y.size() );

    assign( expression );

    return *this;
}

Vector& Vector::operator=( const Expression_SMV& expression )
{
    SCAI_LOG_INFO( logger, "this = alpha * matrix * vectorX -> this = alpha * matrix * vectorX + 0.0 * this" )

    const Scalar& beta = 0.0;

    Expression_SV exp2( beta, *this );

    Expression_SMV_SV tmpExp( expression, exp2 );

    const Vector& vectorX = expression.getArg2().getArg2();

    if( &vectorX != this )
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

Vector& Vector::operator=( const Expression_SVM& expression )
{
    SCAI_LOG_INFO( logger, "this = alpha * vectorX * matrix -> this = alpha * vectorX * matrix + 0.0 * this" )

    const Scalar& beta = 0.0;

    Expression_SV exp2( beta, *this );

    Expression_SVM_SV tmpExp( expression, exp2 );

    const Vector& vectorX = expression.getArg2().getArg1();

    if( &vectorX != this )
    {
        // so this is not aliased to the vector on the rhs
        // as this will be used on rhs we do allocate it here
        // distribution is given by the row distribution of the matrix

        const Matrix& matrix = expression.getArg2().getArg2();

        DistributionPtr dist = matrix.getColDistributionPtr();

        resize( dist );

        // values remain uninitialized as we assume that 0.0 * this (undefined) will
        // never be executed as an operation
    }

    return operator=( tmpExp );
}

Vector& Vector::operator=( const Expression_SMV_SV& expression )
{
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SMV_SV )" )

    const Expression_SMV& exp1 = expression.getArg1();
    const Expression_SV& exp2 = expression.getArg2();
    const Scalar& alpha = exp1.getArg1();
    const Expression<Matrix,Vector,Times>& matrixTimesVectorExp = exp1.getArg2();
    const Scalar& beta = exp2.getArg1();
    const Vector& vectorY = exp2.getArg2();

    const Matrix& matrix = matrixTimesVectorExp.getArg1();
    const Vector& vectorX = matrixTimesVectorExp.getArg2();

    Vector* resultPtr = this;

    common::shared_ptr<Vector> tmpResult;

    if( &vectorX == this )
    {
        SCAI_LOG_DEBUG( logger, "Temporary for X required" )
        tmpResult = common::shared_ptr<Vector>( this->clone( getDistributionPtr() ) );
        resultPtr = tmpResult.get();
    }

    SCAI_LOG_DEBUG( logger, "call matrixTimesVector with matrix = " << matrix )

    matrix.matrixTimesVector( *resultPtr, alpha, vectorX, beta, vectorY );

    if( resultPtr != this )
    {
        swap( *tmpResult );
    }

    return *this;
}

Vector& Vector::operator=( const Expression_SVM_SV& expression )
{
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SVM_SV )" )

    const Expression_SVM& exp1 = expression.getArg1();
    const Expression_SV& exp2 = expression.getArg2();
    const Scalar& alpha = exp1.getArg1();
    const Expression<Vector,Matrix,Times>& vectorTimesMatrixExp = exp1.getArg2();
    const Scalar& beta = exp2.getArg1();
    const Vector& vectorY = exp2.getArg2();

    const Vector& vectorX = vectorTimesMatrixExp.getArg1();
    const Matrix& matrix = vectorTimesMatrixExp.getArg2();

    Vector* resultPtr = this;

    common::shared_ptr<Vector> tmpResult;

    if( &vectorX == this )
    {
        SCAI_LOG_DEBUG( logger, "Temporary for X required" )
        tmpResult = common::shared_ptr<Vector>( this->clone( getDistributionPtr() ) );
        resultPtr = tmpResult.get();
    }

    SCAI_LOG_DEBUG( logger, "call vectorTimesMatrix with matrix = " << matrix )

    matrix.vectorTimesMatrix( *resultPtr, alpha, vectorX, beta, vectorY );

    if( resultPtr != this )
    {
        swap( *tmpResult );
    }

    return *this;
}

Vector& Vector::operator=( const Expression_SV& expression )
{
    SCAI_LOG_DEBUG( logger, "a * vector1 -> a * vector1 + 0.0 * vector1" )

    Expression_SV_SV tmpExp( expression, Expression_SV( Scalar( 0 ), expression.getArg2() ) );

    // calling operator=( tmpExp ) would imply unnecessary checks, so call assign directly

    assign( tmpExp );
    return *this;
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

/* ---------------------------------------------------------------------------------------*/
/*   Compound assignments                                                                 */
/* ---------------------------------------------------------------------------------------*/

Vector& Vector::operator*=( const Scalar value )
{
    return operator=( Expression_SV( value, *this ) );
}

Vector& Vector::operator/=( const Scalar value )
{
    Expression<Scalar,Vector,Times> exp1( Scalar( 1.0 ) / value, *this );
    return operator=( exp1 );
}

Vector& Vector::operator+=( const Vector& other )
{
    return operator=( Expression_SV_SV( Expression_SV( Scalar( 1 ), other ), Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator+=( const Expression_SV& exp )
{
    return operator=( Expression_SV_SV( exp, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator-=( const Expression_SV& exp )
{
    Expression_SV minusExp( -exp.getArg1(), exp.getArg2() );

    return operator=( Expression_SV_SV( minusExp, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator+=( const Expression_SMV& expression )
{
    return operator=( Expression_SMV_SV( expression, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator+=( const Expression_SVM& expression )
{
    return operator=( Expression_SVM_SV( expression, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator-=( const Expression_SMV& exp )
{
    Expression_SMV minusExp( -exp.getArg1(), exp.getArg2() );

    return operator=( Expression_SMV_SV( minusExp, Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator-=( const Vector& other )
{
    return operator=( Expression_SV_SV( Expression_SV( Scalar( 1 ), *this ), Expression_SV( Scalar( -1 ), other ) ) );
}

/* ---------------------------------------------------------------------------------------*/
/*   Miscellaneous                                                                        */
/* ---------------------------------------------------------------------------------------*/

const Scalar Vector::operator()( const IndexType i ) const
{
    return getValue( i );
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

void Vector::setContextPtr( ContextPtr context )
{
    SCAI_ASSERT_DEBUG( context, "NULL context invalid" )

    if( mContext->getType() != context->getType() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": new context = " << *context << ", old context = " << *mContext )
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

} /* end namespace lama */

} /* end namespace scai */
