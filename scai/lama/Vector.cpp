/**
 * @file Vector.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementations of methods for class Vector.
 * @author Jiri Kraus
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/Vector.hpp>

// local library
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/Distribution.hpp>

#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/io/PartitionIO.hpp>

// tracing
#include <scai/tracing.hpp>

// std
#include <map>
#include <ostream>

using namespace scai::common;
using namespace scai::hmemo;
using namespace scai::dmemo;

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( Vector::logger, "Vector" )

/* ---------------------------------------------------------------------------------- */

const char* _Vector::kind2Str( const VectorKind vectorKind )
{
    switch ( vectorKind )
    {
        case DENSE:
            return "DENSE";
            break;

        case SPARSE:
            return "SPARSE";
            break;

        case UNDEFINED:
            return "Undefined";
            break;
    }

    return "<illegal_vector_kind>";
}

_Vector::VectorKind _Vector::str2Kind( const char* str )

{
    for ( int kind = DENSE; kind < UNDEFINED; ++kind )
    {
        if ( strcmp( kind2Str( VectorKind( kind ) ), str ) == 0 )
        {
            return VectorKind( kind );
        }
    }

    return UNDEFINED;
}

/* ---------------------------------------------------------------------------------------*/
/*    VectorKind opertor<<                                                                */
/* ---------------------------------------------------------------------------------------*/

std::ostream& operator<<( std::ostream& stream, const _Vector::VectorKind& kind )
{
    stream << _Vector::kind2Str( kind );
    return stream;
}

/* ---------------------------------------------------------------------------------------*/
/*    Factory to create a vector                                                          */
/* ---------------------------------------------------------------------------------------*/

Vector* Vector::getVector( const VectorKind kind, const common::scalar::ScalarType valueType )
{
    VectorCreateKeyType vectype( kind, valueType );
    return Vector::create( vectype );
}

Vector* Vector::getDenseVector( 
    const common::scalar::ScalarType valueType, 
    DistributionPtr distribution,
    ContextPtr context )
{
    VectorCreateKeyType vectype( Vector::DENSE, valueType );
    Vector* v = Vector::create( vectype );
    v->allocate( distribution );

    if ( context )
    {
        v->setContextPtr( context );
    }

    return v;
}

/* ---------------------------------------------------------------------------------------*/
/*    Constructor / Destructor                                                            */
/* ---------------------------------------------------------------------------------------*/

Vector::Vector( const IndexType size, hmemo::ContextPtr context ) :

    Distributed( shared_ptr<Distribution>( new NoDistribution( size ) ) ),
    mContext( context )
{
    if ( !mContext )
    {
        mContext = Context::getHostPtr();
    }

    SCAI_LOG_INFO( logger, "Vector(" << size << "), replicated, on " << *mContext )
}

Vector::Vector( DistributionPtr distribution, hmemo::ContextPtr context )
    : Distributed( distribution ), mContext( context )
{
    if ( !mContext )
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
/*    Reading vector from a file                                                          */
/* ---------------------------------------------------------------------------------------*/

/* ---------------------------------------------------------------------------------*/

void Vector::readFromSingleFile( const std::string& fileName )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const PartitionId MASTER = 0;
    const PartitionId myRank = comm->getRank();

    // this is a bit tricky stuff, but it avoids an additional copy from array -> vector

    _HArray& localValues = const_cast<_HArray&>( getLocalValues() );

    IndexType vectorSize;

    if ( myRank == MASTER )
    {
        SCAI_LOG_INFO( logger, *comm << ": read array from single file " << fileName )

        FileIO::read( localValues, fileName );

        vectorSize = localValues.size();
    }

    comm->bcast( &vectorSize, 1, MASTER );

    if ( myRank != MASTER )
    {
        localValues.clear();
    }

    DistributionPtr distribution( new CyclicDistribution( vectorSize, vectorSize, comm ) );

    // works fine as assign can deal with alias, i.e. localValues und getLocalValuues() are same

    assign( localValues, distribution );
}

/* ---------------------------------------------------------------------------------*/

void Vector::readFromPartitionedFile( const std::string& myPartitionFileName, DistributionPtr dist )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    // this is a bit tricky stuff, but it provides us with an array of the right type
    // and it avoids an additional copy from array -> vector

    _HArray& localValues = const_cast<_HArray&>( getLocalValues() );

    bool errorFlag = false;

    IndexType localSize = 0;

    try
    {
        FileIO::read( localValues, myPartitionFileName );

        localSize = localValues.size();

        if ( dist.get() )
        {
            // size of storage must match the local size of distribution

            SCAI_ASSERT_EQUAL( dist->getLocalSize(), localSize, "serious mismatch: local matrix in " << myPartitionFileName << " has illegal local size" )
        }
    }
    catch ( common::Exception& e )
    {
        SCAI_LOG_ERROR( logger, *comm << ": failed to read " << myPartitionFileName << ": " << e.what() )
        errorFlag = true;
    }

    errorFlag = comm->any( errorFlag );

    if ( errorFlag )
    {
        COMMON_THROWEXCEPTION( "error reading partitioned matrix" )
    }

    DistributionPtr vectorDist = dist;

    if ( !vectorDist.get() )
    {
        // we have no distribution so assume a general block distribution

        IndexType globalSize = comm->sum( localSize );

        vectorDist.reset( new GenBlockDistribution( globalSize, localSize, comm ) );
    }

    // make sure that all processors have the same number of columns

    assign( localValues, vectorDist );
}

/* ---------------------------------------------------------------------------------*/

void Vector::readFromFile( const std::string& vectorFileName, const std::string& distributionFileName )
{
    // read the distribution

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr distribution = PartitionIO::readDistribution( distributionFileName, comm );

    readFromFile( vectorFileName, distribution );
}

/* ---------------------------------------------------------------------------------*/

void Vector::readFromFile( const std::string& fileName, DistributionPtr distribution )
{
    SCAI_LOG_INFO( logger,
                   *this << ": readFromFile( " << fileName << " )" )

    std::string newFileName = fileName;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();  // take default

    if ( distribution.get() )
    {
        comm = distribution->getCommunicatorPtr();
    }

    bool isPartitioned;

    PartitionIO::getPartitionFileName( newFileName, isPartitioned, *comm );

    if ( !isPartitioned )
    {
        // Alternative solution: each processor reads form the file its local part
     
        readFromSingleFile( newFileName );

        if ( distribution.get() )
        {
            redistribute( distribution );
        }
    }
    else
    {
        readFromPartitionedFile( newFileName, distribution );
    }
}

/* ---------------------------------------------------------------------------------------*/
/*    Assignment operator                                                                 */
/* ---------------------------------------------------------------------------------------*/

Vector& Vector::operator=( const Expression_MV& expression )
{
    SCAI_LOG_DEBUG( logger, "this = matrix * vector1 -> this = 1.0 * matrix * vector1 + 0.0 * this" )
    // expression = A * x, generalized to A * x * 1.0 + 0.0 * this
    // but be careful: this might not be allocated correctly, so we do it here
    const Expression_SMV exp1( Scalar( 1.0 ), expression );
    const Expression_SV exp2( Scalar( 0.0 ), *this );
    const Expression_SMV_SV tempExpression( exp1, exp2 );
    // due to alias of result/vector2 resize already here
    allocate( expression.getArg1().getRowDistributionPtr() );
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
    allocate( expression.getArg1().getDistributionPtr() );
    return *this = tempExpression;
}

Vector& Vector::operator=( const Expression_SV_SV& expression )
{
    SCAI_LOG_DEBUG( logger, "this = a * vector1 + b * vector2, check vector1.size() == vector2.size()" )
    const Vector& x = expression.getArg1().getArg2();
    const Vector& y = expression.getArg2().getArg2();
    SCAI_ASSERT_EQ_ERROR( x.size(), y.size(), "size mismatch for the two vectors in a * x + b * y" );
    assign( expression );
    return *this;
}

Vector& Vector::operator=( const Expression_SMV& expression )
{
    SCAI_LOG_INFO( logger, "this = alpha * matrix * vectorX -> this = alpha * matrix * vectorX + 0.0 * this" )
    const Scalar beta( 0.0 );
    Expression_SV exp2( beta, *this );
    Expression_SMV_SV tmpExp( expression, exp2 );
    const Vector& vectorX = expression.getArg2().getArg2();

    if ( &vectorX != this )
    {
        // so this is not aliased to the vector on the rhs
        // as this will be used on rhs we do allocate it here
        // distribution is given by the row distribution of the matrix
        const Matrix& matrix = expression.getArg2().getArg1();
        DistributionPtr dist = matrix.getRowDistributionPtr();
        allocate( dist );
        // values remain uninitialized as we assume that 0.0 * this (undefined) will
        // never be executed as an operation
    }

    return operator=( tmpExp );
}

Vector& Vector::operator=( const Expression_SVM& expression )
{
    SCAI_LOG_INFO( logger, "this = alpha * vectorX * matrix -> this = alpha * vectorX * matrix + 0.0 * this" )
    const Scalar beta( 0.0 );
    Expression_SV exp2( beta, *this );
    Expression_SVM_SV tmpExp( expression, exp2 );
    const Vector& vectorX = expression.getArg2().getArg1();

    if ( &vectorX != this )
    {
        // so this is not aliased to the vector on the rhs
        // as this will be used on rhs we do allocate it here
        // distribution is given by the row distribution of the matrix
        const Matrix& matrix = expression.getArg2().getArg2();
        DistributionPtr dist = matrix.getColDistributionPtr();
        allocate( dist );
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
    const Expression<Matrix, Vector, Times>& matrixTimesVectorExp = exp1.getArg2();
    const Scalar& beta = exp2.getArg1();
    const Vector& vectorY = exp2.getArg2();
    const Matrix& matrix = matrixTimesVectorExp.getArg1();
    const Vector& vectorX = matrixTimesVectorExp.getArg2();
    Vector* resultPtr = this;
    common::shared_ptr<Vector> tmpResult;

    if ( &vectorX == this )
    {
        SCAI_LOG_DEBUG( logger, "Temporary for X required" )
        tmpResult = common::shared_ptr<Vector>( Vector::create( this->getCreateValue() ) );
        resultPtr = tmpResult.get();
    }

    SCAI_LOG_DEBUG( logger, "call matrixTimesVector with matrix = " << matrix )
    matrix.matrixTimesVector( *resultPtr, alpha, vectorX, beta, vectorY );

    if ( resultPtr != this )
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
    const Expression<Vector, Matrix, Times>& vectorTimesMatrixExp = exp1.getArg2();
    const Scalar& beta = exp2.getArg1();
    const Vector& vectorY = exp2.getArg2();
    const Vector& vectorX = vectorTimesMatrixExp.getArg1();
    const Matrix& matrix = vectorTimesMatrixExp.getArg2();
    Vector* resultPtr = this;
    common::shared_ptr<Vector> tmpResult;

    if ( &vectorX == this )
    {
        SCAI_LOG_DEBUG( logger, "Temporary for X required" )
        tmpResult = common::shared_ptr<Vector>( Vector::create( this->getCreateValue() ) );
        resultPtr = tmpResult.get();
    }

    SCAI_LOG_DEBUG( logger, "call vectorTimesMatrix with matrix = " << matrix )
    matrix.vectorTimesMatrix( *resultPtr, alpha, vectorX, beta, vectorY );

    if ( resultPtr != this )
    {
        swap( *tmpResult );
    }

    return *this;
}

Vector& Vector::operator=( const Expression_SV& expression )
{
    SCAI_LOG_DEBUG( logger, "operator=, SV (  s * vector )  -> SV_SV ( s * vector  + 0 * vector )" )
    Expression_SV_SV tmpExp( expression, Expression_SV( Scalar( 0 ), expression.getArg2() ) );
    // calling operator=( tmpExp ) would imply unnecessary checks, so call assign directly
    assign( tmpExp );
    return *this;
}

Vector& Vector::operator=( const Expression_VV expression )
{
    SCAI_LOG_DEBUG( logger, "operator=, SVV( alpha, x, y) -> x * y" )
    Expression_SVV tmpExp( Scalar( 1.0 ), expression );
    assign( tmpExp );
    return *this;
}

Vector& Vector::operator=( const Expression_SVV expression )
{
    SCAI_LOG_DEBUG( logger, "operator=, SVV( alpha, x, y) -> alpha * x * y" )
    assign( expression );
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

Vector& Vector::operator*=( const Vector& other )
{
    return scale( other );
}

Vector& Vector::operator/=( const Scalar value )
{
    Expression<Scalar, Vector, Times> exp1( Scalar( 1.0 ) / value, *this );
    return operator=( exp1 );
}

Vector& Vector::operator+=( const Vector& other )
{
    return operator=( Expression_SV_SV( Expression_SV( Scalar( 1 ), other ), Expression_SV( Scalar( 1 ), *this ) ) );
}

Vector& Vector::operator+=( const Scalar value )
{
    add( value );
    return *this;
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

    if ( mContext->getType() != context->getType() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": new context = " << *context << ", old context = " << *mContext )
    }

    mContext = context;
}

void Vector::prefetch() const
{
    prefetch( mContext );
}

/* ---------------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
