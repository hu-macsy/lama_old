/**
 * @file Vector.cpp
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
 * @brief Implementation of methods for the abstract class Vector.
 * @author Thomas Brandes
 * @date 31.10.2017
 */

#include <scai/lama/Vector.hpp>
#include <scai/lama/matrix/Matrix.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using common::TypeTraits;

namespace lama
{

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>* Vector<ValueType>::getVector( VectorKind kind )
{
    // create it by factory, untyped vector is cast to typed vector.

    return static_cast<Vector<ValueType>*>( _Vector::getVector( kind, TypeTraits<ValueType>::stype ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>::Vector( const IndexType size, hmemo::ContextPtr context ) :
       
   _Vector( size, context )
    
{
    // context will be set by base class _Vector

    SCAI_LOG_DEBUG( logger, "Vector<" << TypeTraits<ValueType>::id() 
                           << ">( size = " << size << ", ctx = " << getContext() << " )" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>::Vector( dmemo::DistributionPtr distribution, hmemo::ContextPtr context ) :

    _Vector( distribution, context )

{
    // context will be set by base class _Vector

    SCAI_LOG_DEBUG( logger, "Vector<" << TypeTraits<ValueType>::id()
                           << ">( dist = " << getDistribution() << ", ctx = " << getContext() << " )" )
}

template<typename ValueType>
Vector<ValueType>::Vector( const _Vector& other ) : _Vector( other )
{
}

template<typename ValueType>
Vector<ValueType>::Vector( const Vector<ValueType>& other ) : _Vector( other )
{
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Vector<ValueType>::~Vector()
{
    SCAI_LOG_DEBUG( logger, "~Vector<" << TypeTraits<ValueType>::id() << ">" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::ScalarType Vector<ValueType>::getValueType() const
{
    return TypeTraits<ValueType>::stype;
}

/* ========================================================================= */
/*        operator= < vector expression>                                     */
/* ========================================================================= */

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SV_SV<ValueType>& expression )
{
    const Scalar alphaS        = expression.getArg1().getArg1();
    const Scalar betaS         = expression.getArg2().getArg1();
    const ValueType alpha      = alphaS.getValue<ValueType>();
    const ValueType beta       = betaS.getValue<ValueType>();
    const Vector<ValueType>& x = expression.getArg1().getArg2();
    const Vector<ValueType>& y = expression.getArg2().getArg2();

    SCAI_LOG_DEBUG( logger, "this = " << alpha << " * x = " << x << " + " << beta << " * y = " << y )

    // Note: all checks are done the vector specific implementations

    vectorPlusVector( alpha, x, beta, y );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SV_S<ValueType>& expression )
{
    const Expression_SV<ValueType>& exp = expression.getArg1();
    const Scalar& alpha = exp.getArg1();
    const Vector<ValueType>& x = exp.getArg2();
    const Scalar& beta = expression.getArg2();

    vectorPlusScalar( alpha.getValue<ValueType>(), x, beta.getValue<ValueType>() );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SMV<ValueType>& expression )
{
    Scalar alphaS = expression.getArg1(); 
    ValueType alpha = alphaS.getValue<ValueType>();

    const Matrix<ValueType>& matrix = expression.getArg2().getArg1();
    const Vector<ValueType>& vector = expression.getArg2().getArg2();

    SCAI_LOG_INFO( logger, "this = " << alpha << " * matrix * vector" )

    matrix.matrixTimesVector( *this, alpha, vector, ValueType( 0 ), *this, false );

    return *this;
}

template<>
Vector<IndexType>& Vector<IndexType>::operator=( const Expression_SMV<IndexType>& )
{
    COMMON_THROWEXCEPTION( "Matrix<IndexType> not supported" )
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SVM<ValueType>& expression )
{   
    Scalar alphaS = expression.getArg1(); 
    ValueType alpha = alphaS.getValue<ValueType>();

    const Vector<ValueType>& vector = expression.getArg2().getArg1();
    const Matrix<ValueType>& matrix = expression.getArg2().getArg2();

    SCAI_LOG_INFO( logger, "this = " << alpha << " * vector * matrix" )

    matrix.matrixTimesVector( *this, alpha, vector, ValueType( 0 ), *this, true );

    return *this;
}

template<>
Vector<IndexType>& Vector<IndexType>::operator=( const Expression_SVM<IndexType>& )
{
    COMMON_THROWEXCEPTION( "Matrix<IndexType> not supported" )
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SMV_SV<ValueType>& expression )
{
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SMV_SV )" )
    const Expression_SMV<ValueType>& exp1 = expression.getArg1();
    const Expression_SV<ValueType>& exp2 = expression.getArg2();
    Scalar alphaS = exp1.getArg1();
    ValueType alpha = alphaS.getValue<ValueType>();
    const Expression_MV<ValueType> matrixTimesVectorExp = exp1.getArg2();
    const Scalar& betaS = exp2.getArg1();
    const ValueType& beta = betaS.getValue<ValueType>();
    const Vector<ValueType>& vectorY = exp2.getArg2();
    const Matrix<ValueType>& matrix = matrixTimesVectorExp.getArg1();
    const Vector<ValueType>& vectorX = matrixTimesVectorExp.getArg2();

    matrix.matrixTimesVector( *this, alpha, vectorX, beta, vectorY, false );

    return *this;
}

template<>
Vector<IndexType>& Vector<IndexType>::operator=( const Expression_SMV_SV<IndexType>& )
{
    COMMON_THROWEXCEPTION( "Matrix<IndexType> not supported" )
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SVM_SV<ValueType>& expression )

{   
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SVM_SV )" )
    const Expression_SVM<ValueType>& exp1 = expression.getArg1();
    const Expression_SV<ValueType>& exp2 = expression.getArg2();
    const Expression_VM<ValueType>& vectorTimesMatrixExp = exp1.getArg2();
    
    // resolve : result = alhpa * A * x + beta * y
    
    const Scalar& alphaS  = exp1.getArg1();
    const Scalar& betaS   = exp2.getArg1();
    const ValueType alpha = alphaS.getValue<ValueType>();
    const ValueType beta  = betaS.getValue<ValueType>();
    const Vector<ValueType>& vectorY = exp2.getArg2();
    const Vector<ValueType>& vectorX = vectorTimesMatrixExp.getArg1();
    const Matrix<ValueType>& matrix = vectorTimesMatrixExp.getArg2();
    
    matrix.matrixTimesVector( *this, alpha, vectorX, beta, vectorY, true );
    
    return *this;
}

template<>
Vector<IndexType>& Vector<IndexType>::operator=( const Expression_SVM_SV<IndexType>& )
{
    COMMON_THROWEXCEPTION( "Matrix<IndexType> not supported" )
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SV<ValueType>& expression )
{
    const Scalar& alpha = expression.getArg1();
    const Vector<ValueType>& x = expression.getArg2();

    vectorPlusVector( alpha.getValue<ValueType>(), x, 0, x );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_VV<ValueType>& expression )
{
    SCAI_LOG_DEBUG( logger, "operator=, SVV( alpha, x, y) -> x * y" )

    const Vector<ValueType>& x = expression.getArg1();
    const Vector<ValueType>& y = expression.getArg2();

    ValueType alpha = 1;

    vectorTimesVector( alpha, x, y );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SVV<ValueType>& exp )
{
    // extract componennts from alpha * ( x * y )

    Scalar a = exp.getArg1();

    const ValueType alpha = a.getValue<ValueType>();

    const Vector<ValueType>& x = exp.getArg2().getArg1();
    const Vector<ValueType>& y = exp.getArg2().getArg2();

    vectorTimesVector( alpha, x, y );

    return *this;
}

/* ---------------------------------------------------------------------------------------*/
/*   vector [?]= scalar                                                                   */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const ValueType value )
{
    setScalar( value );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator+=( const ValueType value )
{
    binaryOp( *this, common::BinaryOp::ADD, value );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator-=( const ValueType value )
{
    binaryOp( *this, common::BinaryOp::SUB, value );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator*=( const ValueType value )
{
    binaryOp( *this, common::BinaryOp::MULT, value );
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator/=( const ValueType value )
{
    SCAI_ASSERT_NE_ERROR( value, ValueType( 0 ), "Divide by zero for vector" )

    // Note: multiplication is faster than division, so do it right here

    binaryOp( *this, common::BinaryOp::MULT, ValueType( 1 ) / value );

    return *this;
}

/* ---------------------------------------------------------------------------------------*/
/*   vector [?]= scalar * vector                                                          */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator+=( const Expression_SV<ValueType>& exp )
{
    const Scalar b = exp.getArg1();

    ValueType alpha = 1;
    ValueType beta  = b.getValue<ValueType>();

    const Vector<ValueType>& x = *this;
    const Vector<ValueType>& y = exp.getArg2();

    vectorPlusVector( alpha, x, beta, y );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator-=( const Expression_SV<ValueType>& exp )
{
    ValueType alpha = 1;
    Scalar beta  = exp.getArg1();

    const Vector<ValueType>& x = *this;
    const Vector<ValueType>& y = exp.getArg2();

    vectorPlusVector( alpha, x, -beta.getValue<ValueType>(), y );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator+=( const Expression_SMV<ValueType>& expression )
{
    return operator=( Expression_SMV_SV<ValueType>( expression, Expression_SV<ValueType>( ValueType( 1 ), *this ) ) );
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator+=( const Expression_SVM<ValueType>& expression )
{
    return operator=( Expression_SVM_SV<ValueType>( expression, Expression_SV<ValueType>( ValueType( 1 ), *this ) ) );
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator-=( const Expression_SMV<ValueType>& exp )
{
    Expression_SMV<ValueType> minusExp( -exp.getArg1(), exp.getArg2() );
    return operator=( Expression_SMV_SV<ValueType>( minusExp, Expression_SV<ValueType>( ValueType( 1 ), *this ) ) );
}

/* ---------------------------------------------------------------------------------------*/
/*   setRandom, setSparseRandom                                                           */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void Vector<ValueType>::setSparseRandom( const IndexType n, const ValueType& zeroValue, const float fillRate, const IndexType bound )
{
    allocate( n );

    if ( fillRate < 1.0f )
    {
        setScalar( zeroValue );
        fillSparseRandom( fillRate, bound );
    }
    else
    {
        // initialization with zero value not required
        fillRandom( bound );
    }
}

template<typename ValueType>
void Vector<ValueType>::setSparseRandom( dmemo::DistributionPtr dist, const ValueType& zeroValue, const float fillRate, const IndexType bound )
{
    allocate( dist );

    if ( fillRate < 1.0f )
    {
        setScalar( zeroValue );
        fillSparseRandom( fillRate, bound );
    }
    else
    {
        // initialization with zero value not required
        fillRandom( bound );
    }
}

/* ---------------------------------------------------------------------------------------*/
/*   assign concatenation of vectors                                                      */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void Vector<ValueType>::cat( const Vector<ValueType>& v1, const Vector<ValueType>& v2 )
{
    std::vector<const Vector<ValueType>*> vectors;

    vectors.push_back( &v1 );
    vectors.push_back( &v2 );

    dmemo::CommunicatorPtr comm = v1.getDistribution().getCommunicatorPtr();

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( v1.size() + v2.size(), comm ) );

    SCAI_LOG_INFO( logger, "this = " << *this << ", dist of concat vector = " << *dist )

    concatenate( dist, vectors );
}

/* ---------------------------------------------------------------------------------------*/
/*   miscallaneous                                                                        */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void Vector<ValueType>::setRandom( const IndexType n, const IndexType bound )
{
    allocate ( n );
    fillRandom( bound );
}

template<typename ValueType>
void Vector<ValueType>::setRandom( dmemo::DistributionPtr dist, const IndexType bound )
{
    allocate ( dist );
    fillRandom( bound );
}


/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Vector, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
