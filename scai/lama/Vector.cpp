/**
 * @file Vector.cpp
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
    const intern::Scalar alphaS        = expression.getArg1().getArg1();
    const intern::Scalar betaS         = expression.getArg2().getArg1();
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
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_S_SV<ValueType>& expression )
{
    const Expression_SV<ValueType>& exp = expression.getArg2();
    const intern::Scalar& alpha = exp.getArg1();
    const Vector<ValueType>& x = exp.getArg2();
    const intern::Scalar& beta = expression.getArg1();

    vectorPlusScalar( alpha.getValue<ValueType>(), x, beta.getValue<ValueType>() );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SMV<ValueType>& expression )
{
    intern::Scalar alphaS = expression.getArg1(); 
    ValueType alpha = alphaS.getValue<ValueType>();

    const OpMatrix<ValueType>& opMat = expression.getArg2().getArg1();
    const Vector<ValueType>& vector = expression.getArg2().getArg2();

    const Matrix<ValueType>& matrix = opMat.getMatrix();
    const common::MatrixOp op = opMat.getOp();

    SCAI_LOG_INFO( logger, "this = " << alpha << " * matrix " << &matrix << " * vector " << &vector )

    matrix.matrixTimesVector( *this, alpha, vector, ValueType( 0 ), nullptr, op );

    return *this;
}

template<>
Vector<IndexType>& Vector<IndexType>::operator=( const Expression_SMV<IndexType>& )
{
    COMMON_THROWEXCEPTION( "Matrix<IndexType> not supported" )
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SMV_SV<ValueType>& expression )
{
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SMV_SV )" )

    const Expression_SMV<ValueType>& exp1 = expression.getArg1();
    const Expression_SV<ValueType>& exp2  = expression.getArg2();

    // split up exp1 -> alpha * op( matrix ) * vectorX

    const intern::Scalar& alphaS = exp1.getArg1();
    ValueType alpha = alphaS.getValue<ValueType>();
    const Expression_MV<ValueType> matrixTimesVectorExp = exp1.getArg2();
    const common::MatrixOp op = matrixTimesVectorExp.getArg1().getOp();
    const Matrix<ValueType>& matrix = matrixTimesVectorExp.getArg1().getMatrix();
    const Vector<ValueType>& vectorX = matrixTimesVectorExp.getArg2();

    // split up exp2 -> beta * vectorY

    const intern::Scalar& betaS = exp2.getArg1();
    const ValueType& beta = betaS.getValue<ValueType>();
    const Vector<ValueType>& vectorY = exp2.getArg2();

    matrix.matrixTimesVector( *this, alpha, vectorX, beta, &vectorY, op );

    return *this;
}

template<>
Vector<IndexType>& Vector<IndexType>::operator=( const Expression_SMV_SV<IndexType>& )
{
    COMMON_THROWEXCEPTION( "Matrix<IndexType> not supported" )
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SV<ValueType>& expression )
{
    const intern::Scalar& alpha = expression.getArg1();
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

    intern::Scalar a = exp.getArg1();

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
    const intern::Scalar b = exp.getArg1();

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
    intern::Scalar beta  = exp.getArg1();

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
/*   fill assembly data                                                                   */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void Vector<ValueType>::fillFromAssembly( const VectorAssembly<ValueType>& assembly, common::BinaryOp op )
{
    hmemo::HArray<IndexType> nonZeroIndexes;
    hmemo::HArray<ValueType> nonZeroValues;
 
    const dmemo::Distribution& dist = getDistribution();
    const dmemo::Communicator& comm = dist.getCommunicator();

    if ( comm == assembly.getCommunicator() )
    {
         // assembly collected on same set of processors

         assembly.buildLocalData( nonZeroIndexes, nonZeroValues, getDistribution() );
    }
    else if ( comm.getType() == dmemo::Communicator::NO )
    {
         // we can build a replicated vector from assembled data 

         assembly.buildGlobalData( nonZeroIndexes, nonZeroValues, dist.getGlobalSize() );
    }
    else
    {
         COMMON_THROWEXCEPTION( "Processor set of assembled data : " << assembly.getCommunicator()
                                << " does not match processor set of distribution " << dist )
    }

    fillSparseData( nonZeroIndexes, nonZeroValues, op );
}

template<typename ValueType>
void Vector<ValueType>::disassemble( VectorAssembly<ValueType>& assembly, const IndexType offset ) const
{
    const dmemo::Distribution& dist = getDistribution();
    const dmemo::Communicator& comm = dist.getCommunicator();

    if ( comm == assembly.getCommunicator() )
    {
        // that is fine, either replicated or distributed
    } 
    else if ( comm.getType() == dmemo::Communicator::NO )
    {
        // disassemble of replicated matrix will only be done first processor

        if ( comm.getRank() != 0 )
        {
            return;
        }
    }
    else 
    {
        COMMON_THROWEXCEPTION( "Processor set of assembly = " << assembly.getCommunicator() 
                             << " does not match to processor set onto which vector is distributed: " << comm )
    }

    hmemo::HArray<IndexType> ownedIndexes;
    hmemo::HArray<ValueType> localData;

    getDistribution().getOwnedIndexes( ownedIndexes );
    buildLocalValues( localData );

    auto values       = hostReadAccess( localData );
    auto local2Global = hostReadAccess( ownedIndexes );

    ValueType zero = 0;

    for ( IndexType i = 0; i < localData.size(); ++i )
    {
        if ( values[i] == zero )
        {
            continue;  // skip zero values 

        }

        // Be careful: we need global indexes in assembly

        assembly.push( offset + local2Global[i], values[i] );
    }
}


/* ---------------------------------------------------------------------------------------*/
/*   concatenation of vectors                                                             */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void Vector<ValueType>::concatenate( dmemo::DistributionPtr dist, const std::vector<const Vector<ValueType>*>& vectors )
{
    ValueType zero = 0;

    // use an assembly to collect local values from any distribution

    VectorAssembly<ValueType> assembly( dist->getCommunicatorPtr() );

    IndexType offset = 0;

    for ( size_t k = 0; k < vectors.size(); ++k )
    {
        const Vector<ValueType>& v = *vectors[k];

        if ( offset + v.size() > dist->getGlobalSize() )
        {
            COMMON_THROWEXCEPTION( "concatenate fails, exceeds global size of target vector" )
        }

        v.disassemble( assembly, offset );

        offset += v.size();
    }

    setSameValue( dist, zero );
    fillFromAssembly( assembly );
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
