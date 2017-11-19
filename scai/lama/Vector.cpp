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
    return reinterpret_cast<Vector<ValueType>*>( _Vector::getVector( kind, TypeTraits<ValueType>::stype ) );
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

/* ---------------------------------------------------------------------------------------*/
/*   element-wise operations on vector                                                    */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void Vector<ValueType>::cwiseProduct( const _Vector& other )
{
    bool noSwapArgs = false;
    setVector( other, common::BinaryOp::MULT, noSwapArgs );
}

template<typename ValueType>
void Vector<ValueType>::cwiseDivision( const _Vector& other )
{
    bool noSwapArgs = false;
    setVector( other, common::BinaryOp::DIVIDE, noSwapArgs );
}

template<typename ValueType>
void Vector<ValueType>::scale( ValueType value )
{
    bool noSwapArgs = false;
    setScalar( value, common::BinaryOp::MULT, noSwapArgs );
}

/* ========================================================================= */

template<typename ValueType>
Scalar Vector<ValueType>::_l1Norm() const
{
    return Scalar( l1Norm() );
}

template<typename ValueType>
Scalar Vector<ValueType>::_l2Norm() const
{
    return Scalar( l2Norm() );
}

template<typename ValueType>
Scalar Vector<ValueType>::_maxNorm() const
{
    return Scalar( maxNorm() );
}

template<typename ValueType>
Scalar Vector<ValueType>::_maxDiffNorm( const _Vector& other ) const
{
    return Scalar( maxDiffNorm( other ) );
}

template<typename ValueType>
Scalar Vector<ValueType>::_min() const
{
    return Scalar( min() );
}

template<typename ValueType>
Scalar Vector<ValueType>::_max() const
{
    return Scalar( min() );
}

template<typename ValueType>
Scalar Vector<ValueType>::_sum() const
{
    return Scalar( sum() );
}

template<typename ValueType>
Scalar Vector<ValueType>::_dotProduct( const _Vector& other ) const
{
    return Scalar( dotProduct( other ) );
}

/* ========================================================================= */
/*        operator= < vector expression>                                     */
/* ========================================================================= */

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SV_SV<ValueType>& expression )
{
    SCAI_LOG_DEBUG( logger, "this = a * vector1 + b * vector2, check vector1.size() == vector2.size()" )

    const ValueType& alpha     = expression.getArg1().getArg1();
    const ValueType& beta      = expression.getArg2().getArg1();
    const Vector<ValueType>& x = expression.getArg1().getArg2();
    const Vector<ValueType>& y = expression.getArg2().getArg2();

    // Note: all checks are done the vector specific implementations

    vectorPlusVector( alpha, x, beta, y );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SV_S<ValueType>& expression )
{
    const Expression_SV<ValueType>& exp = expression.getArg1();
    const ValueType& alpha = exp.getArg1();
    const Vector<ValueType>& x = exp.getArg2();
    const ValueType& beta = expression.getArg2();

    vectorPlusScalar( alpha, x, beta );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SMV<ValueType>& expression )
{
    SCAI_LOG_INFO( logger, "this = alpha * matrix * vectorX -> this = alpha * matrix * vectorX + 0.0 * this" )

    const ValueType beta = 0;
    Expression_SV<ValueType> exp2( beta, *this );
    Expression_SMV_SV<ValueType> tmpExp( expression, exp2 );
    const Vector<ValueType>& vectorX = expression.getArg2().getArg2();

    if ( &vectorX != this )
    {
        // so this is not aliased to the vector on the rhs
        // as this will be used on rhs we do allocate it here
        // distribution is given by the row distribution of the matrix
        const Matrix<ValueType>& matrix = expression.getArg2().getArg1();
        dmemo::DistributionPtr dist = matrix.getRowDistributionPtr();
        allocate( dist );
        // values remain uninitialized as we assume that 0.0 * this (undefined) will
        // never be executed as an operation
    }

    return operator=( tmpExp );
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SVM<ValueType>& expression )
{   
    SCAI_LOG_INFO( logger, "this = alpha * vectorX * matrix -> this = alpha * vectorX * matrix + 0.0 * this" )
    const ValueType beta = 0;
    Expression_SV<ValueType> exp2( beta, *this );
    Expression_SVM_SV<ValueType> tmpExp( expression, exp2 );
    const Vector<ValueType>& vectorX = expression.getArg2().getArg1();
    
    if ( &vectorX != this )
    {   
        // so this is not aliased to the vector on the rhs
        // as this will be used on rhs we do allocate it here
        // distribution is given by the row distribution of the matrix
        const Matrix<ValueType>& matrix = expression.getArg2().getArg2();
        dmemo::DistributionPtr dist = matrix.getColDistributionPtr();
        allocate( dist );
        // values remain uninitialized as we assume that 0.0 * this (undefined) will
        // never be executed as an operation
    }
    
    return operator=( tmpExp );
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SMV_SV<ValueType>& expression )
{
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SMV_SV )" )
    const Expression_SMV<ValueType>& exp1 = expression.getArg1();
    const Expression_SV<ValueType>& exp2 = expression.getArg2();
    const ValueType& alpha = exp1.getArg1();
    const Expression_MV<ValueType> matrixTimesVectorExp = exp1.getArg2();
    const ValueType& beta = exp2.getArg1();
    const Vector<ValueType>& vectorY = exp2.getArg2();
    const Matrix<ValueType>& matrix = matrixTimesVectorExp.getArg1();
    const Vector<ValueType>& vectorX = matrixTimesVectorExp.getArg2();

    _Vector* resultPtr = this;
    _VectorPtr tmpResult;

    if ( &vectorX == this )
    {
        SCAI_LOG_DEBUG( logger, "Temporary for X required" )
        tmpResult.reset( _Vector::create( this->getCreateValue() ) );
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

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SVM_SV<ValueType>& expression )

{   
    SCAI_LOG_INFO( logger, "Vector::operator=( Expression_SVM_SV )" )
    const Expression_SVM<ValueType>& exp1 = expression.getArg1();
    const Expression_SV<ValueType>& exp2 = expression.getArg2();
    const Expression_VM<ValueType>& vectorTimesMatrixExp = exp1.getArg2();
    
    // resolve : result = alhpa * A * x + beta * y
    
    const ValueType& alpha = exp1.getArg1();
    const ValueType& beta = exp2.getArg1();
    const Vector<ValueType>& vectorY = exp2.getArg2();
    const Vector<ValueType>& vectorX = vectorTimesMatrixExp.getArg1();
    const Matrix<ValueType>& matrix = vectorTimesMatrixExp.getArg2();
    
    matrix.vectorTimesMatrix( *this, alpha, vectorX, beta, vectorY );
    
    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SV<ValueType>& expression )
{
    const ValueType& alpha = expression.getArg1();
    const Vector<ValueType>& x = expression.getArg2();

    vectorPlusVector( alpha, x, 0, x );

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
Vector<ValueType>& Vector<ValueType>::operator=( const Expression_SVV<ValueType>& expression )
{
    const ValueType& alpha = expression.getArg1();

    const Expression_VV<ValueType>& exp = expression.getArg2();
    const Vector<ValueType>& x = exp.getArg1();
    const Vector<ValueType>& y = exp.getArg2();

    vectorTimesVector( alpha, x, y );

    return *this;
}


template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator+=( const Expression_SV<ValueType>& exp )
{
    ValueType alpha = 1;
    ValueType beta  = exp.getArg1();

    const Vector<ValueType>& x = *this;
    const Vector<ValueType>& y = exp.getArg2();

    vectorPlusVector( alpha, x, beta, y );

    return *this;
}

template<typename ValueType>
Vector<ValueType>& Vector<ValueType>::operator-=( const Expression_SV<ValueType>& exp )
{
    ValueType alpha = 1;
    ValueType beta  = exp.getArg1();

    const Vector<ValueType>& x = *this;
    const Vector<ValueType>& y = exp.getArg2();

    vectorPlusVector( alpha, x, -beta, y );

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


/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Vector, SCAI_ARRAY_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
