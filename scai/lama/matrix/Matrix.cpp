/**
 * @file Matrix.cpp
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
 * @brief Implementation of methods for the abstract class Matrix<ValueType>
 * @author Thomas Brandes
 * @date 31.10.2017
 */

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/mepr/TypeList.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using common::TypeTraits;

namespace lama
{

/* ------------------------------------------------------------------------- */
/*    static methods                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
Matrix<ValueType>* Matrix<ValueType>::getMatrix( Format format )
{
    return reinterpret_cast<Matrix<ValueType>*>( _Matrix::getMatrix( format, TypeTraits<ValueType>::stype ) );
}

/* ------------------------------------------------------------------------- */
/*    Constructors / Destructor                                              */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
Matrix<ValueType>::Matrix() : _Matrix()
{
}

template<typename ValueType>
Matrix<ValueType>::Matrix( const IndexType numRows, const IndexType numColumns ) : 

    _Matrix( numRows, numColumns )

{
    SCAI_LOG_DEBUG( logger, "Matrix<" << TypeTraits<ValueType>::id() << "> ( "
                            << _Matrix::getNumRows() << " x " << _Matrix::getNumColumns() << " )" )
}

template<typename ValueType>
Matrix<ValueType>::Matrix( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution ) :

    _Matrix( rowDistribution, colDistribution )
{
}

template<typename ValueType>
Matrix<ValueType>::Matrix( const _Matrix& other, 
                           dmemo::DistributionPtr rowDistribution, 
                           dmemo::DistributionPtr colDistribution ) :

    _Matrix( other, rowDistribution, colDistribution )

{
}

template<typename ValueType>
Matrix<ValueType>::Matrix( const _Matrix& other ) :

    _Matrix( other )
{
}

template<typename ValueType>
Matrix<ValueType>::Matrix( const Matrix<ValueType>& other ) :

    _Matrix( other )
{
}

template<typename ValueType>
Matrix<ValueType>::~Matrix()
{
    SCAI_LOG_DEBUG( logger, "~Matrix<" << TypeTraits<ValueType>::id() << ">" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::scalar::ScalarType Matrix<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
size_t Matrix<ValueType>::getValueTypeSize() const
{
    return sizeof( ValueType );
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::matrixTimesVector(
    _Vector& result,
    const Scalar alpha,
    const _Vector& x,
    const Scalar beta,
    const _Vector& y ) const
{
    SCAI_REGION( "Mat.timesVector" )

    SCAI_LOG_INFO( logger, 
                   "result = " << alpha << " * M<" << this->getValueType() << ">[" << this->getNumRows() << " x " << this->getNumColumns() << "]"
                   << " * x [ " << x.size() << "] + " << beta << " * y[ " << y.size() << "]" )

    if ( &result == &y )
    {
        SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
    }
    else
    {
        // we inherit the row distribution of this matrix to result

        result.allocate( this->getRowDistributionPtr() );
    }

    if ( x.getVectorKind() != VectorKind::DENSE || x.getValueType() != this->getValueType() || &result == &x || x.getDistribution() != this->getColDistribution() )
    {
        SCAI_UNSUPPORTED( "alpha * M * x, x requires temporary DenseVector<" << this->getValueType() << ">" )

        DenseVector<ValueType> tmpX( x, this->getColDistributionPtr() );
        matrixTimesVector( result, alpha, tmpX, beta, y );
        return;
    }

    const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );

    // Note: in case of beta == 0, we might skip this test

    if ( y.getVectorKind() != VectorKind::DENSE || y.getValueType() != this->getValueType() || y.getDistribution() != this->getRowDistribution() )
    {
        SCAI_UNSUPPORTED( "temporary DenseVector<" << this->getValueType() << "> required for y in alpha * M * x + beta * y" )
        DenseVector<ValueType> tmpY( y, this->getRowDistributionPtr() );
        matrixTimesVector( result, alpha, x, beta, tmpY );
        return;
    }

    const DenseVector<ValueType>& denseY = reinterpret_cast<const DenseVector<ValueType>&>( y );

    if ( result.getVectorKind() != VectorKind::DENSE || result.getValueType() != this->getValueType() )
    {
        SCAI_UNSUPPORTED( "temporary DenseVector<" << this->getValueType() << "> required for result in alpha * M * x + beta * y" )
        DenseVector<ValueType> tmpResult( this->getRowDistributionPtr() );
        matrixTimesVector( tmpResult, alpha, x, beta, y );
        result = tmpResult;
        return;
    }

    DenseVector<ValueType>& denseResult = reinterpret_cast<DenseVector<ValueType>&>( result );

    const ValueType alphaV = alpha.getValue<ValueType>();
    const ValueType betaV  = beta.getValue<ValueType>();

    // Now call the typed version implemented by derived class

    matrixTimesVectorImpl( denseResult, alphaV, denseX, betaV, denseY );
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::vectorTimesMatrix(
    _Vector& result,
    const Scalar alpha,
    const _Vector& x,
    const Scalar beta,
    const _Vector& y ) const
{
    SCAI_REGION( "Mat.vectorTimes" )

    SCAI_LOG_INFO( logger, result << " = " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y )

    if ( x.getVectorKind() != VectorKind::DENSE || x.getValueType() != this->getValueType() || &result == &x || x.getDistribution() != this->getRowDistribution() )
    {
        SCAI_UNSUPPORTED( "temporary DenseVector<" << this->getValueType() << "> required for x in alpha * M * x + beta * y" )
        DenseVector<ValueType> tmpX( x, this->getRowDistributionPtr() );
        vectorTimesMatrix( result, alpha, tmpX, beta, y );
        return;
    }

    const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );

    if ( y.getVectorKind() != VectorKind::DENSE || y.getValueType() != this->getValueType() || y.getDistribution() != this->getColDistribution() )
    {
        SCAI_UNSUPPORTED( "temporary DenseVector<" << this->getValueType() << "> required for y in alpha * x * M + beta * y" )
        DenseVector<ValueType> tmpY( y, this->getColDistributionPtr() );
        vectorTimesMatrix( result, alpha, x, beta, tmpY );
        return;
    }

    const DenseVector<ValueType>& denseY = reinterpret_cast<const DenseVector<ValueType>&>( y );

    if ( result.getVectorKind() != VectorKind::DENSE || result.getValueType() != this->getValueType() )
    {
        SCAI_UNSUPPORTED( "temporary DenseVector<" << this->getValueType() << "> required for result in alpha * M * x + beta * y" )
        DenseVector<ValueType> tmpResult( this->getColDistributionPtr() );
        vectorTimesMatrix( tmpResult, alpha, x, beta, y );
        result = tmpResult;
        return;
    }

    if ( &result == &y )
    {
        SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
    }
    else
    {
        result.allocate( this->getColDistributionPtr() );
    }

    DenseVector<ValueType>& denseResult = reinterpret_cast<DenseVector<ValueType>&>( result );

    const ValueType alphaV = alpha.getValue<ValueType>();
    const ValueType betaV  = beta.getValue<ValueType>();

    if ( this->getColDistribution().getCommunicator().getSize() == 1 )
    {
        // Each processor has full columns, resultVector is replicated, communication only needed to sum up results
        // use routine provided by this CRTP

        this->vectorTimesMatrixRepCols( denseResult, alphaV, denseX, betaV, denseY );
    }
    else
    {
        this->vectorTimesMatrixImpl( denseResult, alphaV, denseX, betaV, denseY );
    }
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::setRow( 
    const _Vector& row, 
    const IndexType globalRowIndex,
    const common::binary::BinaryOp op )
{
    SCAI_ASSERT_EQ_ERROR( row.size(), this->getNumColumns(), "row size mismatch" )

    SCAI_LOG_DEBUG( logger, "setRow " << globalRowIndex << ": row = " << row << ", op = " << op )

    bool needsTmp = false;

    if ( row.getValueType() != this->getValueType() )
    {
        needsTmp = true;
        SCAI_UNSUPPORTED( "setRow, matrix has type " << this->getValueType() 
                           << ", row has type " << row.getValueType() << ", use temporary" )
    }
    if ( ! row.getDistribution().isReplicated() )
    {
        needsTmp = true;
        SCAI_UNSUPPORTED( "setRow, row is not replicated, use temporary" )
    }
    if ( row.getVectorKind() != VectorKind::DENSE )
    {
        needsTmp = true;
        SCAI_UNSUPPORTED( "setRow, row is not DENSE vector" )
    }

    if ( needsTmp )
    {
        DenseVector<ValueType> tmpRow( row );
        tmpRow.replicate();
        setRow( tmpRow, globalRowIndex, op );
        return;
    }

    using namespace scai::hmemo;

    SCAI_ASSERT_VALID_INDEX_ERROR( globalRowIndex, this->getNumRows(), "illegal row index" )

    // row should be a DenseVector of same type, otherwise use a temporary

    common::shared_ptr<DenseVector<ValueType> > tmpVector;  // only allocated if needed

    const DenseVector<ValueType>* typedRow = dynamic_cast<const DenseVector<ValueType>*>( &row );

    SCAI_ASSERT_ERROR( typedRow, "illegal dynamic cast" )

    SCAI_ASSERT_ERROR( typedRow->getDistribution().isReplicated(), "cannot set distributed row" )

    SCAI_ASSERT_EQ_ERROR( typedRow->size(), this->getNumColumns(), "row to set has wrong size" )

    // owner sets the row, maybe each processor for replicated row distribution

    IndexType localRowIndex = this->getRowDistribution().global2local( globalRowIndex );

    if ( localRowIndex != nIndex )
    {
        this->setLocalRow( typedRow->getLocalValues(), localRowIndex, op );
    }
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::setColumn( 
    const _Vector& column,
    const IndexType colIndex,
    const common::binary::BinaryOp op )
{
    using namespace scai::hmemo;

    SCAI_ASSERT_VALID_INDEX_ERROR( colIndex, this->getNumColumns(), "illegal col index" )

    // col should be a DenseVector of same type, otherwise use a temporary

    common::shared_ptr<const DenseVector<ValueType> > tmpVector;  // only allocated if needed

    const DenseVector<ValueType>* typedColumn = dynamic_cast<const DenseVector<ValueType>*>( &column );

    if ( !typedColumn )
    {
        // so we create a temporaray DenseVector of same type, has already correct size
        tmpVector.reset( new DenseVector<ValueType>( column ) );
        typedColumn = tmpVector.get();
    }

    SCAI_ASSERT_EQ_ERROR( typedColumn->getDistribution(), this->getRowDistribution(), "distribution mismatch" )

    this->setLocalColumn( typedColumn->getLocalValues(), colIndex, op );
}

/* ========================================================================= */

template<typename ValueType>
void Matrix<ValueType>::vectorTimesMatrixRepCols(
    DenseVector<ValueType>& denseResult,
    const ValueType alphaValue,
    const DenseVector<ValueType>& denseX,
    const ValueType betaValue,
    const DenseVector<ValueType>& denseY ) const
{
    SCAI_REGION( "Mat.vectorTimesMatrixRepCols" )

    const hmemo::HArray<ValueType>& localY = denseY.getLocalValues();
    const hmemo::HArray<ValueType>& localX = denseX.getLocalValues();

    hmemo::HArray<ValueType>& localResult = denseResult.getLocalValues();

    const dmemo::Distribution& colDist = this->getColDistribution();

    // this routine is only for non-replicated columns, i.e. mHaloData is empty

    SCAI_ASSERT( 1, colDist.getNumPartitions() );

    const dmemo::Distribution& rowDist = this->getRowDistribution();
    const dmemo::Communicator& comm = rowDist.getCommunicator();

    const MatrixStorage<ValueType>& localData = reinterpret_cast<const MatrixStorage<ValueType>&>( this->getLocalStorage() );

    if ( comm.getRank() == 0 )
    {
        // only one single processor adds beta * y
        localData.vectorTimesMatrix( localResult, alphaValue, localX, betaValue, localY );
    }
    else
    {
        localData.vectorTimesMatrix( localResult, alphaValue, localX, ValueType( 0 ), localY );
    }

    if ( comm.getSize() >  1 )
    {
        // Sum up all incarnations of localResult

        comm.sumArray( localResult );
    }
}

/* ========================================================================= */

template<typename ValueType>
Scalar Matrix<ValueType>::_l1Norm() const
{
    return Scalar( l1Norm() );
}

template<typename ValueType>
Scalar Matrix<ValueType>::_l2Norm() const
{
    return Scalar( l2Norm() );
}

template<typename ValueType>
Scalar Matrix<ValueType>::_maxNorm() const
{
    return Scalar( maxNorm() );
}

template<typename ValueType>
Scalar Matrix<ValueType>::_maxDiffNorm( const _Matrix& other ) const
{
    return Scalar( maxDiffNorm( other ) );
}

template<typename ValueType>
typename Matrix<ValueType>::RealType Matrix<ValueType>::maxDiffNorm( const _Matrix& other ) const
{
    typedef typename Matrix<ValueType>::RealType RealType;

    IndexType nRows = getNumRows();
    IndexType nCols = getNumColumns();

    SCAI_ASSERT_EQUAL( nRows, other.getNumRows(), "size mismatch" )
    SCAI_ASSERT_EQUAL( nCols, other.getNumColumns(), "size mismatch" )

    DenseVector<ValueType> row;
    DenseVector<ValueType> rowOther;

    RealType diff = 0;

    // now traverse  all rows

    for ( IndexType i = 0; i < nRows; ++i )
    {
        // Note: rows will be broadcast in case of distributed matrices

        getRow( row, i );
        other.getRow( rowOther, i );

        RealType diffRow = row.maxDiffNorm( rowOther );

        if ( diffRow > diff )
        {
            diff = diffRow;
        }
    }

    return diff;
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Matrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */

