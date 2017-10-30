/**
 * @file CRTPMatrix.hpp
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
 * @brief Template class for common methods of SparseMatrix<ValueType> and DenseMatrix<ValueType>
 *        to deal with polymorphism
 * @author Thomas Brandes
 * @date 09.08.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/Matrix.hpp>

// local library
#include <scai/lama/expression/all.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/dmemo/NoDistribution.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/unsupported.hpp>

namespace scai
{

namespace lama
{

/** This template class supports static polymorphism to define
 *  common routines for base classes SparseMatrix<ValueType> and DenseMatrix<ValueType>.
 *
 *  Therefore it uses the Curiously Recurring Template Pattern (CRTP)
 *  as a C++ idiom where a class X derived from this template is
 *  itself a template argument.
 *
 *  @tparam Derived is the derived class of Matrix
 *  @tparam ValueType specifies the type of the matrix values
 *
 *  @todo: create, copy should also be defined here
 */

template<class Derived, typename ValueType>
class CRTPMatrix
{
public:

    void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        const Derived& m = reinterpret_cast<const Derived&>( *this );

        SCAI_REGION( "Mat.timesVector" )

        SCAI_LOG_INFO( Derived::logger, 
                       "result = " << alpha << " * M<" << m.getValueType() << ">[" << m.getNumRows() << " x " << m.getNumColumns() << "]"
                       << " * x [ " << x.size() << "] + " << beta << " * y[ " << y.size() << "]" )

        if ( &result == &y )
        {
            SCAI_LOG_DEBUG( Derived::logger, "alias: result = y is well handled" )
        }
        else
        {
            // we inherit the row distribution of this matrix to result

            result.allocate( m.getRowDistributionPtr() );
        }

        if ( x.getVectorKind() != Vector::DENSE || x.getValueType() != m.getValueType() || &result == &x || x.getDistribution() != m.getColDistribution() )
        {
            SCAI_UNSUPPORTED( "alpha * M * x, x requires temporary DenseVector<" << m.getValueType() << ">" )

            DenseVector<ValueType> tmpX( x, m.getColDistributionPtr() );
            matrixTimesVector( result, alpha, tmpX, beta, y );
            return;
        }

        const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );

        // Note: in case of beta == 0, we might skip this test

        if ( y.getVectorKind() != Vector::DENSE || y.getValueType() != m.getValueType() || y.getDistribution() != m.getRowDistribution() )
        {
            SCAI_UNSUPPORTED( "temporary DenseVector<" << m.getValueType() << "> required for y in alpha * M * x + beta * y" )
            DenseVector<ValueType> tmpY( y, m.getRowDistributionPtr() );
            matrixTimesVector( result, alpha, x, beta, tmpY );
            return;
        }

        const DenseVector<ValueType>& denseY = reinterpret_cast<const DenseVector<ValueType>&>( y );

        if ( result.getVectorKind() != Vector::DENSE || result.getValueType() != m.getValueType() )
        {
            SCAI_UNSUPPORTED( "temporary DenseVector<" << m.getValueType() << "> required for result in alpha * M * x + beta * y" )
            DenseVector<ValueType> tmpResult( m.getRowDistributionPtr() );
            matrixTimesVector( tmpResult, alpha, x, beta, y );
            result = tmpResult;
            return;
        }

        DenseVector<ValueType>& denseResult = reinterpret_cast<DenseVector<ValueType>&>( result );

        const ValueType alphaV = alpha.getValue<ValueType>();
        const ValueType betaV  = beta.getValue<ValueType>();

        m.matrixTimesVectorImpl( denseResult, alphaV, denseX, betaV, denseY );
    }

    void vectorTimesMatrix(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        const Derived& m = reinterpret_cast<const Derived&>( *this );

        SCAI_REGION( "Mat.vectorTimes" )

        SCAI_LOG_INFO( Derived::logger, result << " = " << alpha << " * " << m << " * " << x << " + " << beta << " * " << y )

        if ( x.getVectorKind() != Vector::DENSE || x.getValueType() != m.getValueType() || &result == &x || x.getDistribution() != m.getRowDistribution() )
        {
            SCAI_UNSUPPORTED( "temporary DenseVector<" << m.getValueType() << "> required for x in alpha * M * x + beta * y" )
            DenseVector<ValueType> tmpX( x, m.getRowDistributionPtr() );
            vectorTimesMatrix( result, alpha, tmpX, beta, y );
            return;
        }

        const DenseVector<ValueType>& denseX = reinterpret_cast<const DenseVector<ValueType>&>( x );

        if ( y.getVectorKind() != Vector::DENSE || y.getValueType() != m.getValueType() || y.getDistribution() != m.getColDistribution() )
        {
            SCAI_UNSUPPORTED( "temporary DenseVector<" << m.getValueType() << "> required for y in alpha * x * M + beta * y" )
            DenseVector<ValueType> tmpY( y, m.getColDistributionPtr() );
            vectorTimesMatrix( result, alpha, x, beta, tmpY );
            return;
        }

        const DenseVector<ValueType>& denseY = reinterpret_cast<const DenseVector<ValueType>&>( y );

        if ( result.getVectorKind() != Vector::DENSE || result.getValueType() != m.getValueType() )
        {
            SCAI_UNSUPPORTED( "temporary DenseVector<" << m.getValueType() << "> required for result in alpha * M * x + beta * y" )
            DenseVector<ValueType> tmpResult( m.getColDistributionPtr() );
            vectorTimesMatrix( tmpResult, alpha, x, beta, y );
            result = tmpResult;
            return;
        }

        if ( &result == &y )
        {
            SCAI_LOG_DEBUG( Derived::logger, "alias: result = y is well handled" )
        }
        else
        {
            result.allocate( m.getColDistributionPtr() );
        }

        DenseVector<ValueType>& denseResult = reinterpret_cast<DenseVector<ValueType>&>( result );

        const ValueType alphaV = alpha.getValue<ValueType>();
        const ValueType betaV  = beta.getValue<ValueType>();

        if ( m.getColDistribution().getCommunicator().getSize() == 1 )
        {
            // Each processor has full columns, resultVector is replicated, communication only needed to sum up results
            // use routine provided by this CRTP

            m.vectorTimesMatrixRepCols( denseResult, alphaV, denseX, betaV, denseY );
        }
        else
        {
            m.vectorTimesMatrixImpl( denseResult, alphaV, denseX, betaV, denseY );
        }
    }

    /** Implementation of Matrix::setRow for all typed matrices
     *
     *  Note: all derived classes must provide setLocalRow( rowArray, localRowIndex, op )
     */
    void setRow( const Vector& row, const IndexType globalRowIndex,
                 const common::binary::BinaryOp op )
    {
        Derived& m = reinterpret_cast<Derived&>( *this );

        SCAI_ASSERT_EQ_ERROR( row.size(), m.getNumColumns(), "row size mismatch" )

        bool needsTmp = false;

        if ( row.getValueType() != m.getValueType() )
        {
            needsTmp = true;
            SCAI_UNSUPPORTED( "setRow, matrix has type " << m.getValueType() 
                               << ", row has type " << row.getValueType() << ", use temporary" )
        }
        if ( ! row.getDistribution().isReplicated() )
        {
            needsTmp = true;
            SCAI_UNSUPPORTED( "setRow, row is not replicated, use temporary" )
        }
        if ( row.getVectorKind() != Vector::DENSE )
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

        SCAI_ASSERT_VALID_INDEX_ERROR( globalRowIndex, m.getNumRows(), "illegal row index" )

        // row should be a DenseVector of same type, otherwise use a temporary

        std::shared_ptr<DenseVector<ValueType> > tmpVector;  // only allocated if needed

        const DenseVector<ValueType>* typedRow = dynamic_cast<const DenseVector<ValueType>*>( &row );

        SCAI_ASSERT_ERROR( typedRow, "illegal dynamic cast" )

        SCAI_ASSERT_ERROR( typedRow->getDistribution().isReplicated(), "cannot set distributed row" )

        SCAI_ASSERT_EQ_ERROR( typedRow->size(), m.getNumColumns(), "row to set has wrong size" )

        // owner sets the row, maybe each processor for replicated row distribution

        IndexType localRowIndex = m.getRowDistribution().global2local( globalRowIndex );

        if ( localRowIndex != nIndex )
        {
            m.setLocalRow( typedRow->getLocalValues(), localRowIndex, op );
        }
    }

    /** Implementation of Matrix::setColumn for all typed matrices
     *
     *  The method is implemented by setting the local part of the column on each partition.
     *  All derived classes must provide setLocalColum( colArray, colIndex, op )
     */
    void setColumn( const Vector& column,
                    const IndexType colIndex,
                    const common::binary::BinaryOp op )
    {
        Derived& m = reinterpret_cast<Derived&>( *this );

        using namespace scai::hmemo;

        SCAI_ASSERT_VALID_INDEX_ERROR( colIndex, m.getNumColumns(), "illegal col index" )

        // col should be a DenseVector of same type, otherwise use a temporary

        std::shared_ptr<const DenseVector<ValueType> > tmpVector;  // only allocated if needed

        const DenseVector<ValueType>* typedColumn = dynamic_cast<const DenseVector<ValueType>*>( &column );

        if ( !typedColumn )
        {
            // so we create a temporaray DenseVector of same type, has already correct size
            tmpVector.reset( new DenseVector<ValueType>( column ) );
            typedColumn = tmpVector.get();
        }

        SCAI_ASSERT_EQ_ERROR( typedColumn->getDistribution(), m.getRowDistribution(), "distribution mismatch" )

        m.setLocalColumn( typedColumn->getLocalValues(), colIndex, op );
    }

    /** This method is the same for dense/sparse matrices as column distribution is replicated */

    void vectorTimesMatrixRepCols(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>& denseY ) const
    {
        const Derived& m = reinterpret_cast<const Derived&>( *this );

        SCAI_REGION( "Mat.Sp.vectorTimesMatrixRepCols" )

        const hmemo::HArray<ValueType>& localY = denseY.getLocalValues();
        const hmemo::HArray<ValueType>& localX = denseX.getLocalValues();

        hmemo::HArray<ValueType>& localResult = denseResult.getLocalValues();

        const dmemo::Distribution& colDist = m.getColDistribution();

        // this routine is only for non-replicated columns, i.e. mHaloData is empty

        SCAI_ASSERT( 1, colDist.getNumPartitions() );

        const dmemo::Distribution& rowDist = m.getRowDistribution();
        const dmemo::Communicator& comm = rowDist.getCommunicator();

        const MatrixStorage<ValueType>& localData = static_cast<const Derived*>( this )->getLocalStorage();

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
};

} /* end namespace lama */

} /* end namespace scai */
