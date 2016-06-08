/**
 * @file CRTPMatrix.hpp
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

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>

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
class CRTPMatrix: public Matrix
{
public:

    /** Default constructor. */

    CRTPMatrix()
        : Matrix()
    {
    }

    /** Constructor for a replicated zero-matrix of a given shape. */

    CRTPMatrix( const IndexType numRows, const IndexType numColumns )
        : Matrix( numRows, numColumns )
    {
    }

    /** Constructor for a distributed matrix. */

    CRTPMatrix( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution )
        : Matrix( rowDistribution, colDistribution )
    {
    }

    void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        SCAI_REGION( "Mat.timesVector" )
        SCAI_LOG_INFO( logger, result << " = " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y )

        if ( &result == &y )
        {
            SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
        }
        else if ( &result == &x )
        {
            COMMON_THROWEXCEPTION( "alias: result = x is not handled, use temporary" )
        }
        else
        {
            // we inherit the row distribution of this matrix to result
            result.allocate( getRowDistributionPtr() );
            // no more to check: result.size() == mNumRows, getDistirubtion() == result.getDistribution()
        }

        SCAI_ASSERT_EQ_ERROR( x.getDistribution(), getColDistribution(), "mismatch distribution" )
        SCAI_ASSERT_EQ_ERROR( y.getDistribution(), getRowDistribution(), "mismatch distribution" )
        const DenseVector<ValueType>* denseX = dynamic_cast<const DenseVector<ValueType>*>( &x );
        const DenseVector<ValueType>* denseY = dynamic_cast<const DenseVector<ValueType>*>( &y );
        DenseVector<ValueType>* denseResult = dynamic_cast<DenseVector<ValueType>*>( &result );
        SCAI_ASSERT( denseX, x << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )
        // Note: in case of beta == 0, we might skip this test
        SCAI_ASSERT( denseY, y << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )
        SCAI_ASSERT( denseResult, result << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )
        static_cast<const Derived*>( this )->matrixTimesVectorImpl( *denseResult, alpha.getValue<ValueType>(), *denseX,
                beta.getValue<ValueType>(), *denseY );
    }

    void vectorTimesMatrix(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        SCAI_REGION( "Mat.vectorTimes" )
        SCAI_LOG_INFO( logger, result << " = " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y )

        if ( &result == &y )
        {
            SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
        }
        else if ( &result == &x )
        {
            COMMON_THROWEXCEPTION( "alias: result = x is not handled, use temporary" )
        }
        else
        {
            // we inherit the row distribution of this matrix to result
            result.allocate( getRowDistributionPtr() );
            // no more to check: result.size() == mNumRows, getDistirubtion() == result.getDistribution()
        }

        SCAI_ASSERT_EQ_ERROR( x.getDistribution(), getRowDistribution(), "" )
        SCAI_ASSERT_EQ_ERROR( y.getDistribution(), getColDistribution(), "" )
        const DenseVector<ValueType>* denseX = dynamic_cast<const DenseVector<ValueType>*>( &x );
        const DenseVector<ValueType>* denseY = dynamic_cast<const DenseVector<ValueType>*>( &y );
        DenseVector<ValueType>* denseResult = dynamic_cast<DenseVector<ValueType>*>( &result );
        SCAI_ASSERT( denseX, x << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )
        // Note: in case of beta == 0, we might skip this test
        SCAI_ASSERT( denseY, y << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )
        SCAI_ASSERT( denseResult, result << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )
        static_cast<const Derived*>( this )->vectorTimesMatrixImpl( *denseResult, alpha.getValue<ValueType>(), *denseX,
                beta.getValue<ValueType>(), *denseY );
    }

    /** @brief Get the row of a matrix
     *
     *  @param[out] row will contain the values of the queried row of this matrix
     *  @param[in]  globalRowIndex is the (global) index of the row to access
     */

    void getRow( Vector& row, const IndexType globalRowIndex ) const
    {
        using namespace scai::hmemo;
        SCAI_ASSERT_LT_ERROR( globalRowIndex, getNumRows(), "illegal index" )
        // row should be a DenseVector of same type, otherwise use a temporary
        common::shared_ptr<DenseVector<ValueType> > tmpVector;  // only allocated if needed
        DenseVector<ValueType>* typedRow = dynamic_cast<DenseVector<ValueType>*>( &row );

        if ( !typedRow )
        {
            // so we create a temporaray DenseVector of same type, has already correct size
            tmpVector.reset( new DenseVector<ValueType>() );
            typedRow = tmpVector.get();
            SCAI_LOG_INFO( logger, "temporary vector: " << *tmpVector << ", for row = " << row )
        }

        // if row is not same size and replicated, allocate / replicate it

        if ( !typedRow->getDistribution().isReplicated() || typedRow->size() != getNumColumns() )
        {
            dmemo::DistributionPtr dist( new dmemo::NoDistribution( getNumColumns() ) );
            typedRow->allocate( dist );
            SCAI_LOG_INFO( logger, "allocated vector for row, is now : " << *typedRow )
        }

        // on a replicated matrix each processor can fill the row by its own

        if ( getRowDistribution().isReplicated() )
        {
            SCAI_LOG_INFO( logger, "get local row " << globalRowIndex )
            static_cast<const Derived*>( this )->getLocalRow( *typedRow, globalRowIndex );
        }
        else
        {
            // on a distributed matrix, owner fills row and broadcasts it
            const dmemo::Communicator& comm = getRowDistribution().getCommunicator();
            // owner fills the row
            IndexType localRowIndex = getRowDistribution().global2local( globalRowIndex );
            IndexType owner = 0;
            // context where the row will have its valid data
            ContextPtr contextPtr = Context::getHostPtr();

            if ( localRowIndex != nIndex )
            {
                static_cast<const Derived*>( this )->getLocalRow( *typedRow, localRowIndex );
                owner = comm.getRank() + 1;
                SCAI_LOG_INFO( logger,
                               comm << ": owner of row " << globalRowIndex << ", local index = " << localRowIndex )
            }
            else
            {
                SCAI_LOG_INFO( logger, comm << ": not owner" )
            }

            owner = comm.sum( owner );
            SCAI_ASSERT_GT_ERROR( owner, 0, "Could not find owner of row " << globalRowIndex )
            owner -= 1;  // back to range 0, ..., size-1
            {
                // only owner has to keep a valid copy, others have write-only
                bool keep = owner == comm.getRank();
                WriteAccess<ValueType> rowAccess( typedRow->getLocalValues(), contextPtr, keep );
                SCAI_ASSERT_EQ_DEBUG( rowAccess.size(), getNumColumns(), "mismatch" );
                comm.bcast( rowAccess.get(), getNumColumns(), owner ); // bcast the row
                SCAI_LOG_INFO( logger, comm << ": bcast done, owner = " << owner )
            }
        }

        if ( tmpVector.get() )
        {
            SCAI_LOG_INFO( logger, "copy from tmp vector to result vector" )
            // if we have used a temporary vector, then we copy it back
            row = *tmpVector;   // implicit conversion
            SCAI_LOG_INFO( logger, "copy from tmp vector to result vector done" )
        }
    }

    using Matrix::setIdentity;
    using Matrix::operator=;

protected:

#ifndef SCAI_LOG_LEVEL_OFF

    // here in CRTP we do not use an own logger, just take it from the corresponding matrix class

    using Matrix::logger;

#endif

};

} /* end namespace lama */

} /* end namespace scai */
