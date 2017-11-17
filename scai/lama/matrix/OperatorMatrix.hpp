/**
 * @file OperatorMatrix.hpp
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
 * @brief Typed Matrix class to be used for matrix-free operations
 * @author Thomas Brandes
 * @date 28.06.2017
 */

#pragma once 

#include <scai/lama.hpp>

#include <scai/lama/matrix/Matrix.hpp>

namespace scai
{

namespace lama
{

/** 
 *  OperatorMatrix is an operator matrix class where the entries are not stored explicitly
 *  and mainly provides a a linear mapping equivalent to the matrix vector multiplication.
 *
 *  In contrary to the sparse and dense matrix they do not store explicitly the coordinates 
 *  (matrix-free linear mapping).
 */
template<typename ValueType>
class OperatorMatrix : public Matrix<ValueType>
{

public:

    /** Constructor must provide the distribution of input and output vector for the linear operator.
     *
     *  \param[in] targetDistribution distribution of target vector (row distribution of matrix)
     *  \param[in] sourceDistribution distribution of source vector (col distribution of matrix)
     */
    OperatorMatrix( dmemo::DistributionPtr targetDistribution, dmemo::DistributionPtr sourceDistribution ) : 

        Matrix<ValueType>( targetDistribution, sourceDistribution )

    {
    }

    /** Set distributions, context of the matrix ( base class member variables ) by other matrix */

    OperatorMatrix( const Matrix<ValueType>& other ) : Matrix<ValueType>( other )
    {
    }

    /** Linear mapping to be implemented: $result = alpha * this_matrix * x + beta * y$ */

    virtual void matrixTimesVectorImpl(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>& y ) const = 0;

    using Matrix<ValueType>::logger;

    virtual bool isConsistent() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }

    virtual const char* getTypeName() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }

    virtual void clear()
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }

    virtual void purge()
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }

    virtual void allocate( IndexType, IndexType )
    {
        COMMON_THROWEXCEPTION( "Operator matrix cannot be redefined" )
    }

    virtual void allocate(dmemo::DistributionPtr, dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setIdentity(dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setDenseData(dmemo::DistributionPtr, dmemo::DistributionPtr, const hmemo::_HArray&, Scalar)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setCSRData(dmemo::DistributionPtr, dmemo::DistributionPtr, 
                            IndexType, const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&, const hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setDIAData(dmemo::DistributionPtr, dmemo::DistributionPtr, IndexType, const hmemo::HArray<IndexType>&, const hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void assign(const _Matrix&)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void assignTranspose(const _Matrix&)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void assign(const _MatrixStorage&)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void assign(const _MatrixStorage&, dmemo::DistributionPtr, dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void buildLocalStorage(_MatrixStorage&) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual const _MatrixStorage& getLocalStorage() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void redistribute(dmemo::DistributionPtr, dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void redistribute( const dmemo::Redistributor&, dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setLocalRow( const hmemo::HArray<ValueType>&, IndexType, common::BinaryOp )
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setLocalColumn( const hmemo::HArray<ValueType>&, IndexType, common::BinaryOp )
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void getRow( _Vector&, IndexType ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void getRowLocal( _Vector&, IndexType ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void getColumn( _Vector&, IndexType ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setRow( const _Vector&, IndexType, common::BinaryOp )
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setColumn( const _Vector&, IndexType, common::BinaryOp )
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void getDiagonal( _Vector& ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setDiagonal( const _Vector& )
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setDiagonal(Scalar)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }

    virtual void reduce( _Vector&, 
                         const IndexType,
                         const common::BinaryOp,
                         const common::UnaryOp ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }

    virtual void scale(const _Vector&)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void scale(Scalar)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void conj()
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual ValueType getValue(IndexType, IndexType) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setValue(IndexType, IndexType, ValueType, common::BinaryOp)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual IndexType getNumValues() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }

    virtual void vectorTimesMatrixImpl( DenseVector<ValueType>&, const ValueType, 
                                        const DenseVector<ValueType>&, const ValueType, const DenseVector<ValueType>& ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void matrixTimesScalar(const _Matrix&, Scalar)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void matrixPlusMatrix(Scalar, const _Matrix&, Scalar, const _Matrix&)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void matrixTimesMatrix( _Matrix&, Scalar, const _Matrix&, Scalar, const _Matrix& ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual IndexType getLocalNumValues() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual IndexType getLocalNumRows() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual IndexType getLocalNumColumns() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void setContextPtr(hmemo::ContextPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual hmemo::ContextPtr getContextPtr() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual MatrixKind getMatrixKind() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void prefetch() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void wait() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual void invert(const _Matrix&)
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual ValueType l1Norm() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual ValueType l2Norm() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual ValueType maxNorm() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual ValueType maxDiffNorm( const _Matrix& ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual Matrix<ValueType>* newMatrix() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual Matrix<ValueType>* copy() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
    virtual Format getFormat() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
     virtual size_t getValueTypeSize() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
     virtual bool hasDiagonalProperty() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
     virtual void resetDiagonalProperty()
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
     virtual size_t getMemoryUsage() const
    {
        COMMON_THROWEXCEPTION( "not implemented for operator matrix" )
    }
};

}

}
