/**
 * @file AbstractMatrix.hpp
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
 * @brief Matrix class to be used for matrix-free operations
 * @author Thomas Brandes
 * @date 28.06.2017
 */

#pragma once 

#include <scai/lama.hpp>

#include <scai/lama/matrix/_Matrix.hpp>

namespace scai
{

namespace lama
{

/** This abstract matrix class provides a matrix class where all pure methods are provided 
 *  but just throw an exception.
 *
 *  It might be used for matrix-free solvers where only the matrix * vector operations is
 *  required.
 */
class AbstractMatrix : public Matrix
{

public:

    /** Set distributions, context of the matrix ( base class member variables ) by other matrix */

    AbstractMatrix( const Matrix& other ) : Matrix( other )
    {
    }

    AbstractMatrix( dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist ) : Matrix( rowDist, colDist )
    {
    }

    virtual bool isConsistent() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual const char* getTypeName() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual void clear()
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual void purge()
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual void allocate(IndexType, IndexType)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual void allocate(dmemo::DistributionPtr, dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setIdentity(dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setDenseData(dmemo::DistributionPtr, dmemo::DistributionPtr, const hmemo::_HArray&, Scalar)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setCSRData(dmemo::DistributionPtr, dmemo::DistributionPtr, 
                            IndexType, const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&, const hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setDIAData(dmemo::DistributionPtr, dmemo::DistributionPtr, IndexType, const hmemo::HArray<IndexType>&, const hmemo::_HArray&)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void assign(const Matrix&)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void assignTranspose(const Matrix&)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void assign(const _MatrixStorage&)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void assign(const _MatrixStorage&, dmemo::DistributionPtr, dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void buildLocalStorage(_MatrixStorage&) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual const _MatrixStorage& getLocalStorage() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void redistribute(dmemo::DistributionPtr, dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void redistribute( const dmemo::Redistributor&, dmemo::DistributionPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void getRow( _Vector&, IndexType ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void getRowLocal( _Vector&, IndexType ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void getColumn( _Vector&, IndexType ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setRow( const _Vector&, IndexType, common::binary::BinaryOp )
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setColumn( const _Vector&, IndexType, common::binary::BinaryOp )
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void getDiagonal( _Vector& ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setDiagonal( const _Vector& )
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setDiagonal(Scalar)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual void reduce( _Vector&, 
                         const IndexType,
                         const common::binary::BinaryOp,
                         const common::unary::UnaryOp ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual void scale(const _Vector&)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void scale(Scalar)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void conj()
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual Scalar getValue(IndexType, IndexType) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setValue(IndexType, IndexType, Scalar, common::binary::BinaryOp)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual IndexType getNumValues() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual void matrixTimesVector( _Vector&, const Scalar, const _Vector&, const Scalar, const _Vector& ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }

    virtual void vectorTimesMatrix( _Vector&, Scalar, const _Vector&, Scalar, const _Vector& ) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void matrixTimesScalar(const Matrix&, Scalar)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void matrixPlusMatrix(Scalar, const Matrix&, Scalar, const Matrix&)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void matrixTimesMatrix(Matrix&, Scalar, const Matrix&, Scalar, const Matrix&) const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual IndexType getLocalNumValues() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual IndexType getLocalNumRows() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual IndexType getLocalNumColumns() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void setContextPtr(hmemo::ContextPtr)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual hmemo::ContextPtr getContextPtr() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual MatrixKind getMatrixKind() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void prefetch() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void wait() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual void invert(const Matrix&)
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual Scalar l1Norm() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual Scalar l2Norm() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual Scalar maxNorm() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual Matrix* newMatrix() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
    virtual Matrix* copy() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
     virtual common::scalar::ScalarType getValueType() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
     virtual Format::MatrixStorageFormat getFormat() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
     virtual size_t getValueTypeSize() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
     virtual bool hasDiagonalProperty() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
     virtual void resetDiagonalProperty()
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
     virtual size_t getMemoryUsage() const
    {
        COMMON_THROWEXCEPTION( "not implemented for abstract matrix" )
    }
};

}

}
