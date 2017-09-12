/**
 * @file SparseMatrix.hpp
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
 * @brief Definition of class for distributed sparse matrices.
 * @author Jiri Kraus, Thomas Brandes
 * @date 06.06.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/CRTPMatrix.hpp>

// local library
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/DenseVector.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/common/function.hpp>

namespace scai
{

namespace lama
{

// forward declarations

class Vector;
template<typename > class DenseMatrix;

/**
 * @brief A SparseMatrix represents a distributed 2D sparse matrix with elements of type ValueType.
 *
 * The rows of the matrix are distributed among the processors specified by distribution.
 *
 * The local rows of the matrix are split into a local and a halo part corresponding to the
 * column distribution colDistribution.
 *
 * It is possible to use different storage formats for the local and halo part,
 * but both representations must have the same value type.
 *
 * @tparam ValueType is the value type of the matrix values.
 */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT SparseMatrix: 

    public CRTPMatrix<SparseMatrix<ValueType>, ValueType>,
    public Matrix

{

public:

    typedef ValueType MatrixValueType; //!< This is the type of the matrix values.

    /** Getter for the type name of the class. */

    static const char* typeName();

    /* Implementation of pure method of class Matrix. */

    virtual const char* getTypeName() const
    {
        return typeName();
    }

    /** Getter routine for local part of the sparse matrix. */

    const MatrixStorage<ValueType>& getLocalStorage() const
    {
        return *mLocalData;
    }

    /** Getter routine for halo part of the sparse matrix. */

    const MatrixStorage<ValueType>& getHaloStorage() const
    {
        return *mHaloData;
    }

    /**
     * @brief Implemementation of pure routine
     */
    virtual Format::MatrixStorageFormat getFormat() const
    {
        return mLocalData->getFormat();
    }

    /** This method sets all relevant data of a sparse matrix and checks consistency */

    void set(
        common::shared_ptr<MatrixStorage<ValueType> > localData,
        common::shared_ptr<MatrixStorage<ValueType> > haloData,
        const dmemo::Halo& halo,
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist );

    /**
     * @brief Constructor of a replicated sparse matrix with global storage.
     *
     * @param[in] storage   is a matrix storage of any type
     *
     * The distribution for rows and colums will be replicated, the halo
     * remains empty.
     *
     * Note: this constructor will not create a copy of the storage data.
     *       but just join the ownership
     */
    SparseMatrix( common::shared_ptr<MatrixStorage<ValueType> > storage );

    SparseMatrix( common::shared_ptr<MatrixStorage<ValueType> > storage, dmemo::DistributionPtr rowDist );

    SparseMatrix(
        common::shared_ptr<MatrixStorage<ValueType> > storage,
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist );

    /** Constructor of a sparse matrix with local and halo data available. */

    SparseMatrix(
        common::shared_ptr<MatrixStorage<ValueType> > localData,
        common::shared_ptr<MatrixStorage<ValueType> > haloData,
        const dmemo::Halo& halo,
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist );

    /** Override also the default copy constructor that does not make a
     *  deep copy of the input matrix due to the use of shared pointers.
     */
    SparseMatrix( const SparseMatrix<ValueType>& other );

    /** Implementation of abstract method for sparse matrices. */

    virtual bool isConsistent() const;

    /* Implementation of pure method of class Matrix. */

    virtual void invert( const Matrix& other );

    /* Implementation of pure method of class Matrix. */

    virtual void setContextPtr( const hmemo::ContextPtr context );

    /* Implementation of pure method of class Matrix. */

    virtual hmemo::ContextPtr getContextPtr() const;

    /** Implementation for Matrix::setDenseData */

    virtual void setDenseData(
        dmemo::DistributionPtr rowDistribution,
        dmemo::DistributionPtr colDistribution,
        const hmemo::_HArray& values,
        Scalar eps = Scalar( 0 ) );

    /** Implementation for pure method Matrix::setCSRData. */

    virtual void setCSRData(
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values );

    /** Implementation for pure method Matrix::setDIAData. */

    virtual void setDIAData(
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist,
        const IndexType numDiagonals,
        const hmemo::HArray<IndexType>& offsets,
        const hmemo::_HArray& values );

    /* Implementation of pure method of class Matrix. */

    virtual void clear();

    /* Implementation of pure method of class Matrix. */

    virtual void purge();

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( const IndexType numRows, const IndexType numColumns );

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /* Before overriding the virtual function make the other routine setIdentity( int n ) visible */

    using Matrix::setIdentity;

    virtual void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const
    {
        CRTPMatrix<SparseMatrix<ValueType>, ValueType>::matrixTimesVector( result, alpha, x, beta, y );
    }

    virtual void vectorTimesMatrix(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const 
    {
        CRTPMatrix<SparseMatrix<ValueType>, ValueType>::vectorTimesMatrix( result, alpha, x, beta, y );
    }

    /** @brief Implementation of pure method Matrix::getColumn 
     *
     *  It is recommended to call getColumn with a SparseVector for a sparse matrix.
     */
    virtual void getColumn( Vector& column, const IndexType globalColIndex ) const;

    virtual void setRow( const Vector& row,
                         const IndexType globalRowIndex,
                         const common::binary::BinaryOp op )
    {
        CRTPMatrix<SparseMatrix<ValueType>, ValueType>::setRow( row, globalRowIndex, op );
    }

    virtual void setColumn(
        const Vector& column,
        const IndexType globalColIndex,
        const common::binary::BinaryOp op )
    {
        CRTPMatrix<SparseMatrix<ValueType>, ValueType>::setColumn( column, globalColIndex, op );
    }

    /** Set matrix to a identity square matrix with same row and column distribution. */

    virtual void setIdentity( dmemo::DistributionPtr distribution );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const Matrix& other );

    /** Method that assigns a sparse matrix, specialization of assign( const Matrix& ) */

    void assign( const SparseMatrix<ValueType>& matrix );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const _MatrixStorage& storage );

    /*
     * @brief Setting (distributed) matrix with any local matrix data
     *
     * @param[in] other is local (sparse) matrix data containing all values to be set
     *
     * Size of other matrix must be exactly the same as this matrix. This routine might imply type and storage format
     * changes.
     *
     void assignLocal( const _MatrixStorage&
     {
     COMMON_THROWEXCEPTION( "not available yet" )
     }
     */

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Implementation of of pure method of class Matrix. */

    virtual void buildLocalStorage( _MatrixStorage& storage ) const;

    /**
     * @brief This method will swap all member variables of the two sparse matrices.
     *
     *  This operation might be useful in iteration loops where a sparse matrix
     *  is updated each iteration. It is more convenient than a solution that
     *  is based on using pointers in the application.
     *
     * @param[in] other  other sparse matrix
     */
    void swap( SparseMatrix<ValueType>& other );

    /**
     * @brief Destructor. Releases all allocated resources.
     */
    virtual ~SparseMatrix();

    /* Implementation of pure method of class Matrix. */

    virtual void getDiagonal( Vector& diagonal ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Vector& diagonal );

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Scalar scalar );

    /* Implementation of pure method of class Matrix. */

    virtual void reduce( 
        Vector& v, 
        const IndexType dim, 
        const common::binary::BinaryOp reduceOp, 
        const common::unary::UnaryOp elemOp ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Vector& scaling );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Scalar scaling );

    /* Implementation of pure method of class Matrix. */

    virtual void conj();

    /*
     *  Set local data of the matrix.
     *  The local part of the distributed matrix will be splitted into local / halo part.
     *  corresponding to the column distribution, builds new halo
     */

    /* Implemenation of pure method of class Matrix */

    virtual void matrixTimesScalar( const Matrix& other, const Scalar alpha );

    /**
     * @brief Same as matrixTimesVectorImpl but with multiple results and y
     *
     * Note: result and y must have same shape
     */
    virtual void matrixTimesVectorNImpl(
        DenseMatrix<ValueType>& result,
        const ValueType alpha,
        const DenseMatrix<ValueType>& x,
        const ValueType beta,
        const DenseMatrix<ValueType>& y ) const;

    /* Implementation of method needed for CRTPMatrix */

    void vectorTimesMatrixImpl(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>& y ) const;

    /**
     * @brief Same as matrixTimesVector but with vectors result, x, and y of same value type.
     */
    void matrixTimesVectorImpl(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>& y ) const;

    /**
     * @brief Operation on distributed matrix with halo exchange, sync version
     *
     * @param[out]    localResult is the result array
     * @param[in]     localX is the array with local values
     * @param[in,out] haloX is a temporary array keeping the non-local halo values of X
     * @param[in]     localF is the operation called with local matrix and localX
     * @param[in]     haloF is the operation called with halo matrix and haloX
     *
     * The synchronous version updates haloX by halo communication between the involved
     * processors and calls the functions localF and haloF to build the result array.
     *
     *  - localResult is usually the local part of a distributed vector.
     *  - localX is usually the local part of a distributed vector.
     */
    void haloOperationSync(
        hmemo::HArray<ValueType>& localResult,
        const hmemo::HArray<ValueType>& localX,
        hmemo::HArray<ValueType>& haloX,
        common::function <
        void(
            const MatrixStorage<ValueType>* localMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& localX ) > localF,
        common::function <
        void(
            const MatrixStorage<ValueType>* haloMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& haloX ) > haloF ) const;

    /**
     * @brief Operation on transposed distributed matrix with halo exchange, sync version
     *
     * @param[out]    localResult is the result array
     * @param[in]     localX is the array with local values
     * @param[in,out] haloX is a temporary array keeping the non-local halo values of X
     * @param[in]     localF is the operation called with local matrix and localX
     * @param[in]     haloF is the operation called with halo matrix and haloX
     *
     * Very similiar to haloOperationSync but here the halo computations are done for
     * the other processors (as this processor owns it) and the results are sent to the other
     * processors. Uses inverse halo communication schedule, scatters the received values instead
     * of gathering the send values.
     */
    void invHaloOperationSync(
        hmemo::HArray<ValueType>& localResult,
        const hmemo::HArray<ValueType>& localX,
        hmemo::HArray<ValueType>& haloX,
        common::function <
        void(
            const MatrixStorage<ValueType>* localMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& localX ) > localF,
        common::function <
        void(
            const MatrixStorage<ValueType>* haloMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& haloX ) > haloF ) const;

    /**
     * @brief Operation on distributed matrix with halo exchange, async version
     *
     * The asynchronous version starts the local computation asynchronously so it
     * can overlap with the halo exchange.
     */
    void haloOperationAsync(
        hmemo::HArray<ValueType>& localResult,
        const hmemo::HArray<ValueType>& localX,
        hmemo::HArray<ValueType>& haloX,
        common::function <
        tasking::SyncToken * (
            const MatrixStorage<ValueType>* localMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& localX ) > localAsyncF,
        common::function <
        void(
            const MatrixStorage<ValueType>* haloMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& haloX ) > haloF ) const;

    /* Implemenation of pure method of class Matrix */

    virtual void matrixPlusMatrix( const Scalar alpha, const Matrix& A, const Scalar beta, const Matrix& B );

    /**
     * Override Matrix::cat
     */
    virtual void cat( const IndexType dim, const Matrix* other[], const IndexType n );

    /* Implemenation of pure method of class Matrix */

    virtual void matrixTimesMatrix(
        Matrix& result,
        const Scalar alpha,
        const Matrix& B,
        const Scalar beta,
        const Matrix& C ) const;

    /* Implementation of pure method of class Matrix. */
    virtual Scalar l1Norm() const;

    /* Implementation of pure method of class Matrix. */
    virtual Scalar l2Norm() const;

    /** Implementation of pure method of class Matrix for sparse matrices. */

    virtual Scalar maxNorm() const;

    /** Implementation of pure method of class Matrix for sparse matrices. */

    virtual Scalar maxDiffNorm( const Matrix& other ) const;

    /**
     * @brief Same as maxDiffNorm but with other as sparse matrix of same value type.
     */
    ValueType maxDiffNormImpl( const SparseMatrix<ValueType>& other ) const;

    /* Implementation of pure method of class Matrix. */

    IndexType getLocalNumValues() const;

    /* Implementation of pure method of class Matrix. */

    IndexType getLocalNumRows() const;

    /* Implementation of pure method of class Matrix. */

    IndexType getLocalNumColumns() const;

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getNumValues() const;

    /* Getter routine returning the num values of the local partition (local + halo) of the matrix */

    IndexType getPartitialNumValues() const;

    /* Implementation of pure method of class Matrix. */

    virtual Scalar getValue( IndexType i, IndexType j ) const;

    /** Implementation of pure method Matrix::setValue */

    virtual void setValue(
        const IndexType i,
        const IndexType j,
        const Scalar val,
        const common::binary::BinaryOp op = common::binary::COPY );

    /**
     * @brief Read access to the halo of the distributed matrix.
     *
     * @return   reference to the halo of the distributed matrix
     */
    const dmemo::Halo& getHalo() const;

    /* Implementation of method writeAt for sparse matrix. */

    virtual void writeAt( std::ostream& stream ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void prefetch() const;

    /* Implementation of pure method of class Matrix. */

    virtual void wait() const;

    /* Implementation of pure method of class Matrix. */

    virtual bool hasDiagonalProperty() const;

    /* Implementation of pure method of class Matrix. */

    virtual void resetDiagonalProperty();

    /* Implementation of pure method of class Matrix. */

    virtual common::scalar::ScalarType getValueType() const;

    virtual size_t getValueTypeSize() const;

    SparseMatrix<ValueType>* newMatrix() const;

    /* Implementation of pure method Matrix::copy with covariant return type */

    virtual SparseMatrix<ValueType>* copy() const;

    /* Implementation of pure method of class Matrix. */

    virtual void redistribute( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /* Implementation of pure method of class Matrix. */

    virtual void redistribute( const dmemo::Redistributor& redistributor, dmemo::DistributionPtr colDistribution );

    /**  */
    /**
     * @brief Assign another matrix transposed to this matrix.
     *
     * @param[in] matrix   input matrix that will be transposed
     */
    void assignTranspose( const Matrix& matrix );

    /* Implementation of pure method of class Matrix. */

    virtual size_t getMemoryUsage() const;

    using Matrix::operator=; // make overloaded routines visible before overwriting one

    using Matrix::getColDistribution;
    using Matrix::getColDistributionPtr;
    using Matrix::getRowDistribution;
    using Matrix::getRowDistributionPtr;
    using Matrix::setDistributionPtr;

    using Matrix::getCommunicationKind;

    using Matrix::getNumColumns;
    using Matrix::getNumRows;

    using Matrix::redistribute;

    /** Override the default assignment operator to guarantee deep copy. */

    SparseMatrix<ValueType>& operator=( const SparseMatrix<ValueType>& matrix );

    /**
     * Gives info about the matrix kind (SPARSE).
     */
    virtual Matrix::MatrixKind getMatrixKind() const
    {
        return Matrix::SPARSE;
    }

    /** Get a complete row of this matrix from its local part. */

    void getLocalRowDense( hmemo::HArray<ValueType>& row, const IndexType localRowIndex ) const;

    /** Get a complete row of this matrix from its local part in sparse format */

    void getLocalRowSparse( hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values, const IndexType localRowIndex ) const;

    /** Implementation of pure method Matrix::getRow */

    virtual void getRow( Vector& row, const IndexType globalRowIndex ) const;

    /** Implementation of pure method Matrix::getRowLocal */

    virtual void getRowLocal( Vector& row, const IndexType localRowIndex ) const;

    /** Set a complete row of this matrix in its local part. */

    void setLocalRow( const hmemo::HArray<ValueType>& row,
                      const IndexType localRowIndex,
                      const common::binary::BinaryOp op  );

    /** Get the local part of a col of this matrix */

    void getLocalColumn( hmemo::HArray<ValueType>& col, const IndexType colIndex ) const;

    /** Set the local col of this matrix */

    void setLocalColumn( const hmemo::HArray<ValueType>& column,
                         const IndexType colIndex,
                         const common::binary::BinaryOp op  );

protected:

    /** Test consistency of sparse matrix data, only used if ASSERT_DEBUG is enabled. */

    void checkSettings();

    common::shared_ptr<MatrixStorage<ValueType> > mLocalData; //!< local columns of sparse matrix

    common::shared_ptr<MatrixStorage<ValueType> > mHaloData; //!< local columns of sparse matrix

    dmemo::Halo mHalo; //!< Exchange plans for halo part due to column distribution

    /**
     * @brief Set this matrix = alpha * A + beta * B
     *
     * @param[in]  alpha    scaling of input matrix A
     * @param[in]  A        input matrix
     * @param[in]  beta     scaling of input matrix B
     * @param[in]  B        input matrix
     */
    void matrixPlusMatrixImpl(
        const ValueType alpha,
        const SparseMatrix<ValueType>& A,
        const ValueType beta,
        const SparseMatrix<ValueType>& B );
    /**
     * @brief Set this matrix = alpha * A * B + beta * C
     *
     * @param[in]  alpha    scaling of matrix product
     * @param[in]  A        first matrix of matrix product
     * @param[in]  B        second matrix of matrix product
     * @param[in]  beta     scaling of matrix summand
     * @param[in]  C        matrix to be added to the product
     */
    void matrixTimesMatrixImpl(
        const ValueType alpha,
        const SparseMatrix<ValueType>& A,
        const SparseMatrix<ValueType>& B,
        const ValueType beta,
        const SparseMatrix<ValueType>& C );

    using Matrix::mColDistribution;

public:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

protected:

    /**
     * @brief Default constructor is disabled.
     */
    SparseMatrix();

private:

    /** This method sets a row-distributed matrix corresponding to the distribution of this matrix.
     *  ( no column distribution, no halo ).
     *
     *  @param[in]  otherLocalData  is the local part of the other matrix owned by this partition
     *  @param[in]  otherDist is the row distribution of the other matrix.
     *
     *  The new row / column distribution of this matrix is already set as member variables.
     *  This routine can also handle the case that otherLocalData is a reference to the local
     *  data of this matrix ( helpful to avoid unneccessary copies ).
     */
    void set( const MatrixStorage<ValueType>& otherLocalData, dmemo::DistributionPtr otherDist );

    /** Implementation of transposed assign for sparse matrix of a known value type. */

    void assignTransposeImpl ( const SparseMatrix<ValueType>& matrix );

    static std::string initTypeName();

    mutable hmemo::HArray<ValueType> mTempSendValues; //!< temporary vector for halo communications

};

/***************************************************************************************************/
/* Implementation of inline methods                                                                */
/***************************************************************************************************/

template <typename ValueType>
void SparseMatrix<ValueType>::setContextPtr( const hmemo::ContextPtr context )
{
    SCAI_ASSERT_ERROR( context, "NULL context, cannot be set" )

    mLocalData->setContextPtr( context );
    mHaloData->setContextPtr( context );
}

/* Implementation of pure method of class Matrix. */

template <typename ValueType>
hmemo::ContextPtr SparseMatrix<ValueType>::getContextPtr() const
{
    return mLocalData->getContextPtr();
}

} /* end namespace lama */

} /* end namespace scai */
