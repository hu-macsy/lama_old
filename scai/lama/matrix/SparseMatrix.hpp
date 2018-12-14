/**
 * @file SparseMatrix.hpp
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
 * @brief Definition of class for distributed sparse matrices.
 * @author Jiri Kraus, Thomas Brandes
 * @date 06.06.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/Matrix.hpp>

// local library
#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/HaloPlan.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <functional>

namespace scai
{

namespace lama
{

// forward declarations

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

    public Matrix<ValueType>

{

public:

    /* Using clauses for convenience, avoids using this->... */

    using _Matrix::operator=;
    using _Matrix::setContextPtr; 
    using _Matrix::getNumRows;
    using _Matrix::getNumColumns;

    using _Matrix::getRowDistribution;
    using _Matrix::getRowDistributionPtr;
    using _Matrix::getColDistribution;
    using _Matrix::getColDistributionPtr;

    using _Matrix::redistribute;

    using Matrix<ValueType>::getValueType;
    using Matrix<ValueType>::operator=;
    using Matrix<ValueType>::operator-=;
    using Matrix<ValueType>::operator+=;

    typedef ValueType MatrixValueType; //!< This is the type of the matrix values.

    /** Getter for the type name of the class. */

    static const char* typeName();

    /* Implementation of pure method of class _Matrix. */

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
    virtual Format getFormat() const
    {
        return mLocalData->getFormat();
    }

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
    SparseMatrix( std::shared_ptr<MatrixStorage<ValueType> > storage );

    /**
     *  @brief Constructor of a distributed sparse matrix with local storages
     */
    SparseMatrix( dmemo::DistributionPtr rowDist, std::shared_ptr<MatrixStorage<ValueType> > storage );

    /** Override also the default copy constructor that does not make a
     *  deep copy of the input matrix due to the use of shared pointers.
     */
    SparseMatrix( const SparseMatrix<ValueType>& other );

    /** Move constructor */

    SparseMatrix( SparseMatrix<ValueType>&& other ) noexcept;

    /** Implementation of abstract method for sparse matrices. */

    virtual bool isConsistent() const;

    /* Implementation of pure method of class _Matrix. */

    virtual void invert( const _Matrix& other );

    /* Implementation of pure method of class _Matrix. */

    virtual void setContextPtr( const hmemo::ContextPtr context );

    /* Implementation of pure method of class _Matrix. */

    virtual hmemo::ContextPtr getContextPtr() const;

    /* Implementation of pure method of class _Matrix. */

    virtual void clear();

    /* Implementation of pure method of class _Matrix. */

    virtual void purge();

    /* Implementation of pure method of class _Matrix. */

    virtual void allocate( const IndexType numRows, const IndexType numColumns );

    /* Implementation of pure method of class _Matrix. */

    virtual void allocate( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /**
     * @brief This method removes all zero elements of a sparse matrix, i.e. only entries whose absolute
     *        value is greater than eps are considered to be non-zero.
     *
     * This operation should be executed by all processors as it might recompute the halo.
     *
     * \code
     *    auto diffMatrix = eval<CSRSparseMatrix<ValueType>>( matrix1 - matrix2 );
     *    diffMatrix.compress( 0.0001 );
     * \endcode
     */
    void compress( const RealType<ValueType> eps = 0 );

    /** @brief Implementation of pure method Matrix<ValueType>::getColumn 
     *
     *  It is recommended to call getColumn with a SparseVector for a sparse matrix.
     *  It works also with a dense vector but setting the zero values causes a significant overhead.
     */
    virtual void getColumn( Vector<ValueType>& column, const IndexType globalColIndex ) const;

    /** Set matrix to a identity square matrix with same row and column distribution. */

    virtual void setIdentity( dmemo::DistributionPtr distribution );

    /** Implementation of pure method Matrix<ValueType>::assignDiagonal */

    virtual void assignDiagonal( const Vector<ValueType>& diagonal );

    /* Implementation of pure method of class _Matrix. */

    virtual void assign( const _Matrix& other );

    /** Method that assigns a sparse matrix, specialization of assign( const _Matrix& ) */

    void assign( const SparseMatrix<ValueType>& matrix );

    /* Implementation of pure method of class _Matrix. */

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

    /* Implementation of pure method of class _Matrix. */

    virtual void assignDistribute( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    virtual void assignLocal( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist );

    virtual void assignDistribute( const _Matrix& other, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Implementation of of pure method of class _Matrix. */

    virtual void buildLocalStorage( _MatrixStorage& storage ) const;

    /** Implemenation of Matrix<ValueType>::disassemble */

    virtual void disassemble(
        MatrixAssembly<ValueType>& assembly,
        const IndexType rowOffset = 0,
        const IndexType colOffset = 0 ) const;

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

    /** Implementation of pure method Matrix<ValueType>::getDiagonal */

    virtual void getDiagonal( Vector<ValueType>& diagonal ) const;

    /** Implementation of pure method Matrix<ValueType>::setDiagonal */

    virtual void setDiagonal( const Vector<ValueType>& diagonal );

    /** Implementation of pure method Matrix<ValueType>::setDiagonal */

    virtual void setDiagonal( const ValueType& scalar );

    /* Implementation of pure method Matrix<ValueType>::reduce */

    virtual void reduce( 
        Vector<ValueType>& v, 
        const IndexType dim, 
        const common::BinaryOp reduceOp, 
        const common::UnaryOp elemOp ) const;

    /* Implementation of reduce for dense vector. */

    void reduceImpl(
        DenseVector<ValueType>& v,
        const IndexType dim,
        const common::BinaryOp reduceOp,
        const common::UnaryOp elemOp ) const;

    /* ======================================================================= */
    /*     scaling of matrix entries                                           */
    /* ======================================================================= */

    /* Implementation of pure method Matrix<ValueType>::scale */

    virtual void scale( const ValueType& alpha );

    /* Implementation of pure method Matrix<ValueType>::scaleRows */

    virtual void scaleRows( const DenseVector<ValueType>& scaleY );

    /* Implementation of pure method Matrix<ValueType>::scaleColumns */

    virtual void scaleColumns( const DenseVector<ValueType>& scaleY );

    /* Implementation of pure method of class _Matrix. */

    virtual void conj();

    /* Implemenation of pure method Matrix<ValueType>::matrixTimesScalar */

    virtual void matrixTimesScalar( const Matrix<ValueType>& other, const ValueType alpha );

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
        const DenseVector<ValueType>* y ) const;

    /**
     * @brief Same as matrixTimesVector but with vectors result, x, and y of same value type.
     */
    void matrixTimesVectorImpl(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>* y ) const;

    /**
     * @brief Implementation of pure method Matrix<ValueType>::matrixTimesVectorDense
     */
    void matrixTimesVectorDense(
        DenseVector<ValueType>& result,
        const ValueType alpha,
        const DenseVector<ValueType>& x,
        const ValueType beta,
        const DenseVector<ValueType>* y,
        const common::MatrixOp op ) const;

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
        std::function <
        void(
            const MatrixStorage<ValueType>* localMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& localX ) > localF,
        std::function <
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
        std::function <
        void(
            const MatrixStorage<ValueType>* localMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& localX ) > localF,
        std::function <
        void(
            const MatrixStorage<ValueType>* haloMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& haloX ) > haloF ) const;

    /*
     * @brief Operation on distributed matrix with halo exchange, asynchronous communication
     *
     * This asynchronous version starts the communication asynchronously so it
     * can overlap with the halo exchange.
     */
    void haloOperationAsyncComm(
        hmemo::HArray<ValueType>& localResult,
        const hmemo::HArray<ValueType>& localX,
        hmemo::HArray<ValueType>& haloX,
        std::function <
        void(
            const MatrixStorage<ValueType>* localMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& localX ) > localF,
        std::function <
        void(
            const MatrixStorage<ValueType>* haloMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& haloX ) > haloF ) const;

    /**
     * @brief Operation on distributed matrix with halo exchange, async execution of local computation
     *
     * This asynchronous version starts the local computation asynchronously so it
     * can overlap with the halo exchange.
     */
    void haloOperationAsyncLocal(
        hmemo::HArray<ValueType>& localResult,
        const hmemo::HArray<ValueType>& localX,
        hmemo::HArray<ValueType>& haloX,
        std::function <
        tasking::SyncToken * (
            const MatrixStorage<ValueType>* localMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& localX ) > localAsyncF,
        std::function <
        void(
            const MatrixStorage<ValueType>* haloMatrix,
            hmemo::HArray<ValueType>& localResult,
            const hmemo::HArray<ValueType>& haloX ) > haloF ) const;

    /* Implementation of pure method Matrix<ValueType>::binaryOp */

    virtual void binaryOp( const Matrix<ValueType>& matrixA, const common::BinaryOp op, const Matrix<ValueType>& matrixB );

    /* Implemenation of pure method Matrix<ValueType>::matrixPlusMatrix */

    virtual void matrixPlusMatrix( const ValueType alpha, const Matrix<ValueType>& A, const ValueType beta, const Matrix<ValueType>& B );

    /* Implemenation of pure method Matrix<ValueType>::matrixTimesMatrix */

    virtual void matrixTimesMatrix(
        Matrix<ValueType>& result,
        const ValueType alpha,
        const Matrix<ValueType>& B,
        const ValueType beta,
        const Matrix<ValueType>& C ) const;

    /* Implementation of pure method of class Matrix. */

    virtual RealType<ValueType> l1Norm() const;

    /* Implementation of pure method of class _Matrix. */
    virtual RealType<ValueType> l2Norm() const;

    /** Implementation of pure method of class _Matrix for sparse matrices. */

    virtual RealType<ValueType> maxNorm() const;

    /** Implementation of pure method of class _Matrix for sparse matrices. */

    virtual RealType<ValueType> maxDiffNorm( const Matrix<ValueType>& other ) const;

    /**
     * @brief Same as maxDiffNorm but with other as sparse matrix of same value type.
     */
    RealType<ValueType> maxDiffNormImpl( const SparseMatrix<ValueType>& other ) const;

    /* Implementation of pure method of class _Matrix. */

    IndexType getLocalNumValues() const;

    /* Implementation of pure method of class _Matrix. */

    IndexType getLocalNumRows() const;

    /* Implementation of pure method of class _Matrix. */

    IndexType getLocalNumColumns() const;

    /* Implementation of pure method of class _Matrix. */

    virtual IndexType getNumValues() const;

    /* Getter routine returning the num values of the local partition (local + halo) of the matrix */

    IndexType getPartitialNumValues() const;

    /* Implementation of pure method of class _Matrix. */

    virtual ValueType getValue( IndexType i, IndexType j ) const;

    /** Implementation of pure method _Matrix::setValue */

    virtual void setValue(
        const IndexType i,
        const IndexType j,
        const ValueType val,
        const common::BinaryOp op = common::BinaryOp::COPY );

    /**
     * @brief Read access to the halo of the distributed matrix.
     *
     * @return   reference to the halo of the distributed matrix
     */
    const dmemo::HaloPlan& getHaloPlan() const;

    /* Implementation of method writeAt for sparse matrix. */

    virtual void writeAt( std::ostream& stream ) const;

    /* Implementation of pure method of class _Matrix. */

    virtual void prefetch() const;

    /* Implementation of pure method of class _Matrix. */

    virtual void wait() const;

    /* Implementation of pure method _Matrix::newMatrix with covariant return type */

    SparseMatrix<ValueType>* newMatrix() const;

    /* Implementation of pure method _Matrix::copy with covariant return type */

    virtual SparseMatrix<ValueType>* copy() const;

    /* Implementation of pure method _Matrix::redistribute */

    virtual void redistribute( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /* Implementation of pure method of _Matrix::redistribute */

    virtual void redistribute( const dmemo::Redistributor& redistributor, dmemo::DistributionPtr colDistribution );

    /* Implementation of pure method _Matrix::resize */

    virtual void resize( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /**
     * @brief Assign another matrix transposed to this matrix.
     *
     * @param[in] matrix   input matrix that will be transposed
     */
    void assignTranspose( const _Matrix& matrix );

    /* Implementation of pure method of class _Matrix. */

    virtual size_t getMemoryUsage() const;

    /** Override the default assignment operator to guarantee deep copy. */

    SparseMatrix<ValueType>& operator=( const SparseMatrix<ValueType>& matrix );

    SparseMatrix<ValueType>& operator=( SparseMatrix<ValueType>&& matrix );

    /**
     * Gives info about the matrix kind (SPARSE).
     */
    virtual MatrixKind getMatrixKind() const
    {
        return MatrixKind::SPARSE;
    }

    /** Get a complete row of this matrix from its local part. */

    void getLocalRowDense( hmemo::HArray<ValueType>& row, const IndexType localRowIndex ) const;

    /** Get a complete row of this matrix from its local part in sparse format */

    void getLocalRowSparse( hmemo::HArray<IndexType>& indexes, hmemo::HArray<ValueType>& values, const IndexType localRowIndex ) const;

    /** Implementation of pure method Matrix<ValueType>::getRow */

    virtual void getRow( Vector<ValueType>& row, const IndexType globalRowIndex ) const;

    /** Implementation of pure method Matrix<ValueType>::getRowLocal */

    virtual void getRowLocal( Vector<ValueType>& row, const IndexType localRowIndex ) const;

    /** Set a complete row of this matrix in its local part. */

    void setLocalRow( const hmemo::HArray<ValueType>& row,
                      const IndexType localRowIndex,
                      const common::BinaryOp op  );

    /** Get the local part of a col of this matrix */

    void getLocalColumn( hmemo::HArray<ValueType>& col, const IndexType colIndex ) const;

    /** Set the local col of this matrix */

    void setLocalColumn( const hmemo::HArray<ValueType>& column,
                         const IndexType colIndex,
                         const common::BinaryOp op  );

protected:

    /** Test consistency of sparse matrix data, only used if ASSERT_DEBUG is enabled. */

    void checkSettings();

    std::shared_ptr<MatrixStorage<ValueType> > mLocalData; //!< local columns of sparse matrix

    std::shared_ptr<MatrixStorage<ValueType> > mHaloData; //!< local columns of sparse matrix

    dmemo::HaloPlan mHaloPlan; //!< Exchange plans for halo part due to column distribution

    /**
     * @brief Set this matrix = alpha * A + beta * B
     *
     * @param[in]  alpha    scaling of input matrix A
     * @param[in]  A        input matrix
     * @param[in]  beta     scaling of input matrix B
     * @param[in]  B        input matrix
     */
    void matrixPlusMatrixSparse(
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

    /**
     *  @brief element-wise binary operation for two sparse matrices.
     */
    virtual void binaryOpSparse( 
        const SparseMatrix<ValueType>& matrixA, 
        const common::BinaryOp op, 
        const SparseMatrix<ValueType>& matrixB );

    /** Implementation of pure method Matrix<ValueType>::selectComplexPart */

    virtual void selectComplexPart( Matrix<RealType<ValueType> >& x, common::ComplexPart kind ) const;

    /** Implementation of pure method Matrix<ValueType>::buildComplex */

    virtual void buildComplex( const Matrix<RealType<ValueType> >& x, const Matrix<RealType<ValueType> >& y );

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

/* Implementation of pure method of class _Matrix. */

template <typename ValueType>
hmemo::ContextPtr SparseMatrix<ValueType>::getContextPtr() const
{
    return mLocalData->getContextPtr();
}

} /* end namespace lama */

} /* end namespace scai */
