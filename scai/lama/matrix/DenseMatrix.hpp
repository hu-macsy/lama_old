/**
 * @file DenseMatrix.hpp
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
 * @brief Definition of matrix class for distributed matrixes in Dense format.
 * @author Michael Drost
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/Matrix.hpp>

// local library
#include <scai/lama/matrix/SparseMatrix.hpp>

#include <scai/lama/storage/DenseStorage.hpp>

#include <memory>

namespace scai
{

namespace lama
{

template<typename ValueType> class DenseVector;
// forward declaration

/** Class for dense matrices where rows are distributed among rows and columns
 *  are splitted according to a column distribution.
 *
 *  The local rows are splitted according to a column distribution into different
 *  blocks that are stored in a vector of shared pointers. A copy of this vector
 *  will not make deep copies of the blocks so default copy constructor and
 *  assignment operator must be overridden.
 *
 *  @tparam ValueType is the value type of the matrix values.
 */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT DenseMatrix:

    public Matrix<ValueType>,
    public _Matrix::Register<DenseMatrix<ValueType> >    // register at factory
{

public:

    /* Using clauses for convenience, avoids using this->... */

    using _Matrix::operator=;
    using _Matrix::setContextPtr; 
    using _Matrix::getNumRows;
    using _Matrix::getNumColumns;
    using _Matrix::setIdentity;   

    using _Matrix::getRowDistribution;
    using _Matrix::getRowDistributionPtr;
    using _Matrix::getColDistribution;
    using _Matrix::getColDistributionPtr;

    using _Matrix::redistribute;

    using Matrix<ValueType>::operator=;
    using Matrix<ValueType>::operator+=;
    using Matrix<ValueType>::operator-=;
    using Matrix<ValueType>::getValueType;

    typedef ValueType MatrixValueType; //!< This is the type of the matrix values.

    typedef DenseStorage<ValueType> StorageType;

    typedef std::shared_ptr<DenseStorage<ValueType> > DenseStoragePtr;

    /** Getter for the type name of the class. */

    static const char* typeName();

    /** Default constructor. */

    DenseMatrix( hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Constructor of a replicated dense matrix.
     *
     * @param[in] numRows      number of rows, must be non-negative.
     * @param[in] numColumns   number of columns, must be non-negative.
     * @param[in] ctx          context for the new matrix.
     */
    DenseMatrix( const IndexType numRows, const IndexType numColumns, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /**
     * Constructor of a distributed dense matrix.
     *
     * @param[in] rowDist   size and distribution of rows
     * @param[in] colDist   size and distribution of columns
     * @param[in] ctx       context for memory/operations of this matrix
     *
     * For consistency with the constructors of sparse matrices the values
     * of the dense matrix are initialized with 0 here.
     */
    DenseMatrix( dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist, hmemo::ContextPtr ctx = hmemo::Context::getContextPtr() );

    /** Overwrites default copy constructor so it uses other copy constructor.
     *
     *  Note: Default copy constructor would not make deep copies of the
     *        storages, so it must be overridden
     */
    DenseMatrix( const DenseMatrix<ValueType>& other );

    /** Override default move constructor */

    DenseMatrix( DenseMatrix<ValueType>&& other ) noexcept;

    /** Constructor of a (replicated) sparse matrix by global storage.
     *
     *  @param[in] globalStorage  contains the full storage, must be of same format and type as matrix
     */
    explicit DenseMatrix( DenseStorage<ValueType> globalStorage );

    /** Constructor of a distributed dense matrix by local storage
     *
     *  @param[in] rowDist       is distribution of localData
     *  @param[in] localStorage  contains local rows of the distributed matrix
     *
     *  The number of rows for the local storage must be rowDist->getLocalSize(), and the 
     *  number of columns must be the same on all processors.
     *
     *  The context of the matrix is the same as the context of localStorage.
     */
    DenseMatrix( dmemo::DistributionPtr rowDist, DenseStorage<ValueType> localStorage );

    /**
     * Destructor, releases all allocated resources.
     */
    virtual ~DenseMatrix();

    /** Implementation of abstract method for dense matrices. */

    virtual bool isConsistent() const;

    /** Overrides the default assignment operator to guarantee deep copy. */

    DenseMatrix& operator=( const DenseMatrix<ValueType>& matrix );

    /** Override implicitly declared move assignment operator */

    DenseMatrix& operator=( DenseMatrix<ValueType>&& matrix );

    /** Implementation for _Matrix::getTypeName() */

    const char* getTypeName() const;

    /**
     * Gives info about the matrix kind (DENSE).
     */
    virtual MatrixKind getMatrixKind() const
    {
        return MatrixKind::DENSE;
    }

    /* Implementation of pure method of class _Matrix. */

    virtual void setContextPtr( const hmemo::ContextPtr context );

    /* Implementation of pure method of class _Matrix. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        if ( mData.size() >  0 )
        {
            return mData[0]->getContextPtr();
        }
        else
        {
            return hmemo::Context::getHostPtr();
        }
    }

    /** Implementation of pure method _Matrix::setIdentity. */

    virtual void setIdentity( dmemo::DistributionPtr distribution );

    /** Implementation of pure method Matrix<ValueType>::assignDiagonal */

    virtual void assignDiagonal( const Vector<ValueType>& diagonal );

    /** Implementation for pure method _Matrix::setCSRData. */

    virtual void setCSRData(
        dmemo::DistributionPtr rowDist,
        dmemo::DistributionPtr colDist,
        const IndexType numValues,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::_HArray& values );

    /** Implementation of pure method for the dense storage format. */

    virtual void buildCSRData( hmemo::HArray<IndexType>& rowIA, hmemo::HArray<IndexType>& rowJA, hmemo::_HArray& rowValues ) const;

    /** Implementation of pure method. */

    virtual void setCSRData(
        const hmemo::HArray<IndexType>& rowIA,
        const hmemo::HArray<IndexType>& rowJA,
        const hmemo::_HArray& rowValues,
        dmemo::DistributionPtr rowDistribution,
        dmemo::DistributionPtr colDistribution );

    /** Local version of setCSRData . */

    void setCSRDataLocal(
        const hmemo::HArray<IndexType>& rowIA,
        const hmemo::HArray<IndexType>& rowJA,
        const hmemo::_HArray& rowValues ) const;

    /* Implementation of pure method of class _Matrix. */

    virtual void clear();

    /* Implementation of pure method _Matrix::purge. */

    virtual void purge();

    /* Implementation of pure method of class _Matrix. */

    virtual void allocate( const IndexType numRows, const IndexType numColumns );

    /* Implementation of pure method of class _Matrix. */

    virtual void allocate( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /* Implementation of pure method of class _Matrix. */

    virtual void assign( const _Matrix& other );

    /* Implementation of pure method of class _Matrix. */

    virtual void assignTranspose( const _Matrix& other );

    void assignTransposeImpl( const DenseMatrix<ValueType>& Mat );

    /** @brief Swap will swap all member variables of the two dense matrices.
     *
     *  This operation might be useful in iteration loops where a dense matrix
     *  is updated each iteration. It is more convenient than a solution that
     *  is based on using pointers in the application.
     *
     * @param[in] other   TODO[doxy] Complete Description.
     */
    void swap( DenseMatrix<ValueType>& other );

    /** Method that assigns a sparse matrix, specialization of assign( const _Matrix& ) */

    template<typename OtherValueType>
    void assignSparse( const SparseMatrix<OtherValueType>& other );

    /* Implementation of pure method of class _Matrix. */

    virtual void assign( const _MatrixStorage& storage );

    template<typename OtherValueType>
    void assignImpl( const Matrix<OtherValueType>& other );

    /* Implementation of pure method of class _Matrix. */

    virtual void assignDistribute( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    virtual void assignLocal( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist );

    virtual void assignDistribute( const _Matrix& other, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** @brief TODO[doxy] Complete Description.
     *
     * @param[in] other   TODO[doxy] Complete Description.
     */
    void assignLocal( const _MatrixStorage& other );

    /** Implementation of _Matrix::buildLocalStorage. */

    virtual void buildLocalStorage( _MatrixStorage& storage ) const;

    /** Implemenation of Matrix<ValueType>::disassemble */

    virtual void disassemble(
        MatrixAssembly<ValueType>& assembly,
        const IndexType rowOffset = 0,
        const IndexType colOffset = 0 ) const;

    /* Implementation of pure method of class _Matrix. */

    virtual void redistribute( dmemo::DistributionPtr rowDistributionPtr, dmemo::DistributionPtr colDistributionPtr );

    /* Implementation of pure method of class _Matrix. */

    virtual void redistribute( const dmemo::Redistributor& redistributor, dmemo::DistributionPtr colDistributionPtr );

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

    /* Implementation of pure method Matrix<ValueType>::scale */

    virtual void scaleRows( const DenseVector<ValueType>& scaleY );

    /* Implementation of pure method of class _Matrix. */

    virtual void conj();

    /* Implementation of pure method of class _Matrix. */

    virtual ValueType getValue( IndexType i, IndexType j ) const;

    /** Implementation of pure method _Matrix::setValue */

    virtual void setValue(
        const IndexType i,
        const IndexType j,
        const ValueType val,
        const common::BinaryOp op = common::BinaryOp::COPY );

    /* Implemenation of pure method Matrix<ValueType>::matrixTimesScalar */

    virtual void matrixTimesScalar( const Matrix<ValueType>& other, const ValueType alpha );

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
     *  @brief Provide method matrixTimesVector where all vectors are now dense
     *
     * @param[out] denseResult   result vector 
     * @param[in]  alphaValue    scaling factor for matrix * vector
     * @param[in]  denseX        vector that is used for multiplication
     * @param[in]  betaValue     scaling factor for additional summand
     * @param[in]  denseY        additional summand ( beta = 0 if not available )
     *
     *  Note: _Matrix::matrixTimesMatrix is implemented in the CRTPMatrix class.
     *        that requires this method.
     *
     *  Note: all vectors must have the right distribution.
     */
    void matrixTimesVectorImpl(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>* denseY ) const;

    void vectorTimesMatrixImpl(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>* denseY ) const;

    /* Implementation of pure method of class _Matrix. */

    virtual void invert( const _Matrix& other );

    /** Invert in place */

    void invert()
    {
        this->invert( *this );
    }

    /* Implementation of pure method of class _Matrix. */
    virtual RealType<ValueType> l1Norm() const;

    /* Implementation of pure method of class _Matrix. */
    virtual RealType<ValueType> l2Norm() const;

    /** Implementation of _Matrix::maxNorm for dense matrices. */

    virtual RealType<ValueType> maxNorm() const;

    /** Implementation of _Matrix::maxDiffNorm for dense matrices. */

    virtual RealType<ValueType> maxDiffNorm( const Matrix<ValueType>& other ) const;

    /** Get the maximal difference between two elements for dense matrices of same type. */

    ValueType maxDiffNormImpl( const DenseMatrix<ValueType>& other ) const;

    /* Implemenation of pure method Matrix<ValueType>::binaryOp */

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

    /* Implementation of pure method of class _Matrix. */

    virtual void prefetch() const;

    /** @brief TODO[doxy] Complete Description.
     *
     * @param[in] loc   TODO[doxy] Complete Description.
     */
    void prefetch( hmemo::ContextPtr loc ) const;

    /* Implementation of pure method of class _Matrix. */

    void wait() const;

    /** This method returns the storage containing the local data regarding row/col distribution. */

    const DenseStorage<ValueType>& getLocalStorage() const;

    DenseStorage<ValueType>& getLocalStorage();

    /* Implementation of pure method of class _Matrix. */

    virtual IndexType getLocalNumValues() const;

    /* Implementation of pure method of class _Matrix. */

    virtual IndexType getLocalNumRows() const;

    /* Implementation of pure method of class _Matrix. */

    virtual IndexType getLocalNumColumns() const;

    /* Implementation of pure method of class _Matrix. */

    virtual IndexType getNumValues() const;

    /* Implementation of method writeAt for dense matrix. */

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Implementation of pure function _Matrix::copy with covariant return type.
     */
    virtual DenseMatrix<ValueType>* newMatrix() const;

    /**
     * @brief Implementation of pure function _Matrix::copy with covariant return type.
     */
    virtual DenseMatrix<ValueType>* copy() const;

    /* Implementation of pure method _Matrix::getFormat */

    virtual Format getFormat() const
    {
        return Format::DENSE;
    }

    /* Implementation of pure method of class _Matrix. */

    virtual size_t getMemoryUsage() const;

    /** local data is allocated in chunks according to column distribution */

    std::vector<std::shared_ptr<DenseStorage<ValueType> > > mData;

    /** Implementation of pure methode Matrix<ValueType>::getRow */

    virtual void getRow( Vector<ValueType>& row, const IndexType globalRowIndex ) const;

    /** Implementation of pure methode Matrix<ValueType>::getRowLocal */

    virtual void getRowLocal( Vector<ValueType>& row, const IndexType globalRowIndex ) const;

    /** Implementation of pure method Matrix<ValueType>::getColumn */

    virtual void getColumn( Vector<ValueType>& col, const IndexType globalColIndex ) const;

    /** Get a complete row of the local storage */

    void getLocalRow( hmemo::HArray<ValueType>& row, const IndexType iLocal ) const;

    void setLocalRow( const hmemo::HArray<ValueType>& row,
                      const IndexType localRowIndex,
                      const common::BinaryOp op  );

    void setLocalColumn( const hmemo::HArray<ValueType>& column,
                         const IndexType colIndex,
                         const common::BinaryOp op  );

    /** Copy a dense matrix with different data type; inherits sizes and distributions */

    template<typename otherT>
    void assignDense( const DenseMatrix<otherT>& other );

    /** Optimized implementation for dense vectors as diagonal. */

    template<typename OtherT>
    void getDiagonalImpl( DenseVector<OtherT>& diagonal ) const;

protected:

    /**
     * @brief Same as matrixPlusMatrix but now input matrices A and B are really dense
     */
    void matrixPlusMatrixDense(
        const ValueType alpha,
        const DenseMatrix<ValueType>& A,
        const ValueType beta,
        const DenseMatrix<ValueType>& B );

    /** Implementation of pure method Matrix<ValueType>::selectComplexPart */

    virtual void selectComplexPart( Matrix<RealType<ValueType> >& x, common::ComplexPart kind ) const;

    /** Implementation of pure method Matrix<ValueType>::buildComplex */

    virtual void buildComplex( const Matrix<RealType<ValueType> >& x, const Matrix<RealType<ValueType> >& y );

private:

    /** Allocate of storage for the column blocks. */

    void allocateData( hmemo::ContextPtr ctx );

    /** Join column data of column distributed dense data
     *
     *  @param[out]  result     will be the joined data
     *  @param[in]   firstRow   first local row
     *  @param[in]   nRows      number of rows to join
     */
    void joinColumnData( hmemo::HArray<ValueType>& result, const IndexType firstRow, const IndexType nRows ) const;

    /** binaryOp but with dense matrices */

    void binaryOpDense( const DenseMatrix<ValueType>& matrixA, const common::BinaryOp op, const DenseMatrix<ValueType>& matrixB );

    /***************************************************************************
     *  Static Methods for dense storage                                        *
     ***************************************************************************/

    /** Split dense storage data according to a distribution into chunks.
     *
     *  @param[out]  chunks        vector of shared pointer to the new allocated chunks
     *  @param[in]   columnData    is the dense storage to split
     *  @param[in]   columnDist    is the distribution used for splitting global data to local data
     *
     *  Note:  columnData.getNumColumns() == columnDist.getGlobalSize()
     *  After: chunks.size() == columnDist.getNumPartitions()
     */
    static void splitColumnData(
        std::vector<std::shared_ptr<DenseStorage<ValueType> > >& chunks,
        const DenseStorage<ValueType>& columnData,
        const dmemo::Distribution& columnDist );

    /** Restrict dense storage of a replicated matrix to its local part according to row distribution.
     *
     * @param[out] local is the local part of the dense storage for this partition
     * @param[in] global is the replicated dense storage
     * @param[in] rowDistribution distribution used to localize the dense storage.
     */

    static void localize(
        DenseStorage<ValueType>& local,
        const DenseStorage<ValueType>& global,
        const dmemo::Distribution& rowDistribution );

    void redistributeRows( dmemo::DistributionPtr rowDistribution );

    /** Split the replicated columns into chunks according to the column distribution. */

    void splitColumns( dmemo::DistributionPtr colDistribution );

    mutable hmemo::HArray<ValueType> mSendValues;
    mutable hmemo::HArray<ValueType> mReceiveValues;

    void computeOwners();

    /** Special implementation of invert in place for a cyclic distributed matrix. */

    void invertReplicated();

    static std::string initTypeName();

public:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    // static methods, variables to register create routine in _Matrix factory of base class.

    static _Matrix* create();

    // key for factory

    static MatrixCreateKeyType createValue();

    MatrixCreateKeyType getCreateValue() const;
};

} /* end namespace lama */

} /* end namespace scai */
