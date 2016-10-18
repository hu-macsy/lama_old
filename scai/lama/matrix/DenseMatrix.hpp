/**
 * @file DenseMatrix.hpp
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
 * @brief Definition of matrix class for distributed matrixes in Dense format.
 * @author Michael Drost
 * @date 22.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/CRTPMatrix.hpp>

// local library
#include <scai/lama/matrix/SparseMatrix.hpp>

#include <scai/lama/storage/DenseStorage.hpp>

// internal scai libraries
#include <scai/common/shared_ptr.hpp>

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

    public CRTPMatrix<DenseMatrix<ValueType>, ValueType>,
    public Matrix::Register<DenseMatrix<ValueType> >    // register at factory
{

public:

    typedef ValueType MatrixValueType; //!< This is the type of the matrix values.

    typedef DenseStorage<ValueType> StorageType;

    typedef common::shared_ptr<DenseStorage<ValueType> > DenseStoragePtr;

    /** Getter for the type name of the class. */

    static const char* typeName();

    /** Default constructor. */

    DenseMatrix();

    /** Constructor of a replicated dense matrix.
     *
     * @param[in] numRows      number of rows, must be non-negative.
     * @param[in] numColumns   number of columns, must be non-negative.
     */
    DenseMatrix( const IndexType numRows, const IndexType numColumns );

    /**
     * Constructor of a distributed dense matrix.
     *
     * @param[in] rowDist   size and distribution of rows
     * @param[in] colDist   size and distribution of columns
     *
     * For consistency with the constructors of sparse matrices the values 
     * of the dense matrix are initialized with 0 here.
     */
    DenseMatrix( dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Overwrites default copy constructor so it uses other copy constructor.
     *
     *  Note: Default copy constructor would not make deep copies of the
     *        storages, so it must be overridden
     */
    DenseMatrix( const DenseMatrix<ValueType>& other );

    /** Constructs a dense matrix from any other matrix that can be of a different type.
     *
     *  @param[in] other   input matrix.
     *  @param[in] transposeFlag if true the input matrix will be transposed
     */
    DenseMatrix( const Matrix& other, bool transposeFlag = false );

    /** Constructor of a (replicated) dense matrix by global storage.
     *
     *  @param[in] globalData  contains the matrix storage
     */
    explicit DenseMatrix( const _MatrixStorage& globalData );

    /** Constructs a dense matrix from any other matrix with new distributions.
     *
     *  @param[in] other             input matrix.
     *  @param[in] rowDistribution   new distribution of rows among processors
     *  @param[in] colDistribution   new distribution of columns for blocking
     *
     *  The following codes are equivalent:
     *
     *  \code
     *      DenseMatrix dense( other, rowDist, colDist );
     *      // same as
     *      DenseMatrix dense( other );
     *      dense->redistribute( rowDist, colDist );
     *  \endcode
     *
     *  The constructor with distributions is more convenient and might be more efficient
     *  due to less memory allocations as less temporary data is needed.
     */
    DenseMatrix( const Matrix& other, dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /** Constructs a dense matrix from another dense matrix with new distributions.
     *
     *  @param[in] other             input matrix.
     *  @param[in] rowDistribution   new distribution of rows among processors
     *  @param[in] colDistribution   new distribution of columns for blocking
     *
     */
    DenseMatrix(
        const DenseMatrix<ValueType>& other,
        dmemo::DistributionPtr rowDistribution,
        dmemo::DistributionPtr colDistribution );

    /** Constructor of a dense matrix by local storage.
     *
     *  @param[in] localData   contains local rows of the distributed matrix
     *  @param[in] rowDist     is distribution of localData
     *  @param[in] colDist     specifies how to split columns of local rows
     *
     *  This constructor works also fine if localData is the full global matrix;
     *  in this case only local rows will be taken on this processor.
     */
    DenseMatrix( const _MatrixStorage& localData, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /**
     * Constructor of a replicated dense matrix from the passed csr sparse matrix.
     *
     * @param[in] numRows       the number of rows of the matrix
     * @param[in] numColumns    the number of columns of the matrix
     * @param[in] numNoneZeros  the number of none zeros of the matrix
     * @param[in] ia            row pointer of the input csr sparse matrix
     * @param[in] ja            column indexes of the input csr sparse matrix
     * @param[in] values        the none zero values of the input csr sparse matrix
     */
    template<typename OtherValueType>
    DenseMatrix(
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numNoneZeros,
        const IndexType* const ia,
        const IndexType* const ja,
        const OtherValueType* const values );

    /**
     * Contructor of a dense matrix by matrix expression alpha * A * B + beta * C
     *
     * @param[in] expression  matrix expression alpha * A * B + beta * C
     */
    DenseMatrix( const Expression_SMM_SM& expression );

    /**
     * Constructor of a dense matrix by matrix expression alpha * A * B
     *
     * @param[in] expression   matrix espression alpha * A * B
     */
    DenseMatrix( const Expression_SMM& expression );

    /**
     * Constructor of a dense matrix by matrix expression alpha * A + beta * b
     *
     * @param[in] expression   matrix espression scalar * matrix + scalar * matrix
     */
    DenseMatrix( const Expression_SM_SM& expression );

    /**
     * Constructor of a dense matrix by matrix expression alpha * A
     *
     * @param[in] expression   matrix expression alpha * A where alpha is a Scalar and A a matrix
     */
    DenseMatrix( const Expression_SM& expression );

    /** Constructor of a replicated dense matrix by reading the matrix
     *  data from a file.
     *
     *  @param[in] filename   Name of the file with matrix data.

     *  Next releases will also support distributed/parallel I/O. In the
     *  meantime this constructor should be used with a following call of
     *  the redistribute method.
     */
    DenseMatrix( const std::string& filename );

    /**
     * Destructor, releases all allocated resources.
     */
    virtual ~DenseMatrix();

    /** Implementation of abstract method for dense matrices. */

    virtual bool isConsistent() const;

    /** Make overloaded operator= available before overriding the default one. */

    using Matrix::operator=;

    /** Overrides the default assignment operator to guarantee deep copy. */

    DenseMatrix& operator=( const DenseMatrix& matrix );

    /** Implementation for Matrix::getTypeName() */

    const char* getTypeName() const;

    /**
     * Gives info about the matrix kind (DENSE).
     */
    virtual Matrix::MatrixKind getMatrixKind() const
    {
        return Matrix::DENSE;
    }

    /* Implementation of pure method of class Matrix. */

    virtual void setContextPtr( const hmemo::ContextPtr context );

    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::setContextPtr; // setContextPtr( localContext, haloContext )

    /* Implementation of pure method of class Matrix. */

    virtual hmemo::ContextPtr getContextPtr() const
    {
        return mData[0]->getContextPtr();
    }

    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::setIdentity; // setIdentity( const IndexType n )

    /** Implementation of pure method Matrix::setIdentity. */

    virtual void setIdentity( dmemo::DistributionPtr distribution );

    /** Implementation of pure Matrix::setDenseData */

    virtual void setDenseData(
        dmemo::DistributionPtr rowDistribution,
        dmemo::DistributionPtr colDistribution,
        const hmemo::_HArray& values,
        const Scalar eps );

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
        const hmemo::HArray<IndexType>& offset,
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

    /* Implementation of pure method of class Matrix. */

    virtual void clear();

    /* Implementation of pure method Matrix::purge. */

    virtual void purge();

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( const IndexType numRows, const IndexType numColumns );

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const Matrix& other );

    /* Implementation of pure method of class Matrix. */

    virtual void assignTranspose( const Matrix& other );

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

    /** Method that assigns a sparse matrix, specialization of assign( const Matrix& ) */

    void assignSparse( const CRTPMatrix<SparseMatrix<ValueType>, ValueType>& other );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const _MatrixStorage& storage );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const _MatrixStorage& storage, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** @brief TODO[doxy] Complete Description.
     *
     * @param[in] other   TODO[doxy] Complete Description.
     */
    void assignLocal( const _MatrixStorage& other );

    /** Implementation of Matrix::buildLocalStorage. */

    virtual void buildLocalStorage( _MatrixStorage& storage ) const;

    /* Implementation of pure method of class Matrix. */

    void redistribute( dmemo::DistributionPtr rowDistribution, dmemo::DistributionPtr colDistribution );

    /* Implementation of pure method of class Matrix. */

    virtual void getDiagonal( Vector& diagonal ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Vector& diagonal );

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Scalar diagonalValue );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Vector& values );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Scalar value );

    /* Implementation of pure method of class Matrix. */

    virtual void conj();

    /* Implementation of pure method of class Matrix. */

    virtual Scalar getValue( IndexType i, IndexType j ) const;

    /* Implemenation of pure method of class Matrix */

    virtual void matrixTimesScalar( const Matrix& other, const Scalar alpha );

    /**
     *  @brief Matrix times vector with same value types and correct distributions.
     *
     * @param[out] denseResult   TODO[doxy] Complete Description.
     * @param[in]  alphaValue    TODO[doxy] Complete Description.
     * @param[in]  denseX        TODO[doxy] Complete Description.
     * @param[in]  betaValue     TODO[doxy] Complete Description.
     * @param[in]  denseY        TODO[doxy] Complete Description.
     *
     *  Note: Matrix::matrixTimesMatrix is implemented in the CRTPMatrix class.
     *        that requires this method.
     *
     *  Note: all vectors must have the right distribution.
     */
    void matrixTimesVectorImpl(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>& denseY ) const;

    void vectorTimesMatrixImpl(
        DenseVector<ValueType>& denseResult,
        const ValueType alphaValue,
        const DenseVector<ValueType>& denseX,
        const ValueType betaValue,
        const DenseVector<ValueType>& denseY ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void invert( const Matrix& other );

    /** Invert in place */

    void invert()
    {
        this->invert( *this );
    }

    /* Implementation of pure method of class Matrix. */
    virtual Scalar l1Norm() const;

    /* Implementation of pure method of class Matrix. */
    virtual Scalar l2Norm() const;

    /** Implementation of Matrix::maxNorm for dense matrices. */

    virtual Scalar maxNorm() const;

    /** Implementation of Matrix::maxDiffNorm for dense matrices. */

    virtual Scalar maxDiffNorm( const Matrix& other ) const;

    /** Get the maximal difference between two elements for dense matrices of same type. */

    ValueType maxDiffNormImpl( const DenseMatrix<ValueType>& other ) const;

    /* Implemenation of pure method of class Matrix */

    virtual void matrixPlusMatrix( const Scalar alpha, const Matrix& A, const Scalar beta, const Matrix& B );

    /** Implementation of pure method Matrix::matrixTimesMatrix */

    void matrixTimesMatrix(
        Matrix& result,
        const Scalar alpha,
        const Matrix& x,
        const Scalar beta,
        const Matrix& y ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void prefetch() const;

    /** @brief TODO[doxy] Complete Description.
     *
     * @param[in] loc   TODO[doxy] Complete Description.
     */
    void prefetch( hmemo::ContextPtr loc ) const;

    /* Implementation of pure method of class Matrix. */

    void wait() const;

    /** This method returns the storage containing the local data regarding row/col distribution. */

    const DenseStorage<ValueType>& getLocalStorage() const;

    DenseStorage<ValueType>& getLocalStorage();

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getLocalNumValues() const;

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getLocalNumRows() const;

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getLocalNumColumns() const;

    /* Implementation of pure method of class Matrix. */

    virtual IndexType getNumValues() const;

    /* Implementation of pure method of class Matrix. */

    virtual bool hasDiagonalProperty() const;

    /* Implementation of pure method of class Matrix. */

    virtual void resetDiagonalProperty();

    /* Implementation of method writeAt for dense matrix. */

    virtual void writeAt( std::ostream& stream ) const;

    /* Implementation of pure method of class Matrix. */

    virtual common::scalar::ScalarType getValueType() const;

    virtual size_t getValueTypeSize() const;

    /**
     * @brief Implementation of pure function Matrix::copy with covariant return type.
     */
    virtual DenseMatrix<ValueType>* newMatrix() const;

    /**
     * @brief Implementation of pure function Matrix::copy with covariant return type.
     */
    virtual DenseMatrix<ValueType>* copy() const;

    /* Implementation of pure method Matrix::getFormat */

    virtual Format::MatrixStorageFormat getFormat() const
    {
        return Format::DENSE;
    }

    /* Implementation of pure method of class Matrix. */

    virtual size_t getMemoryUsage() const;

    /** local data is allocated in chunks according to column distribution */

    std::vector<common::shared_ptr<DenseStorage<ValueType> > > mData;

    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::getNumRows;
    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::getNumColumns;

    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::getRowDistribution;
    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::getRowDistributionPtr;
    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::getColDistribution;
    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::getColDistributionPtr;

    const utilskernel::LArray<PartitionId>& getOwners() const
    {
        return mOwners;
    }

    /** Get a complete row of the local storage, used by getRow in CRTPMatrix */

    void getLocalRow( hmemo::HArray<ValueType>& row, const IndexType iLocal ) const;

    void setLocalRow( const hmemo::HArray<ValueType>& row, 
                      const IndexType localRowIndex,
                      const utilskernel::binary::BinaryOp op  );

    void getLocalColumn( hmemo::HArray<ValueType>& col, const IndexType colIndex ) const;

    void setLocalColumn( const hmemo::HArray<ValueType>& column, 
                         const IndexType colIndex,
                         const utilskernel::binary::BinaryOp op  );

    /** Copy a dense matrix with different data type; inherits sizes and distributions */

    template<typename otherT>
    void copyDenseMatrix( const DenseMatrix<otherT>& other );

    /** Optimized implementation for dense vectors as diagonal. */

    template<typename OtherT>
    void getDiagonalImpl( DenseVector<OtherT>& diagonal ) const;

protected:

    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::mNumRows;
    using CRTPMatrix<DenseMatrix<ValueType>, ValueType>::mNumColumns;

    utilskernel::LArray<PartitionId> mOwners;

    /**
     * @brief Set this matrix = alpha * A + beta * B
     *
     * @param[in]  alpha    TODO[doxy] Complete Description.
     * @param[in]  A        TODO[doxy] Complete Description.
     * @param[in]  beta     TODO[doxy] Complete Description.
     * @param[in]  B        TODO[doxy] Complete Description.
     */
    void matrixPlusMatrixImpl(
        const ValueType alpha,
        const DenseMatrix<ValueType>& A,
        const ValueType beta,
        const DenseMatrix<ValueType>& B );

private:

    /** Allocate of storage for the column blocks. */

    void allocateData();

    /** Join column data of column distributed dense data
     *
     *  @param[out]  result     will be the joined data
     *  @param[in]   firstRow   first local row
     *  @param[in]   nRows      number of rows to join
     *
     *  Note: the column distribution is taken from the values of mOwners
     */
    void joinColumnData( scai::hmemo::HArray<ValueType>& result, const IndexType firstRow, const IndexType nRows ) const;

    /***************************************************************************
     *  Static Methods for dense storage                                        *
     ***************************************************************************/

    /** Split dense storage data according to a distribution into chunks.
     *
     *  @param[out]  chunks        vector of shared pointer to the new allocated chunks
     *  @param[in]   columnData    is the dense storage to split
     *  @param[in]   numChunks     is the number of chunks, same as number of partitions of distribution
     *  @param[in]   columnOwners  owner for each column
     */
    static void splitColumnData(
        std::vector<common::shared_ptr<DenseStorage<ValueType> > >& chunks,
        const DenseStorage<ValueType>& columnData,
        const PartitionId numChunks,
        const hmemo::HArray<PartitionId>& columnOwners );

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

    //TODO: no implementation: implement or delete
    //void initChunks();  // common initialization for constructors

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    void    computeOwners();

    /** Special implementation of invert in place for a cyclic distributed matrix. */

    void invertReplicated();

    static std::string initTypeName();

public:

    // static methods, variables to register create routine in Matrix factory of base class.

    static Matrix* create();

    // key for factory

    static MatrixCreateKeyType createValue();

    MatrixCreateKeyType getCreateValue() const;
};

/*  template methods implementations */

template<typename ValueType>
template<typename OtherValueType>
void DenseMatrix<ValueType>::copyDenseMatrix( const DenseMatrix<OtherValueType>& other )
{
    // check for valid pointer, might be dynamic cast went wrong somewhere else
    //SCAI_ASSERT_ERROR( &other, "NULL matrix in assignment operator" )
    SCAI_LOG_INFO( logger, "copy dense, this = " << this << ", other = " << &other )
    // inherit size and distributions
    Matrix::setDistributedMatrix( other.getRowDistributionPtr(), other.getColDistributionPtr() );
    mData.resize( other.mData.size() );
    IndexType n = static_cast<IndexType>( other.mData.size() );

    for ( IndexType i = 0; i < n; ++i )
    {
        SCAI_LOG_DEBUG( logger, "copy block " << i << " of " << n << " = " << *other.mData[i] )
        mData[i].reset( new DenseStorage<ValueType>( *other.mData[i] ) );
    }

    mOwners = other.getOwners();
}

template<typename ValueType>
template<typename OtherValueType>
DenseMatrix<ValueType>::DenseMatrix(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numNoneZeros,
    const IndexType* const ia,
    const IndexType* const ja,
    const OtherValueType* const values )

    : CRTPMatrix<DenseMatrix<ValueType>, ValueType>( numRows, numColumns )
{
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( mNumRows, mNumColumns ) );
    mData[0]->setCSRData( numNoneZeros, ia, ja, values );
    computeOwners();
}

} /* end namespace lama */

} /* end namespace scai */
