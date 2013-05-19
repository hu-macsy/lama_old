/**
 * @file DenseMatrix.hpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Definition of matrix class for distributed matrixes in Dense format.
 * @author Michael Drost
 * @date 22.02.2011
 * $Id$
 */
#ifndef LAMA_DENSEMATRIX_HPP_
#define LAMA_DENSEMATRIX_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/matrix/CRTPMatrix.hpp>

//others
#include <lama/matrix/SparseMatrix.hpp>

#include <lama/storage/DenseStorage.hpp>

//boost
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

//boost
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

namespace lama
{

template<typename T> class DenseVector;
// forward declaration

/** Class for dense matrices where rows are distributed among rows and columns
 *  are splitted according to a column distribution.
 *
 *  The local rows are splitted according to a column distribution into different
 *  blocks that are stored in a vector of shared pointers. A copy of this vector
 *  will not make deep copies of the blocks so default copy constructor and
 *  assignment operator must be overridden.
 *
 *  @tparam T is the value type of the matrix values.
 */

template<typename T>
class LAMA_DLL_IMPORTEXPORT DenseMatrix: public CRTPMatrix<DenseMatrix<T>,T>
{

public:

    typedef T ValueType;

    typedef boost::shared_ptr<DenseStorage<T> > DenseStoragePtr;

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
     * @param[in] rowDist   TODO[doxy] Complete Description.
     * @param[in] colDist   TODO[doxy] Complete Description.
     */
    DenseMatrix( DistributionPtr rowDist, DistributionPtr colDist );

    /**
     * Constructor of a square unity matrix.
     *
     * @param[in] dist   TODO[doxy] Complete Description.
     */
    explicit DenseMatrix( DistributionPtr dist );

    /** Overwrites default copy constructor so it uses other copy constructor.
     *
     *  Note: Default copy constructor would not make deep copies of the
     *        storages, so it must be overridden
     */
    DenseMatrix( const DenseMatrix<ValueType>& other );

    /** Constructs a dense matrix from any other matrix that can be of a different type.
     *
     *  @param[in] other   input matrix.
     *
     *  New dense matrix has the same size and the same distribution.
     */
    DenseMatrix( const Matrix& other );

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
    DenseMatrix( const Matrix& other, DistributionPtr rowDistribution, DistributionPtr colDistribution );

    /** Constructs a dense matrix from another dense matrix with new distributions.
     *
     *  @param[in] matrix            input matrix.
     *  @param[in] rowDistribution   new distribution of rows among processors
     *  @param[in] colDistribution   new distribution of columns for blocking
     *
     */
    DenseMatrix(
        const DenseMatrix<ValueType>& matrix,
        DistributionPtr rowDistribution,
        DistributionPtr colDistribution );

    /** Constructor of a dense matrix by local storage.
     *
     *  @param[in] localData   contains local rows of the distributed matrix
     *  @param[in] rowDist     is distribution of localData
     *  @param[in] colDist     specifies how to split columns of local rows
     *
     *  This constructor works also fine if localData is the full global matrix;
     *  in this case only local rows will be taken on this processor.
     */
    DenseMatrix( const _MatrixStorage& localData, DistributionPtr rowDist, DistributionPtr colDist );

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
     * Contructor of a dense matrix by matrix expression alhpa * A * B + beta * C
     *
     * @param[in] expression  matrix expression alhpa * A * B + beta * C
     */
    DenseMatrix(
        const Expression<Expression<Scalar, Expression<Matrix, Matrix, Times>, Times>,
                         Expression<Scalar, Matrix, Times>, 
                         Plus> expression );

    /**
     * Contructor of a dense matrix by matrix expression A * B 
     *
     * @param[in] expression   matrix expression A * B
     */
    DenseMatrix( const Expression<Matrix, Matrix, Times> expression );

    /**
     * Constructor of a dense matrix by matrix espression alhpa * A * B
     *
     * @param[in] expression   matrix espression alhpa * A * B
     */
    DenseMatrix( const Expression<Scalar, Expression<Matrix, Matrix, Times>, Times> expression );

    /**
     * Constructor of a dense matrix by matrix expression alhpa * A
     *
     * @param[in] expression   matrix expression alpha * A where alpha is a Scalar and A a matrix
     */
    DenseMatrix( const Expression<Scalar, Matrix, Times> expression );

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

    /** Overrides the default assignment operator to guarantee deep copy. */

    DenseMatrix& operator=( const DenseMatrix& matrix );

    /** Redefines assignment operator to get the correct return value; implementation is same as for base classes. */

    DenseMatrix& operator=( const Matrix& matrix );

    DenseMatrix& operator=( const Expression<Scalar, Matrix, Times> expression );

    DenseMatrix& operator=( const Expression<Matrix, Matrix, Times> expression );

    DenseMatrix& operator=( const Expression<Scalar, Expression<Matrix, Matrix, Times>, Times> expression );

    DenseMatrix& operator=( const Expression<Expression<Scalar, Expression<Matrix, Matrix, Times>, Times>,
                                             Expression<Scalar, Matrix, Times>,
                                             Plus> expression );

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

    virtual void setContext( const ContextPtr context );

    using CRTPMatrix<DenseMatrix<T>, T>::setContext; // setContext( localContext, haloContext )

    /* Implementation of pure method of class Matrix. */

    virtual ContextPtr getContextPtr() const
    {
        return mData[0]->getContextPtr();
    }

    using CRTPMatrix<DenseMatrix<T>, T>::setIdentity;  // setIdentity( const IndexType n )

    /** Implementation of pure method Matrix::setIdentity. */

    virtual void setIdentity( DistributionPtr distribution );

    /** Implementation of pure Matrix::setDenseData */

    virtual void setDenseData( 
        DistributionPtr rowDistribution,
        DistributionPtr colDistribution,
        const _LAMAArray& values, 
        const double eps );

    /** Implementation for pure method Matrix::setCSRData. */

    virtual void setCSRData(
        DistributionPtr rowDist,
        DistributionPtr colDist,
        const IndexType numValues,
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const _LAMAArray& values );

    /** Implementation of pure method for the dense storage format. */

    virtual void buildCSRData( LAMAArray<IndexType>& rowIA, LAMAArray<IndexType>& rowJA, _LAMAArray& rowValues ) const;

    /** Implementation of pure method. */

    virtual void setCSRData(
        const LAMAArray<IndexType>& rowIA,
        const LAMAArray<IndexType>& rowJA,
        const _LAMAArray& rowValues,
        DistributionPtr rowDistribution,
        DistributionPtr colDistribution );

    /** Local version of setCSRData . */

    void setCSRDataLocal(
        const LAMAArray<IndexType>& rowIA,
        const LAMAArray<IndexType>& rowJA,
        const _LAMAArray& rowValues ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void clear();

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( const IndexType numRows, const IndexType numColumns );

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const Matrix& other );

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

    void assignSparse( const _SparseMatrix& other );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const _MatrixStorage& storage );

    /* Implementation of pure method of class Matrix. */

    virtual void assign( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist );

    /** @brief TODO[doxy] Complete Description.
     *
     * @param[in] other   TODO[doxy] Complete Description.
     */
    void assignLocal( const _MatrixStorage& other );

    /** Implementation of Matrix::buildLocalStorage. */

    virtual void buildLocalStorage( _MatrixStorage& storage ) const;

    /* Implementation of pure method of class Matrix. */

    void redistribute( DistributionPtr rowDistribution, DistributionPtr colDistribution );

    /* Implementation of pure method of class Matrix. */

    virtual void getDiagonal( Vector& diagonal ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void getRow( Vector& row, const IndexType globalRowIndex ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Vector& diagonal );

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Scalar diagonalValue );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Vector& values );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Scalar value );

    /* Implementation of pure method of class Matrix. */

    virtual Scalar getValue( lama::IndexType i, lama::IndexType j ) const;

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

    /* Implementation of pure method of class Matrix. */

    virtual void invert( const Matrix& other );

    /** Invert in place */

    void invert()
    {
        this->invert( *this );
    }

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
    void prefetch( lama::ContextPtr loc ) const;

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

    std::vector<DenseStoragePtr>& getCyclicLocalValues();

    const std::vector<DenseStoragePtr>& getCyclicLocalValues() const;

    /* Implementation of pure method of class Matrix. */

    virtual bool hasDiagonalProperty() const;

    /* Implementation of pure method of class Matrix. */

    virtual void resetDiagonalProperty();

    /* Implementation of method writeAt for dense matrix. */

    virtual void writeAt( std::ostream& stream ) const;

    /* Implementation of pure method of class Matrix. */

    virtual Scalar::ScalarType getValueType() const;

    /** Method writes dense matrix to a file.
     *
     *  Writing is only supported for a replicated matrix.
     */

    void writeToFile(
        const std::string& fileName,
        const File::FileType fileType = File::BINARY,
        const File::DataType dataType = File::INTERNAL,
        const File::IndexDataType indexDataTypeIA = File::INT,
        const File::IndexDataType indexDataTypeJA = File::INT ) const;
    /**
     * @brief Implementation of pure function Matrix::create with covariant return type. 
     */
    virtual DenseMatrix<ValueType>* create() const;

    /**
     * @brief Implementation of pure function Matrix::copy with covariant return type. 
     */
    virtual DenseMatrix<ValueType>* copy() const;

    /* Implementation of pure method of class Matrix. */

    virtual size_t getMemoryUsage() const;

    /** local data is allocated in chunks according to column distribution */

    std::vector<boost::shared_ptr<DenseStorage<ValueType> > > mData;

    using CRTPMatrix<DenseMatrix<T>,T>::getNumRows;
    using CRTPMatrix<DenseMatrix<T>,T>::getNumColumns;

    using CRTPMatrix<DenseMatrix<T>,T>::getDistribution;
    using CRTPMatrix<DenseMatrix<T>,T>::getDistributionPtr;
    using CRTPMatrix<DenseMatrix<T>,T>::getColDistribution;
    using CRTPMatrix<DenseMatrix<T>,T>::getColDistributionPtr;

    const std::vector<PartitionId>& getOwners() const
    {
        return mOwners;
    }

protected:

    using CRTPMatrix<DenseMatrix<T>,T>::mNumRows;
    using CRTPMatrix<DenseMatrix<T>,T>::mNumColumns;

    std::vector<PartitionId> mOwners;

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

    /***************************************************************************
     *  Static Methods for dense storage                                        *
     ***************************************************************************/

    /** Join dense storage of column distributed data
     *
     *  @param[out]  result        will be the joined storage
     *  @param[in]   chunks        is a vector of all chunks, one chunk for each partition
     *  @param[in]   columnOwners  vector with owner for each column
     *
     *  Note: the column distribution itself is given implicitly by the vector of owners
     */
    static void joinColumnData(
        DenseStorage<ValueType>& result,
        const std::vector<boost::shared_ptr<DenseStorage<ValueType> > >& chunks,
        const std::vector<IndexType>& columnOwners );

    /** Split dense storage data according to a distribution into chunks.
     *
     *  @param[out]  chunks        vector of shared pointer to the new allocated chunks
     *  @param[in]   columnData    is the dense storage to split
     *  @param[in]   numChunks     is the number of chunks, same as number of partitions of distribution
     *  @param[in]   columnOwners  owner for each column
     */
    static void splitColumnData(
        std::vector<boost::shared_ptr<DenseStorage<ValueType> > >& chunks,
        const DenseStorage<ValueType>& columnData,
        const PartitionId numChunks,
        const std::vector<IndexType>& columnOwners );

    /** Restrict dense storage of a replicated matrix to its local part according to row distribution.
     *
     * @param[out] local is the local part of the dense storage for this partition
     * @param[in] global is the replicated dense storage
     * @param[in] rowDistribution distribution used to localize the dense storage.
     */

    static void localize(
        DenseStorage<ValueType>& local,
        const DenseStorage<ValueType>& global,
        const Distribution& rowDistribution );

    /** Copy a dense matrix with different data type; inherits sizes and distributions */

    template<typename otherT>
    void copyDenseMatrix( const DenseMatrix<otherT>& other );

    /** Optimized implementation for dense vectors as diagonal. */

    template<typename OtherT>
    void getDiagonalImpl( DenseVector<OtherT>& diagonal ) const;

    void redistributeRows( DistributionPtr rowDistribution );

    /** Split the replicated columns into chunks according to the column distribution. */

    void splitColumns( DistributionPtr colDistribution );

    void getRow( DenseVector<ValueType>& row, const IndexType i ) const;

    mutable LAMAArray<ValueType> mSendValues;
    mutable LAMAArray<ValueType> mReceiveValues;

    //TODO: no implementation: implement or delete
    //void initChunks();  // common initialization for constructors

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    void computeOwners();

    /** @brief Predicate to check if SCALapack is supported via LAMAInterface. */

    bool hasScalaPack();

    /** Special implementation of invert in place for a cyclic distributed matrix. */

    void invertCyclic();

    void invertReplicated();
};

/*  template methods implementations */

template<typename ValueType>
template<typename OtherValueType>
void DenseMatrix<ValueType>::copyDenseMatrix( const DenseMatrix<OtherValueType>& other )
{
    // check for valid pointer, might be dynamic cast went wrong somewhere else

    LAMA_ASSERT_ERROR( &other, "NULL matrix in assignment operator" )

    LAMA_LOG_INFO( logger, "copy dense, this = " << this << ", other = " << &other )

    // inherit size and distributions

    Matrix::setDistributedMatrix( other.getDistributionPtr(), other.getColDistributionPtr() );

    mData.resize( other.mData.size() );

    IndexType n = static_cast<IndexType>( other.mData.size() );

    for ( IndexType i = 0; i < n; ++i )
    {
        LAMA_LOG_DEBUG( logger, "copy block " << i << " of " << n << " = " << *other.mData[i] )

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

    : CRTPMatrix<DenseMatrix<ValueType>,ValueType>( numRows, numColumns )
{
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( mNumRows, mNumColumns ) );
    mData[0]->setCSRData( numNoneZeros, ia, ja, values );
    computeOwners();
}

} // namespace lama

#endif // LAMA_DENSEMATRIX_HPP_
