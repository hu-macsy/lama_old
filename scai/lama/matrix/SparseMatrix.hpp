/**
 * @file SparseMatrix.hpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Definition of class for distributed sparse matrices.
 * @author Jiri Kraus, Thomas Brandes
 * @date 06.06.2011
 * @since 1.0.0
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

using tasking::SyncToken;

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
class COMMON_DLL_IMPORTEXPORT SparseMatrix: public CRTPMatrix<SparseMatrix<ValueType>,ValueType>
{

    friend class SpecializedJacobi;

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
    ;

    /** Getter routine for halo part of the sparse matrix. */

    const MatrixStorage<ValueType>& getHaloStorage() const
    {
        return *mHaloData;
    }
    ;

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

    SparseMatrix( common::shared_ptr<MatrixStorage<ValueType> > storage, DistributionPtr rowDist );

    SparseMatrix(
        common::shared_ptr<MatrixStorage<ValueType> > storage,
        DistributionPtr rowDist,
        DistributionPtr colDist );

    /** Constructor of a sparse matrix with local and halo data available. */

    SparseMatrix(
        common::shared_ptr<MatrixStorage<ValueType> > localData,
        common::shared_ptr<MatrixStorage<ValueType> > haloData,
        const Halo& halo,
        DistributionPtr rowDist,
        DistributionPtr colDist );

    SparseMatrix( const Matrix& matrix, const bool transposeFlag = false );

    SparseMatrix( const Matrix& other, DistributionPtr rowDist, DistributionPtr colDist );

    /** Override also the default copy constructor that does not make a
     *  deep copy of the input matrix due to the use of shared pointers.
     */
    SparseMatrix( const SparseMatrix<ValueType>& other );

    /** Implementation of abstract method for sparse matrices. */

    virtual bool isConsistent() const;

    /* Implementation of pure method of class Matrix. */

    virtual void invert( const Matrix& other );

    /* Implementation of pure method of class Matrix. */

    virtual void setContext( const ContextPtr context )
    {
        setContext( context, context );
    }

    /* Implementation of pure method of class Matrix. */

    virtual void setContext( const ContextPtr localContext, const ContextPtr haloContext )
    {
        mLocalData->setContext( localContext );
        mHaloData->setContext( haloContext );
    }

    /* Implementation of pure method of class Matrix. */

    virtual ContextPtr getContextPtr() const
    {
        return mLocalData->getContextPtr();
    }

    /** Implementation for Matrix::setDenseData */

    virtual void setDenseData(
        DistributionPtr rowDistribution,
        DistributionPtr colDistribution,
        const ContextArray& values,
        Scalar eps = Scalar( 0 ) );

    /** Implementation for pure method Matrix::setCSRData. */

    virtual void setCSRData(
        DistributionPtr rowDist,
        DistributionPtr colDist,
        const IndexType numValues,
        const LAMAArray<IndexType>& ia,
        const LAMAArray<IndexType>& ja,
        const ContextArray& values );

    /* Implementation of pure method of class Matrix. */

    virtual void clear();

    /* Implementation of pure method of class Matrix. */

    virtual void purge();

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( const IndexType numRows, const IndexType numColumns );

    /* Implementation of pure method of class Matrix. */

    virtual void allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution );

    /* Before overriding the virtual function make the other routine setIdentity( int n ) visible */

    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::setIdentity;

    /** Set matrix to a identity square matrix with same row and column distribution. */

    virtual void setIdentity( DistributionPtr distribution );

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

    virtual void assign( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist );

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

    virtual void getRow( Vector& row, const IndexType globalRowIndex ) const;

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Vector& diagonal );

    /* Implementation of pure method of class Matrix. */

    virtual void setDiagonal( const Scalar scalar );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Vector& scaling );

    /* Implementation of pure method of class Matrix. */

    virtual void scale( const Scalar scaling );

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

    /* Implemenation of pure method of class Matrix */

    virtual void matrixTimesVector(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const;

    void vectorTimesMatrix(
        Vector& result,
        const Scalar alpha,
        const Vector& x,
        const Scalar beta,
        const Vector& y ) const;

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
        LAMAArray<ValueType>& localResult,
        const LAMAArray<ValueType>& localX,
        LAMAArray<ValueType>& haloX,
        common::function<
        void(
            const MatrixStorage<ValueType>* localMatrix,
            LAMAArray<ValueType>& localResult,
            const LAMAArray<ValueType>& localX )> localF,
        common::function<
        void(
            const MatrixStorage<ValueType>* haloMatrix,
            LAMAArray<ValueType>& localResult,
            const LAMAArray<ValueType>& haloX )> haloF ) const;

    /**
     * @brief Operation on distributed matrix with halo exchange, async version
     *
     * The asynchronous version starts the local computation asynchronously so it
     * can overlap with the halo exchange.
     */
    void haloOperationAsync(
        LAMAArray<ValueType>& localResult,
        const LAMAArray<ValueType>& localX,
        LAMAArray<ValueType>& haloX,
        common::function<
        SyncToken*(
            const MatrixStorage<ValueType>* localMatrix,
            LAMAArray<ValueType>& localResult,
            const LAMAArray<ValueType>& localX )> localAsyncF,
        common::function<
        void(
            const MatrixStorage<ValueType>* haloMatrix,
            LAMAArray<ValueType>& localResult,
            const LAMAArray<ValueType>& haloX )> haloF ) const;

    void vectorHaloOperationSync(
        LAMAArray<ValueType>& localResult,
        const LAMAArray<ValueType>& localX,
        const LAMAArray<ValueType>& localY,
        common::function<
        void(
            const MatrixStorage<ValueType>* localMatrix,
            LAMAArray<ValueType>& localResult,
            const LAMAArray<ValueType>& localX )> calcF,
        common::function<
        void(
            LAMAArray<ValueType>& localResult,
            const LAMAArray<ValueType>& localX,
            const LAMAArray<ValueType>& localY )> addF ) const;

    void vectorHaloOperationAsync(
        LAMAArray<ValueType>& localResult,
        const LAMAArray<ValueType>& localX,
        const LAMAArray<ValueType>& localY,
        common::function<
        SyncToken*(
            const MatrixStorage<ValueType>* localMatrix,
            LAMAArray<ValueType>& localResult,
            const LAMAArray<ValueType>& localX )> calcF,
        common::function<
        /*SyncToken**/void(
            LAMAArray<ValueType>& localResult,
            const LAMAArray<ValueType>& localX,
            const LAMAArray<ValueType>& localY )> addF ) const;

    /* Implemenation of pure method of class Matrix */

    virtual void matrixPlusMatrix( const Scalar alpha, const Matrix& A, const Scalar beta, const Matrix& B );

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

    /**
     * @brief Read access to the halo of the distributed matrix.
     *
     * @return   reference to the halo of the distributed matrix
     */
    const Halo& getHalo() const;

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

    /* Implementation of pure method Matrix::create with covariant return type */

    virtual SparseMatrix<ValueType>* clone() const;

    /* Implementation of pure method Matrix::copy with covariant return type */

    virtual SparseMatrix<ValueType>* copy() const;

    /* Implementation of pure method of class Matrix. */

    void redistribute( DistributionPtr rowDistribution, DistributionPtr colDistribution );

    /**  */
    /**
     * @brief Assign another matrix transposed to this matrix.
     *
     * @param[in] matrix   input matrix that will be transposed
     */
    void assignTranspose( const Matrix& matrix );

    /* Implementation of pure method of class Matrix. */

    virtual size_t getMemoryUsage() const;

    /** Writes this sparse matrix to a file in CSR format. */

    void writeToFile(
        const std::string& fileName,
        const File::FileType fileType = File::BINARY,
        const File::DataType dataType = File::INTERNAL,
        const File::IndexDataType indexDataTypeIA = File::INT,
        const File::IndexDataType indexDataTypeJA = File::INT ) const;

    /**
     * @brief Assigns this matrix with a replicated sparse matrix read from file.
     *
     * Creates a replicated sparse matrix read from file. Currently supported is
     * the matrix market format, XDR, formatted, unformatted, binary.
     *
     * TODO: set reference to description in StorageIO.
     *
     * @param[in] filename      the filename to read from
     *
     * Note: Derived classes might use this routine within a constructor for convenience.
     *       This class does not support such a constructor as no file format is known.
     */
    void readFromFile( const std::string& filename );

    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::operator=; // make overloaded routines visible before overwriting one

    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::getColDistribution;
    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::getColDistributionPtr;
    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::getDistribution;
    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::getDistributionPtr;
    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::setDistributionPtr;

    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::getCommunicationKind;

    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::getNumColumns;
    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::getNumRows;

    /** Override the default assignment operator to guarantee deep copy. */

//    SparseMatrix& operator=( const SparseMatrix& matrix );
    SparseMatrix<ValueType>& operator=( const SparseMatrix<ValueType>& matrix );

    /**
     * Gives info about the matrix kind (SPARSE).
     */
    virtual Matrix::MatrixKind getMatrixKind() const
    {
        return Matrix::SPARSE;
    }

protected:

    /** Test consistency of sparse matrix data, only used if ASSERT_DEBUG is enabled. */

    void checkSettings();

    common::shared_ptr<MatrixStorage<ValueType> > mLocalData; //!< local columns of sparse matrix

    common::shared_ptr<MatrixStorage<ValueType> > mHaloData; //!< local columns of sparse matrix

    Halo mHalo; //!< Exchange plans for halo part due to column distribution

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

    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::mNumRows;
    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::mNumColumns;
    using CRTPMatrix<SparseMatrix<ValueType>,ValueType>::mColDistribution;
private:

    /**
     * @brief Default constructor is disabled.
     */
    SparseMatrix();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

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
    void    set( const MatrixStorage<ValueType>& otherLocalData, DistributionPtr otherDist );

    /** Implementation of transposed assign for sparse matrix of a known value type. */

    void assignTransposeImpl ( const SparseMatrix<ValueType>& matrix );

    /** Get a complete row of local part only. */

    void getLocalRow( DenseVector<ValueType>& row, const IndexType iLocal ) const;

    mutable LAMAArray<ValueType> mTempSendValues; //!< temporary vector for halo communications
};

} /* end namespace lama */

} /* end namespace scai */
