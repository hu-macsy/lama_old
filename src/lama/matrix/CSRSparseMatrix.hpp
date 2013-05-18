/**
 * @file CSRSparseMatrix.hpp
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
 * @brief Definition of matrix class for distributed sparse matrixes in CSR format.
 * @author Jiri Kraus, Thomas Brandes
 * @date 22.02.2011
 * $Id$
 */
#ifndef LAMA_CSR_SPARSE_MATRIX_HPP_
#define LAMA_CSR_SPARSE_MATRIX_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/matrix/SparseMatrix.hpp>

// others
#include <lama/storage/CSRStorage.hpp>

#include <lama/distribution/GeneralDistribution.hpp>
#include <lama/distribution/NoDistribution.hpp>

namespace lama
{

/** Definition of a derived class for SparseMatrix that uses the CSR storage
 *  format for the local and halo data of the distributed sparse matrix.
 *
 *  As the storage format is known here this class can offer more advanced
 *  constructors that do not exist for SparseMatrix as there the storage
 *  format is not fixed.
 */

template<typename T>
class LAMA_DLL_IMPORTEXPORT CSRSparseMatrix: public SparseMatrix<T>
{

public:

    /** type definition of the value type. */

    typedef T ValueType;

    typedef CSRStorage<T> StorageType;

    /** Static method that returns the name of the matrix class. */

    static const char* typeName();

    /** Default constructor, creates a replicated matrix of size 0 x 0 */

    CSRSparseMatrix();

    /** Constructor, creates a replicated zero-matrix of size numRows x numColums */

    CSRSparseMatrix( const IndexType numRows, const IndexType numColumns );

    /** Constructor, creates a distributed zero-matrix by given row and column distribution */

    CSRSparseMatrix( DistributionPtr rowDist, DistributionPtr colDist )

        : SparseMatrix<ValueType>( createStorage( rowDist->getLocalSize(), colDist->getGlobalSize() ),
                                   rowDist, colDist )
    {
        // Note: splitting of local rows to local + halo part is done by SparseMatrix constructor
    }

    /** Override default constructor, make sure that deep copies are created. */

    CSRSparseMatrix( const CSRSparseMatrix& other )

        : SparseMatrix<ValueType>( createStorage() )

    {
        this->setCommunicationKind( other.getCommunicationKind() );
        this->setContext( other.getContextPtr() );
        SparseMatrix<ValueType>::assign( other );
    }

    CSRSparseMatrix( const Matrix& other, bool transposeFlag = false )

        : SparseMatrix<ValueType>( createStorage() )

    {
        this->setCommunicationKind( other.getCommunicationKind() );

        if ( transposeFlag )
        {
            SparseMatrix<ValueType>::assignTranspose( other );
        }
        else
        {
            SparseMatrix<ValueType>::assign( other );
        }
    }

    /** Constructor of a sparse matrix by another input matrix.
     *
     * @param[in] other     is the input matrix.
     * @param[in] rowDist   TODO[doxy] Complete Description.
     * @param[in] colDist   TODO[doxy] Complete Description.
     */
    CSRSparseMatrix( const Matrix& other, DistributionPtr rowDist, DistributionPtr colDist )

        : SparseMatrix<ValueType>( createStorage() )

    {
        this->setCommunicationKind( other.getCommunicationKind() );

        // ToDo: this could be done better to avoid intermediate copies

        SparseMatrix<ValueType>::assign( other );
        this->redistribute( rowDist, colDist );
    }

    /** Constructor of a (replicated) sparse matrix by global storage.
     *
     *  @param[in] globalData  contains local rows of the distributed matrix
     */
    CSRSparseMatrix( const _MatrixStorage& globalData )

        : SparseMatrix<ValueType>( createStorage() )

    {
        DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
        DistributionPtr colDist( new NoDistribution( globalData.getNumRows() ) );

        SparseMatrix<ValueType>::assign( globalData, rowDist, colDist );
    }

    /** Constructor of a sparse matrix by local storage.
     *
     *  @param[in] localData   contains local rows of the distributed matrix
     *  @param[in] rowDist     is distribution of localData
     *  @param[in] colDist     specifies how to split local rows for halo
     *
     *  This constructor works also fine if localData is the full global matrix;
     *  in this case only local rows will be taken on this processor.
     */
    CSRSparseMatrix( const _MatrixStorage& localData, DistributionPtr rowDist, DistributionPtr colDist )

        : SparseMatrix<ValueType>( createStorage() )

    {
        SparseMatrix<ValueType>::assign( localData, rowDist, colDist );
    }

    /** Constructor of a replicated sparse matrix by reading the matrix
     *  data from a file.
     *
     *  @param[in] filename   name of the file where the matrix is read from
     *
     *  Next releases will also support distributed/parallel I/O. In the
     *  meantime this constructor should be used with a following call of
     *  the redistribute method.
     */
    explicit CSRSparseMatrix( const std::string& filename )

        : SparseMatrix<ValueType>( createStorage() )

    {
        this->readFromFile( filename );
    }

    // Expression constructors

    explicit CSRSparseMatrix( const Expression<Matrix,Matrix,Times>& expression )

        : SparseMatrix<ValueType>( createStorage() )
    {
        Matrix::operator=( expression );
    }

    explicit CSRSparseMatrix( const Expression<Scalar,Matrix,Times>& expression )

        : SparseMatrix<ValueType>( createStorage() )
    {
        Matrix::operator=( expression );
    }

    explicit CSRSparseMatrix( const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& expression )

        : SparseMatrix<ValueType>( createStorage() )
    {
        Matrix::operator=( expression );
    }

    /** @brief Constructor of CSRSparseMatrix by sum of two matrices.
     *
     *  @param expression is alpha * matA + beta * matB
     *
     */
    explicit CSRSparseMatrix(
        const Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus> expression )

        : SparseMatrix<ValueType>( createStorage() )
    {
        // inherit context from matA

        SparseMatrix<ValueType>::setContext( expression.getArg1().getArg2().getContextPtr() );
        Matrix::operator=( expression );
    }

    /** @brief Constructor of a CSR sparse matrix with distributed CSR storage data.
     *
     * @param[in] numLocalRows       the number of rows of the matrix
     * @param[in] numLocalNonZeros   the number of local none zeros of the matrix
     * @param[in] numHaloNonZeros    the number of halo none zeros of the matrix
     * @param[in] localIA            row pointer of the input csr sparse matrix (local)
     * @param[in] localJA            column indexes of the input csr sparse matrix (local)
     * @param[in] localValues        the none zero values of the input csr sparse matrix (local)
     * @param[in] haloIA             row pointer of the input csr sparse matrix (halo)
     * @param[in] haloJA             column indexes of the input csr sparse matrix (halo)
     * @param[in] haloValues         the none zero values of the input csr sparse matrix (halo)
     * @param[in] ownedIndexes       the global Indexes of the local rows
     * @param[in] communicator       communicator of the distribution
     */
    template<typename LocalValueType,typename HaloValueType>
    CSRSparseMatrix(
        const IndexType numLocalRows,
        const IndexType numLocalNonZeros,
        const IndexType numHaloNonZeros,
        const IndexType localIA[],
        const IndexType localJA[],
        const LocalValueType localValues[],
        const IndexType haloIA[],
        const IndexType haloJA[],
        const HaloValueType haloValues[],
        const std::vector<IndexType>& ownedIndexes,
        const CommunicatorPtr communicator )

        : SparseMatrix<ValueType>( createStorage() )

    {
        LAMA_LOG_INFO( logger,
                       communicator << ": construct distributed matrix " << numLocalRows << " by local and halo data + owned indexes" );

        // For the distribution we need the global number of rows, not available as arg, so compute it

        IndexType numGlobalRows = communicator->sum( numLocalRows );

        mLocalData->setRawCSRData( numLocalRows, numLocalRows, numLocalNonZeros, localIA, localJA, localValues );
        mHaloData->setRawCSRData( numLocalRows, numGlobalRows, numHaloNonZeros, haloIA, haloJA, haloValues );

        DistributionPtr dist( new GeneralDistribution( numGlobalRows, ownedIndexes, communicator ) );

        // Halo is already splitted, but still contains the global indexes

        mHaloData->buildHalo( mHalo, *dist ); // build halo, maps global indexes to halo indexes

        Matrix::setDistributedMatrix( dist, dist );
    }

    /**
     * @brief Destructor. Releases all allocated resources.
     */
    ~CSRSparseMatrix();

    /** Override the default assignment operator that would not make deep copies. */

    CSRSparseMatrix& operator=( const CSRSparseMatrix& matrix );

    /** Redefine assignment operator to get the correct return value; implementation is same as for base classes. */

    CSRSparseMatrix& operator=( const Matrix& matrix );

    CSRSparseMatrix& operator=( const Expression<Matrix,Matrix,Times>& expression );

    CSRSparseMatrix& operator=( const Expression<Scalar,Matrix,Times>& expression );

    CSRSparseMatrix& operator=( const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& expression );

    CSRSparseMatrix& operator=(
        const Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> exp );

    CSRSparseMatrix& operator=(
        const Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus> exp );

    /** Override MatrixStorage<ValueType>::getLocalStorage with covariant return type. */

    virtual const StorageType& getLocalStorage() const;

    /** @todo this getter should be removed as write access to local strage is dangerous */

    virtual StorageType& getLocalStorage();

    /** Override MatrixStorage<ValueType>::getHaloStorage with covariant return type. */

    virtual const StorageType& getHaloStorage() const;

    /** Swap local storage data, allows consistent write access to local storage. */

    virtual void swapLocalStorage( StorageType& localStorage );

    // TODO: create, copy with covariant return types

    /* Implementation of pure method of class Matrix. */

    virtual const char* getTypeName() const;

protected:

    using SparseMatrix<ValueType>::mLocalData;
    using SparseMatrix<ValueType>::mHaloData;
    using SparseMatrix<ValueType>::mHalo;

private:

    /** This private routine provides empty CSR storage for a CSRSparseMatrix. */

    boost::shared_ptr<MatrixStorage<ValueType> > createStorage()
    {
        return boost::shared_ptr<MatrixStorage<ValueType> >( new StorageType() );
    }

    boost::shared_ptr<MatrixStorage<ValueType> > createStorage( const IndexType numRows, const IndexType numColumns )
    {
        boost::shared_ptr<MatrixStorage<ValueType> > storage( new StorageType() );
        storage->allocate( numRows, numColumns );
        return storage;
    }

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace lama

#endif // LAMA_CSR_SPARSE_MATRIX_HPP_
