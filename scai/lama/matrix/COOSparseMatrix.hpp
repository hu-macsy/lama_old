/**
 * @file COOSparseMatrix.hpp
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
 * @brief Definition of matrix class for distributed sparse matrixes in COO format.
 * @author Jiri Kraus, Thomas Brandes
 * @date 22.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/SparseMatrix.hpp>

// local library
#include <scai/lama/storage/COOStorage.hpp>

#include <scai/lama/distribution/GeneralDistribution.hpp>
#include <scai/lama/distribution/NoDistribution.hpp>

namespace scai
{

namespace lama
{

/** Definition of a derived class for SparseMatrix that uses the COO storage
 *  format for the local and halo data of the distributed sparse matrix.
 *
 *  As the storage format is known here this class can offer more advanced
 *  constructors that do not exist for SparseMatrix as there the storage
 *  format is not fixed.
 */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT COOSparseMatrix: 

    public SparseMatrix<ValueType>,
    public Matrix::Register<COOSparseMatrix<ValueType> >    // register at factory
{

public:

    typedef ValueType MatrixValueType; //!< This is the type of the matrix values.

    /** Type definition of the storage type for this sparse matrix. */

    typedef COOStorage<ValueType> StorageType;

    /** Static method that returns the name of the matrix class. */

    static const char* typeName();

    /** Default constructor, creates a replicated matrix of size 0 x 0 */

    COOSparseMatrix();

    /** Constructor, creates a replicated zero-matrix of size numRows x numColums */

    COOSparseMatrix( const IndexType numRows, const IndexType numColumns );

    /** Constructor, creates a distributed zero-matrix by given row and column distribution */

    COOSparseMatrix( DistributionPtr rowDist, DistributionPtr colDist );

    /** Override default constructor, make sure that deep copies are created. */

    COOSparseMatrix( const COOSparseMatrix& other );

    /** Most general copy constrcuctor with possibility of transpose. */

    COOSparseMatrix( const Matrix& other, bool transposeFlag = false );

    /** Constructor of a sparse matrix by another input matrix with redistribution.
     *
     * @param[in] other     is the input matrix.
     * @param[in] rowDist   row distribution of the new matrix
     * @param[in] colDist   column distribution of the new matrix
     */
    COOSparseMatrix( const Matrix& other, DistributionPtr rowDist, DistributionPtr colDist );

    /** Constructor of a (replicated) sparse matrix by global storage.
     *
     *  @param[in] globalData  contains local rows of the distributed matrix
     */
    explicit COOSparseMatrix( const _MatrixStorage& globalData );

    /** Constructor of a sparse matrix by local storage.
     *
     *  @param[in] localData   contains local rows of the distributed matrix
     *  @param[in] rowDist     is distribution of localData
     *  @param[in] colDist     specifies how to split local rows for halo
     *
     *  This constructor works also fine if localData is the full global matrix;
     *  in this case only local rows will be taken on this processor.
     */
    COOSparseMatrix( const _MatrixStorage& localData, DistributionPtr rowDist, DistributionPtr colDist );

    /** Constructor of a replicated sparse matrix by reading the matrix
     *  data from a file.
     *
     *  @param[in] filename   name of the file where the matrix is read from
     *
     *  Next releases will also support distributed/parallel I/O. In the
     *  meantime this constructor should be used with a following call of
     *  the redistribute method.
     */
    explicit COOSparseMatrix( const std::string& filename );

    // Expression constructors

    explicit COOSparseMatrix( const Expression_SM& expression );

    explicit COOSparseMatrix( const Expression_SMM& expression );

    explicit COOSparseMatrix( const Expression_SM_SM& expression );

    /** @brief Constructor of a COO sparse matrix with distributed COO storage data.
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
    COOSparseMatrix(
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
        const CommunicatorPtr communicator );

    /**
     * @brief Destructor. Releases all allocated resources.
     */
    ~COOSparseMatrix();

    // Make all assignment operators of base class visible before overwriting one

    using SparseMatrix<ValueType>::operator=;

    /** Override the default assignment operator that would not make deep copies. */

    COOSparseMatrix& operator=( const COOSparseMatrix& matrix );

    /** Redefine assignment operator to get the correct return value; implementation is same as for base classes. */

    /*
     COOSparseMatrix& operator=( const Matrix& matrix );

     COOSparseMatrix& operator=( const Expression<Matrix,Matrix,Times>& expression );

     COOSparseMatrix& operator=( const Expression<Scalar,Matrix,Times>& expression );

     COOSparseMatrix& operator=( const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& expression );

     COOSparseMatrix& operator=(
     const Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> exp );

     COOSparseMatrix& operator=(
     const Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus> exp );
     */

    /** Override MatrixStorage<ValueType>::getLocalStorage with covariant return type. */

    virtual const StorageType& getLocalStorage() const;

    /** @todo this getter should be removed as write access to local strage is dangerous */

    virtual StorageType& getLocalStorage();

    /** Override MatrixStorage<ValueType>::getHaloStorage with covariant return type. */

    virtual const StorageType& getHaloStorage() const;

    /** Swap local storage data, allows consistent write access to local storage. */

    virtual void swapLocalStorage( StorageType& localStorage );

    /* Implementation of pure method Matrix::clone with covariant return type */

    virtual COOSparseMatrix<ValueType>* clone() const;

    /* Implementation of pure method Matrix::copy with covariant return type */

    virtual COOSparseMatrix<ValueType>* copy() const;

    /* Implementation of pure method of class Matrix. */

    virtual const char* getTypeName() const;

    using SparseMatrix<ValueType>::setContextPtr;

protected:

    using SparseMatrix<ValueType>::mLocalData;
    using SparseMatrix<ValueType>::mHaloData;
    using SparseMatrix<ValueType>::mHalo;

private:

    /** This private routine provides empty COO storage for a COOSparseMatrix. */

    common::shared_ptr<MatrixStorage<ValueType> > createStorage();

    common::shared_ptr<MatrixStorage<ValueType> > createStorage( const IndexType numRows, const IndexType numColumns );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

public:

    // static create method that will be used to register at Matrix factory

    static Matrix* create();

    // key for factory 

    static std::pair<MatrixStorageFormat, common::scalar::ScalarType> createValue();
};

template<typename ValueType>
template<typename LocalValueType,typename HaloValueType>
COOSparseMatrix<ValueType>::COOSparseMatrix(
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
    SCAI_LOG_INFO( logger,
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

} /* end namespace lama */

} /* end namespace scai */
