/**
 * @file lama/matrix/ELLSparseMatrix.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Definition of matrix class for distributed sparse matrixes in ELL format.
 * @author Jiri Kraus, Thomas Brandes
 * @date 22.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/SparseMatrix.hpp>

// local library
#include <scai/lama/storage/ELLStorage.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

namespace scai
{

namespace lama
{

/** Definition of a derived class for SparseMatrix that uses the ELL storage
 *  format for the local and halo data of the distributed sparse matrix.
 *
 *  As the storage format is known here this class can offer more advanced
 *  constructors that do not exist for SparseMatrix as there the storage
 *  format is not fixed.
 */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT ELLSparseMatrix: 

    public SparseMatrix<ValueType>,
    public Matrix::Register<ELLSparseMatrix<ValueType> >    // register at factory
{

public:

    typedef ValueType MatrixValueType; //!< This is the type of the matrix values.

    /** Type definition of the storage type for this sparse matrix. */

    typedef ELLStorage<ValueType> StorageType;

    /** Static method that returns the name of the matrix class. */

    static const char* typeName();

    /** Default constructor, creates a replicated matrix of size 0 x 0 */

    ELLSparseMatrix();

    /** Constructor, creates a replicated zero-matrix of size numRows x numColums */

    ELLSparseMatrix( const IndexType numRows, const IndexType numColumns );

    /** Constructor, creates a distributed zero-matrix by given row and column distribution */

    ELLSparseMatrix( dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Override default constructor, make sure that deep copies are created. */

    ELLSparseMatrix( const ELLSparseMatrix& other );

    /** Most general copy constrcuctor with possibility of transpose. */

    ELLSparseMatrix( const Matrix& other, bool transposeFlag = false );

    /** Constructor of a sparse matrix by another input matrix with redistribution.
     *
     * @param[in] other     is the input matrix.
     * @param[in] rowDist   row distribution of the new matrix
     * @param[in] colDist   column distribution of the new matrix
     */
    ELLSparseMatrix( const Matrix& other, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Constructor of a (replicated) sparse matrix by global storage.
     *
     *  @param[in] globalData  contains local rows of the distributed matrix
     */
    explicit ELLSparseMatrix( const _MatrixStorage& globalData );

    /** Constructor of a sparse matrix by local storage.
     *
     *  @param[in] localData   contains local rows of the distributed matrix
     *  @param[in] rowDist     is distribution of localData
     *  @param[in] colDist     specifies how to split local rows for halo
     *
     *  This constructor works also fine if localData is the full global matrix;
     *  in this case only local rows will be taken on this processor.
     */
    ELLSparseMatrix( const _MatrixStorage& localData, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Constructor of a replicated sparse matrix by reading the matrix
     *  data from a file.
     *
     *  @param[in] filename   name of the file where the matrix is read from
     *
     *  Next releases will also support distributed/parallel I/O. In the
     *  meantime this constructor should be used with a following call of
     *  the redistribute method.
     */
    explicit ELLSparseMatrix( const std::string& filename );

    // Expression constructors

    explicit ELLSparseMatrix( const Expression_SM& expression );

    explicit ELLSparseMatrix( const Expression_SMM& expression );

    explicit ELLSparseMatrix( const Expression_SM_SM& expression );

    /** @brief Constructor of a ELL sparse matrix with distributed ELL storage data.
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
    ELLSparseMatrix(
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
        const dmemo::CommunicatorPtr communicator );

    /**
     * @brief Destructor. Releases all allocated resources.
     */
    ~ELLSparseMatrix();

    // Make all assignment operators of base class visible before overwriting one

    using SparseMatrix<ValueType>::operator=;

    /** Override the default assignment operator that would not make deep copies. */

    ELLSparseMatrix& operator=( const ELLSparseMatrix& matrix );

    /** Redefine assignment operator to get the correct return value; implementation is same as for base classes. */

    /*
     ELLSparseMatrix& operator=( const Matrix& matrix );

     ELLSparseMatrix& operator=( const Expression<Matrix,Matrix,Times>& expression );

     ELLSparseMatrix& operator=( const Expression<Scalar,Matrix,Times>& expression );

     ELLSparseMatrix& operator=( const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& expression );

     ELLSparseMatrix& operator=(
     const Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> exp );

     ELLSparseMatrix& operator=(
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

    /* Implementation of pure method Matrix::newMatrix with covariant return type */

    virtual ELLSparseMatrix<ValueType>* newMatrix() const;

    /* Implementation of pure method Matrix::copy with covariant return type */

    virtual ELLSparseMatrix<ValueType>* copy() const;

    /* Implementation of pure method Matrix::getFormat */

    virtual Format::MatrixStorageFormat getFormat() const
    {
        return Format::ELL;
    }

    /* Implementation of pure method of class Matrix. */

    virtual const char* getTypeName() const;

    using SparseMatrix<ValueType>::setContextPtr;

protected:

    using SparseMatrix<ValueType>::mLocalData;
    using SparseMatrix<ValueType>::mHaloData;
    using SparseMatrix<ValueType>::mHalo;

private:

    /** This private routine provides empty ELL storage for a ELLSparseMatrix. */

    common::shared_ptr<MatrixStorage<ValueType> > createStorage();

    common::shared_ptr<MatrixStorage<ValueType> > createStorage( const IndexType numRows, const IndexType numColumns );

    static std::string initTypeName();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

public:

    // static create method that will be used to register at Matrix factory

    static Matrix* create();

    // key for factory 

    static MatrixCreateKeyType createValue();
};

template<typename ValueType>
template<typename LocalValueType,typename HaloValueType>
ELLSparseMatrix<ValueType>::ELLSparseMatrix(
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
    const dmemo::CommunicatorPtr communicator )

    : SparseMatrix<ValueType>( createStorage() )

{
    SCAI_LOG_INFO( logger,
                   communicator << ": construct distributed matrix " << numLocalRows << " by local and halo data + owned indexes" );

    // For the distribution we need the global number of rows, not available as arg, so compute it

    IndexType numGlobalRows = communicator->sum( numLocalRows );

    mLocalData->setRawCSRData( numLocalRows, numLocalRows, numLocalNonZeros, localIA, localJA, localValues );
    mHaloData->setRawCSRData( numLocalRows, numGlobalRows, numHaloNonZeros, haloIA, haloJA, haloValues );

    dmemo::DistributionPtr dist( new dmemo::GeneralDistribution( numGlobalRows, ownedIndexes, communicator ) );

    // Halo is already splitted, but still contains the global indexes

    mHaloData->buildHalo( mHalo, *dist ); // build halo, maps global indexes to halo indexes

    Matrix::setDistributedMatrix( dist, dist );
}

} /* end namespace lama */

} /* end namespace scai */
