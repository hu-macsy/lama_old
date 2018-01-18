/**
 * @file lama/matrix/DIASparseMatrix.hpp
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
 * @brief Definition of matrix class for distributed sparse matrixes in DIA format.
 * @author Jiri Kraus, Thomas Brandes
 * @date 22.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/lama/matrix/SparseMatrix.hpp>

// local library
#include <scai/lama/storage/DIAStorage.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

namespace scai
{

namespace lama
{

/** Definition of a derived class for SparseMatrix that uses the DIA storage
 *  format for the local and halo data of the distributed sparse matrix.
 *
 *  As the storage format is known here this class can offer more advanced
 *  constructors that do not exist for SparseMatrix as there the storage
 *  format is not fixed.
 */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT DIASparseMatrix:

    public SparseMatrix<ValueType>,
    public _Matrix::Register<DIASparseMatrix<ValueType> >    // register at factory
{

public:

    /** Type definition of the storage type for this sparse matrix. */

    typedef DIAStorage<ValueType> StorageType;

    /** Static method that returns the name of the matrix class. */

    static const char* typeName();

    /** Default constructor, creates a replicated matrix of size 0 x 0 */

    DIASparseMatrix();

    /** Constructor, creates a replicated zero-matrix of size numRows x numColums */

    DIASparseMatrix( const IndexType numRows, const IndexType numColumns );

    /** Constructor, creates a distributed zero-matrix by given row and column distribution */

    DIASparseMatrix( dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Override default constructor, make sure that deep copies are created. */

    DIASparseMatrix( const DIASparseMatrix<ValueType>& other );

    /** Most general copy constrcuctor with possibility of transpose. */

    DIASparseMatrix( const Matrix<ValueType>& other, bool transposeFlag = false );

    /** Constructor of a sparse matrix by another input matrix with redistribution.
     *
     * @param[in] other     is the input matrix.
     * @param[in] rowDist   row distribution of the new matrix
     * @param[in] colDist   column distribution of the new matrix
     */
    DIASparseMatrix( const Matrix<ValueType>& other, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Constructor of a (replicated) sparse matrix by global storage.
     *
     *  @param[in] globalData  contains local rows of the distributed matrix
     */
    explicit DIASparseMatrix( const MatrixStorage<ValueType>& globalData );

    /** Constructor of a sparse matrix by local storage.
     *
     *  @param[in] localData   contains local rows of the distributed matrix
     *  @param[in] rowDist     is distribution of localData
     *  @param[in] colDist     specifies how to split local rows for halo
     *
     *  This constructor works also fine if localData is the full global matrix;
     *  in this case only local rows will be taken on this processor.
     */
    DIASparseMatrix( const MatrixStorage<ValueType>& localData, dmemo::DistributionPtr rowDist, dmemo::DistributionPtr colDist );

    /** Constructor of a replicated sparse matrix by reading the matrix
     *  data from a file.
     *
     *  @param[in] filename   name of the file where the matrix is read from
     *
     *  Next releases will also support distributed/parallel I/O. In the
     *  meantime this constructor should be used with a following call of
     *  the redistribute method.
     */
    explicit DIASparseMatrix( const std::string& filename );

    // Expression constructors

    explicit DIASparseMatrix( const Expression_SM<ValueType>& expression );

    explicit DIASparseMatrix( const Expression_SMM<ValueType>& expression );

    explicit DIASparseMatrix( const Expression_SM_SM<ValueType>& expression );

    template<typename OtherValueType>
    explicit DIASparseMatrix( const CastMatrixExpression<ValueType, OtherValueType>& expression );

    explicit DIASparseMatrix( const ComplexBuildMatrixExpression<RealType<ValueType> >& expression );

    template<common::ComplexSelection kind, typename OtherValueType>
    explicit DIASparseMatrix( const ComplexSelectionMatrixExpression<OtherValueType, kind>& expression );

    /** @brief Constructor of a DIA sparse matrix with distributed DIA storage data.
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
    template<typename LocalValueType, typename HaloValueType>
    DIASparseMatrix(
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
    ~DIASparseMatrix();

    // Make all assignment operators of base class visible before overwriting one

    using SparseMatrix<ValueType>::operator=;

    /** Override the default assignment operator that would not make deep copies. */

    DIASparseMatrix& operator=( const DIASparseMatrix& matrix );

    /** Override MatrixStorage<ValueType>::getLocalStorage with covariant return type. */

    virtual const StorageType& getLocalStorage() const;

    /** @todo this getter should be removed as write access to local strage is dangerous */

    virtual StorageType& getLocalStorage();

    /** Override MatrixStorage<ValueType>::getHaloStorage with covariant return type. */

    virtual const StorageType& getHaloStorage() const;

    /** Swap local storage data, allows consistent write access to local storage. */

    virtual void swapLocalStorage( StorageType& localStorage );

    /* Implementation of pure method _Matrix::newMatrix with covariant return type */

    virtual DIASparseMatrix<ValueType>* newMatrix() const;

    /* Implementation of pure method _Matrix::copy with covariant return type */

    virtual DIASparseMatrix<ValueType>* copy() const;

    /* Implementation of pure method _Matrix::getFormat */

    virtual Format getFormat() const
    {
        return Format::DIA;
    }

    /* Implementation of pure method of class _Matrix. */

    virtual const char* getTypeName() const;

    using SparseMatrix<ValueType>::setContextPtr;

protected:

    using SparseMatrix<ValueType>::mLocalData;
    using SparseMatrix<ValueType>::mHaloData;
    using SparseMatrix<ValueType>::mHalo;

private:

    /** This private routine provides empty DIA storage for a DIASparseMatrix. */

    std::shared_ptr<MatrixStorage<ValueType> > createStorage();

    std::shared_ptr<MatrixStorage<ValueType> > createStorage( const IndexType numRows, const IndexType numColumns );

    static std::string initTypeName();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

public:

    // static create method that will be used to register at _Matrix factory

    static _Matrix* create();

    // key for factory

    static MatrixCreateKeyType createValue();
};

template<typename ValueType>
template<common::ComplexSelection kind, typename OtherValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const ComplexSelectionMatrixExpression<OtherValueType, kind>& expression ) :

   SparseMatrix<ValueType>( createStorage() )

{
    const Matrix<OtherValueType>& m = expression.getArg();
    this->setContextPtr( m.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( m.getCommunicationKind() );
    Matrix<ValueType>::operator=( expression );
}

template<typename ValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const ComplexBuildMatrixExpression<RealType<ValueType> >& expression ) :

   SparseMatrix<ValueType>( createStorage() )

{
    const Matrix<RealType<ValueType> >& m = expression.getRealArg();
    this->setContextPtr( m.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( m.getCommunicationKind() );
    Matrix<ValueType>::operator=( expression );
}

template<typename ValueType>
template<typename OtherValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix( const CastMatrixExpression<ValueType, OtherValueType>& expression ) :

   SparseMatrix<ValueType>( createStorage() )

{
    const Matrix<OtherValueType>& m = expression.getArg();
    this->setContextPtr( m.getContextPtr() );
    SparseMatrix<ValueType>::setCommunicationKind( m.getCommunicationKind() );
    Matrix<ValueType>::operator=( expression );
}

template<typename ValueType>
template<typename LocalValueType, typename HaloValueType>
DIASparseMatrix<ValueType>::DIASparseMatrix(
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
    hmemo::HArrayRef<IndexType> myIndexes( static_cast<IndexType>( ownedIndexes.size() ), &ownedIndexes[0] );
    dmemo::DistributionPtr dist( new dmemo::GeneralDistribution( numGlobalRows, myIndexes, communicator ) );
    // Halo is already splitted, but still contains the global indexes
    mHaloData->buildHalo( mHalo, *dist ); // build halo, maps global indexes to halo indexes
    _Matrix::setDistributedMatrix( dist, dist );
}

} /* end namespace lama */

} /* end namespace scai */
